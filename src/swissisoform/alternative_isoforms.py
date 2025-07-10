"""Alternative isoform data handling module.

This module provides functionality for loading and processing alternative
isoform data from BED format files, with support for visualizing
alternative start sites and intelligently merging truncation regions.
"""

import pandas as pd
from typing import Optional, Dict, List, Tuple


class AlternativeIsoform:
    """Handles alternative isoform data from BED format files.

    The name field in the BED file is expected to contain gene information in the format:
    ENSG00000260916.6_CCPG1_AUG_TruncationToAnno
    
    This class intelligently merges truncation regions:
    - Same gene + start_codon + overlapping regions -> merged into one
    - Same gene + start_codon + non-overlapping regions -> kept separate
    - Different start_codons -> always kept separate
    """

    def __init__(self):
        """Initialize an empty isoform handler."""
        self.isoforms = pd.DataFrame()

    def load_bed(self, file_path: str) -> None:
        """Load data from BED format file.

        Args:
            file_path: Path to the BED file
        """
        # Read BED format with standard columns
        self.isoforms = pd.read_csv(
            file_path,
            sep="\t",
            names=["chrom", "start", "end", "name", "score", "strand"],
        )

        # Parse the name field which contains gene information
        def parse_name(name: str) -> Dict[str, str]:
            parts = name.split("_")
            if len(parts) >= 4:
                return {
                    "gene_id": parts[0],
                    "gene_name": parts[1],
                    "start_codon": parts[2],
                    "isoform_type": parts[3],
                }
            # Handle malformed names gracefully
            return {
                "gene_id": name,
                "gene_name": name,
                "start_codon": "UNK",
                "isoform_type": "UNK",
            }

        # Parse name field and add as new columns
        info_df = pd.DataFrame(self.isoforms["name"].apply(parse_name).tolist())
        self.isoforms = pd.concat([self.isoforms, info_df], axis=1)

        # Add columns needed for visualization
        self.isoforms["feature_type"] = "alternative_start"
        self.isoforms["source"] = "truncation"

        # Convert position columns to integers
        self.isoforms["start"] = self.isoforms["start"].astype(int)
        self.isoforms["end"] = self.isoforms["end"].astype(int)

    def _regions_overlap(self, region1: Tuple[int, int], region2: Tuple[int, int], max_gap: int = 0) -> bool:
        """Check if two genomic regions overlap or are within max_gap of each other.
        
        Args:
            region1: (start, end) of first region
            region2: (start, end) of second region
            max_gap: Maximum gap between regions to still consider them "overlapping"
            
        Returns:
            True if regions overlap or are within max_gap
        """
        start1, end1 = region1
        start2, end2 = region2
        
        # Check for overlap or adjacency within max_gap
        return not (end1 + max_gap < start2 or end2 + max_gap < start1)

    def _merge_overlapping_regions(self, regions: List[Tuple[int, int]], max_gap: int = 0) -> List[Tuple[int, int]]:
        """Merge overlapping or adjacent regions.
        
        Args:
            regions: List of (start, end) tuples
            max_gap: Maximum gap between regions to merge them
            
        Returns:
            List of merged (start, end) tuples
        """
        if not regions:
            return []
            
        # Sort regions by start position
        sorted_regions = sorted(regions)
        merged = []
        
        current_start, current_end = sorted_regions[0]
        
        for start, end in sorted_regions[1:]:
            if start <= current_end + max_gap + 1:  # Overlapping or adjacent
                current_end = max(current_end, end)
            else:
                # No overlap, save current and start new
                merged.append((current_start, current_end))
                current_start, current_end = start, end
        
        # Don't forget the last region
        merged.append((current_start, current_end))
        
        return merged

    def _find_connected_components(self, regions: List[Tuple[int, int]], max_gap: int = 0) -> List[List[int]]:
        """Find connected components of overlapping regions using Union-Find.
        
        Args:
            regions: List of (start, end) tuples
            max_gap: Maximum gap between regions to consider them connected
            
        Returns:
            List of lists, where each inner list contains indices of regions that should be merged
        """
        n = len(regions)
        if n <= 1:
            return [[i] for i in range(n)]
            
        # Union-Find data structure
        parent = list(range(n))
        
        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]
        
        def union(x, y):
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py
        
        # Check all pairs for overlap
        for i in range(n):
            for j in range(i + 1, n):
                if self._regions_overlap(regions[i], regions[j], max_gap):
                    union(i, j)
        
        # Group indices by their root parent
        components = {}
        for i in range(n):
            root = find(i)
            if root not in components:
                components[root] = []
            components[root].append(i)
        
        return list(components.values())

    def _merge_truncation_groups(self, features: pd.DataFrame, max_gap: int = 10000) -> pd.DataFrame:
        """Intelligently merge truncation regions based on overlap patterns.
        
        This function:
        1. Groups truncations by gene_name + start_codon + strand + chromosome
        2. For each group, finds connected components of overlapping regions
        3. Merges each connected component into a single truncation
        4. Keeps non-overlapping components as separate truncations
        
        Args:
            features: DataFrame with truncation features
            max_gap: Maximum gap between regions to consider them part of the same truncation
            
        Returns:
            DataFrame with intelligently merged truncation regions
        """
        if features.empty:
            return features
            
        merged_features = []
        
        # Group by the truncation identity (same biological truncation)
        group_columns = ['gene_name', 'start_codon', 'strand', 'chrom']
        
        for group_key, group_df in features.groupby(group_columns):
            gene_name, start_codon, strand, chrom = group_key
            
            if len(group_df) == 1:
                # Single region - no merging needed
                merged_features.append(group_df.iloc[0])
                continue
            
            # Extract regions as (start, end) tuples
            regions = [(row['start'], row['end']) for _, row in group_df.iterrows()]
            
            # Find connected components (groups of overlapping regions)
            components = self._find_connected_components(regions, max_gap)
            
            # Process each connected component
            for comp_idx, component_indices in enumerate(components):
                # Get the regions for this component
                component_regions = [regions[i] for i in component_indices]
                
                # Merge overlapping regions within this component
                merged_regions = self._merge_overlapping_regions(component_regions, max_gap)
                
                # Create merged features for each merged region in this component
                for merged_start, merged_end in merged_regions:
                    # Use the first region in the component as template
                    template_idx = component_indices[0]
                    template_row = group_df.iloc[template_idx]
                    
                    # Create merged feature
                    merged_feature = template_row.copy()
                    merged_feature['start'] = merged_start
                    merged_feature['end'] = merged_end
                    
                    # Create descriptive name
                    num_original = len(component_indices)
                    merged_feature['name'] = f"{gene_name}_{start_codon}_merged_{num_original}regions_{merged_start}_{merged_end}"
                    
                    merged_features.append(merged_feature)
        
        if merged_features:
            result = pd.DataFrame(merged_features)
            return result
        else:
            return features

    def get_visualization_features(self, gene_name: str, max_gap: int = 10000) -> pd.DataFrame:
        """Get features formatted for visualization with intelligent truncation merging.

        Args:
            gene_name: Name of the gene to get features for
            max_gap: Maximum gap between truncation regions to merge them (default: 10kb)

        Returns:
            DataFrame: Features formatted for visualization with merged truncations

        Raises:
            ValueError: If no data has been loaded
        """
        if self.isoforms.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        # Get all features for this gene
        features = self.isoforms[self.isoforms["gene_name"] == gene_name].copy()

        if features.empty:
            return pd.DataFrame()

        # Intelligently merge truncation regions
        merged_features = self._merge_truncation_groups(features, max_gap)

        # Format features for visualization
        viz_features = pd.DataFrame(
            {
                "chromosome": merged_features["chrom"],
                "source": "truncation",
                "feature_type": "alternative_start",
                "start": merged_features["start"],
                "end": merged_features["end"],
                "score": merged_features["score"],
                "strand": merged_features["strand"],
                "frame": ".",
                "gene_id": merged_features["gene_id"],
                "transcript_id": merged_features["gene_id"] + "_alt",
                "gene_name": merged_features["gene_name"],
                "start_codon": merged_features["start_codon"],
                "name": merged_features["name"],
            }
        )

        return viz_features

    def get_gene_list(self) -> List[str]:
        """Get list of all genes in the dataset.

        Returns:
            List of gene names

        Raises:
            ValueError: If no data has been loaded
        """
        if self.isoforms.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        return sorted(self.isoforms["gene_name"].unique().tolist())

    def get_stats(self) -> Dict:
        """Get basic statistics about the loaded isoforms.

        Returns:
            Dictionary containing various statistics

        Raises:
            ValueError: If no data has been loaded
        """
        if self.isoforms.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        stats = {
            "total_sites": len(self.isoforms),
            "unique_genes": len(self.get_gene_list()),
            "chromosomes": sorted(self.isoforms["chrom"].unique().tolist()),
            "start_codons": sorted(self.isoforms["start_codon"].unique().tolist()),
        }

        return stats

    def get_alternative_starts(self, gene_name: Optional[str] = None, max_gap: int = 10000) -> pd.DataFrame:
        """Get alternative start sites with intelligent merging, optionally filtered by gene name.

        Args:
            gene_name: Gene name to filter by
            max_gap: Maximum gap between truncation regions to merge them

        Returns:
            DataFrame: Alternative start sites with positions and metadata

        Raises:
            ValueError: If no data has been loaded
        """
        if self.isoforms.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        if gene_name is not None:
            features = self.isoforms[self.isoforms["gene_name"] == gene_name].copy()
            return self._merge_truncation_groups(features, max_gap)
        else:
            # Apply merging to all genes
            all_merged = []
            for gene in self.get_gene_list():
                gene_features = self.isoforms[self.isoforms["gene_name"] == gene].copy()
                merged_gene_features = self._merge_truncation_groups(gene_features, max_gap)
                if not merged_gene_features.empty:
                    all_merged.append(merged_gene_features)
            
            if all_merged:
                return pd.concat(all_merged, ignore_index=True)
            else:
                return pd.DataFrame()