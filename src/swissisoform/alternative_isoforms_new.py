"""Alternative start sites data handling module.

This module provides functionality for loading and processing alternative
start sites from BED format files, calculating truncation/extension regions,
and supporting efficiency-based filtering.

Expects BED name fields in format: GENE_ENSEMBL_TYPE_CODON_EFFICIENCY
"""

import pandas as pd
import numpy as np
from typing import Optional, Dict, List, Tuple, Set
from itertools import product


class AlternativeIsoform:
    """Handles alternative start sites data from BED format files.

    The name field in the BED file is expected to contain gene information in the format:
    GENE_ENSEMBL_TYPE_CODON_EFFICIENCY (e.g., TIMM10B_ENSG00000132286.11_Annotated_AUG_15.15592)

    This class processes start sites to:
    - Group start sites by gene and type (Annotated, Truncated, Extended)
    - Calculate truncation/extension regions relative to canonical starts
    - Support efficiency-based filtering for top-performing start sites
    - Generate features for visualization and downstream analysis
    """

    def __init__(self, debug=False):
        """Initialize an empty start sites handler."""
        self.start_sites = pd.DataFrame()
        self.debug = debug

    def _debug_print(self, message: str):
        """Print debug message if debug mode is enabled."""
        if self.debug:
            print(f"DEBUG: {message}")

    def load_bed(self, file_path: str) -> None:
        """Load start sites data from BED format file.

        Args:
            file_path: Path to the BED file
        """
        # Read BED format with standard columns
        self.start_sites = pd.read_csv(
            file_path,
            sep="\t",
            names=["chrom", "start", "end", "name", "score", "strand"],
        )

        # Parse the name field which contains gene information in format:
        # GENE_ENSEMBL_TYPE_CODON_EFFICIENCY (e.g., TIMM10B_ENSG00000132286.11_Annotated_AUG_15.15592)
        def parse_name(name: str) -> Dict[str, str]:
            parts = name.split("_")

            # Expected format: TIMM10B_ENSG00000132286.11_Annotated_AUG_15.15592
            if len(parts) >= 5 and parts[1].startswith("ENSG"):
                return {
                    "gene_name": parts[0],
                    "gene_id": parts[1],
                    "start_type": parts[2],  # Annotated, Truncated, Extended
                    "start_codon": parts[3],
                    "efficiency": float(parts[4])
                    if parts[4].replace(".", "").isdigit()
                    else 0.0,
                }
            # Handle malformed names gracefully
            else:
                return {
                    "gene_id": name,
                    "gene_name": name,
                    "start_type": "Unknown",
                    "start_codon": "UNK",
                    "efficiency": 0.0,
                }

        # Parse name field and add as new columns
        info_df = pd.DataFrame(self.start_sites["name"].apply(parse_name).tolist())
        self.start_sites = pd.concat([self.start_sites, info_df], axis=1)

        # Add columns needed for visualization
        self.start_sites["feature_type"] = "alternative_start"
        self.start_sites["source"] = "ribosome_profiling"

        # Convert position columns to integers
        self.start_sites["start"] = self.start_sites["start"].astype(int)
        self.start_sites["end"] = self.start_sites["end"].astype(int)

        # Calculate start sites
        self.start_sites["start_position"] = self.start_sites["start"]

        self._debug_print(f"Loaded {len(self.start_sites)} start sites")
        self._debug_print(
            f"Start types: {self.start_sites['start_type'].value_counts().to_dict()}"
        )

    def get_gene_start_sites(self, gene_name: str) -> pd.DataFrame:
        """Get all start sites for a specific gene.

        Args:
            gene_name: Name of the gene

        Returns:
            DataFrame of start sites for the gene
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        gene_sites = self.start_sites[self.start_sites["gene_name"] == gene_name].copy()
        return gene_sites.sort_values("start_position")

    def calculate_truncation_regions(
        self, gene_name: str, top_n_efficiency: Optional[int] = None
    ) -> pd.DataFrame:
        """Calculate truncation/extension regions for a gene.

        Args:
            gene_name: Name of the gene
            top_n_efficiency: If specified, only use top N most efficient alternative starts

        Returns:
            DataFrame containing truncation regions with metadata
        """
        gene_sites = self.get_gene_start_sites(gene_name)

        if gene_sites.empty:
            return pd.DataFrame()

        # Find canonical (Annotated) start site
        canonical_sites = gene_sites[gene_sites["start_type"] == "Annotated"]
        if canonical_sites.empty:
            self._debug_print(f"Warning: No Annotated start site found for {gene_name}")
            return pd.DataFrame()

        # Use the most efficient canonical start if multiple exist
        canonical_start = canonical_sites.loc[canonical_sites["efficiency"].idxmax()]
        strand = canonical_start["strand"]
        canonical_pos = canonical_start["start"]

        self._debug_print(
            f"Gene {gene_name}: canonical start at {canonical_pos} (strand {strand})"
        )

        # Get alternative start sites (Truncated and Extended)
        alternative_sites = gene_sites[
            gene_sites["start_type"].isin(["Truncated", "Extended"])
        ]

        if alternative_sites.empty:
            self._debug_print(f"No alternative start sites found for {gene_name}")
            return pd.DataFrame()

        # Apply efficiency filtering if requested
        if top_n_efficiency is not None and len(alternative_sites) > top_n_efficiency:
            alternative_sites = alternative_sites.nlargest(
                top_n_efficiency, "efficiency"
            )
            self._debug_print(
                f"Filtered to top {top_n_efficiency} most efficient alternative starts"
            )

        # Calculate truncation regions
        truncation_regions = []

        for _, alt_site in alternative_sites.iterrows():
            alt_pos = alt_site["start"]

            # Then update the region calculations:
            if strand == "+":
                if alt_pos > canonical_pos:
                    # Downstream alternative start = N-terminal truncation
                    region_start = canonical_pos
                    region_end = alt_pos - 1  # BASE PAIR PRECEDING the truncated start
                    region_type = "truncation"
                    region_length = (alt_pos - 1) - canonical_pos + 1
                else:
                    # Upstream alternative start = N-terminal extension
                    region_start = alt_pos
                    region_end = (
                        canonical_pos - 1
                    )  # BASE PAIR PRECEDING canonical start
                    region_type = "extension"
                    region_length = (canonical_pos - 1) - alt_pos + 1
            else:  # strand == "-"
                if alt_pos < canonical_pos:
                    # Upstream alternative start (on - strand) = N-terminal truncation
                    region_start = alt_pos
                    region_end = (
                        canonical_pos - 1
                    )  # BASE PAIR PRECEDING canonical start
                    region_type = "truncation"
                    region_length = (canonical_pos - 1) - alt_pos + 1
                else:
                    # Downstream alternative start (on - strand) = N-terminal extension
                    region_start = canonical_pos
                    region_end = alt_pos - 1  # BASE PAIR PRECEDING alternative start
                    region_type = "extension"
                    region_length = (alt_pos - 1) - canonical_pos + 1

            truncation_regions.append(
                {
                    "chromosome": alt_site["chrom"],
                    "region_start": region_start,
                    "region_end": region_end,
                    "region_type": region_type,
                    "region_length": region_length,
                    "strand": strand,
                    "gene_name": gene_name,
                    "gene_id": alt_site["gene_id"],
                    "canonical_start_pos": canonical_pos,
                    "alternative_start_pos": alt_pos,
                    "alternative_start_type": alt_site["start_type"],
                    "alternative_start_codon": alt_site["start_codon"],
                    "efficiency": alt_site["efficiency"],
                    "canonical_efficiency": canonical_start["efficiency"],
                    # Create unique identifier for this truncation
                    "truncation_id": f"{gene_name}_{region_type}_{alt_site['start_codon']}_{alt_pos}",
                    "track_id": f"{gene_name}_{region_type}_{alt_pos}",
                }
            )

        return pd.DataFrame(truncation_regions)

    def get_visualization_features(
        self,
        gene_name: str,
        top_n_efficiency: Optional[int] = None,
        top_n_per_type: Optional[int] = None,
    ) -> pd.DataFrame:
        """Get features formatted for visualization.

        Args:
            gene_name: Name of the gene to get features for
            top_n_efficiency: If specified, only use top N most efficient alternative starts (excludes Annotated)
            top_n_per_type: If specified, use top N per start type (Annotated, Truncated, Extended independently)

        Returns:
            DataFrame: Features formatted for visualization
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        truncation_regions = self.calculate_truncation_regions(
            gene_name, top_n_efficiency
        )

        if truncation_regions.empty:
            return pd.DataFrame()

        # Convert to visualization format
        viz_features = []

        for _, region in truncation_regions.iterrows():
            feature = {
                "chromosome": region["chromosome"],
                "source": "ribosome_profiling",
                "feature_type": "alternative_start_region",
                "start": region["region_start"],
                "end": region["region_end"],
                "score": region["efficiency"],
                "strand": region["strand"],
                "frame": ".",
                "gene_id": region["gene_id"],
                "transcript_id": region["track_id"],
                "gene_name": region["gene_name"],
                "start_codon": region["alternative_start_codon"],
                "name": region["truncation_id"],
                "region_type": region["region_type"],
                "region_length": region["region_length"],
                "efficiency": region["efficiency"],
                "canonical_efficiency": region["canonical_efficiency"],
                "efficiency_ratio": region["efficiency"]
                / region["canonical_efficiency"]
                if region["canonical_efficiency"] > 0
                else 0,
            }
            viz_features.append(feature)

        return pd.DataFrame(viz_features)

    def get_gene_list(self) -> List[str]:
        """Get list of all genes in the dataset.

        Returns:
            List of gene names
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        return sorted(self.start_sites["gene_name"].unique().tolist())

    def get_stats(self) -> Dict:
        """Get basic statistics about the loaded start sites.

        Returns:
            Dictionary containing various statistics
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        stats = {
            "total_start_sites": len(self.start_sites),
            "unique_genes": len(self.get_gene_list()),
            "chromosomes": sorted(self.start_sites["chrom"].unique().tolist()),
            "start_types": self.start_sites["start_type"].value_counts().to_dict(),
            "start_codons": sorted(self.start_sites["start_codon"].unique().tolist()),
            "efficiency_stats": {
                "mean": self.start_sites["efficiency"].mean(),
                "median": self.start_sites["efficiency"].median(),
                "min": self.start_sites["efficiency"].min(),
                "max": self.start_sites["efficiency"].max(),
            },
        }

        return stats

    def filter_by_efficiency(
        self,
        min_efficiency: Optional[float] = None,
        top_n_per_gene: Optional[int] = None,
    ) -> "AlternativeIsoform":
        """Create a filtered copy based on efficiency criteria.

        Args:
            min_efficiency: Minimum efficiency threshold
            top_n_per_gene: Keep only top N most efficient alternative starts per gene

        Returns:
            New AlternativeIsoform instance with filtered data
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        filtered_data = self.start_sites.copy()

        # Apply minimum efficiency filter
        if min_efficiency is not None:
            filtered_data = filtered_data[filtered_data["efficiency"] >= min_efficiency]
            self._debug_print(
                f"Filtered by min efficiency {min_efficiency}: {len(filtered_data)} sites remaining"
            )

        # Apply top N per gene filter
        if top_n_per_gene is not None:
            # Keep all Annotated sites, filter only alternatives
            annotated_sites = filtered_data[filtered_data["start_type"] == "Annotated"]
            alternative_sites = filtered_data[
                filtered_data["start_type"] != "Annotated"
            ]

            # Get top N alternatives per gene
            top_alternatives = (
                alternative_sites.groupby("gene_name")
                .apply(lambda x: x.nlargest(top_n_per_gene, "efficiency"))
                .reset_index(drop=True)
            )

            filtered_data = pd.concat(
                [annotated_sites, top_alternatives], ignore_index=True
            )
            self._debug_print(
                f"Filtered to top {top_n_per_gene} alternatives per gene: {len(filtered_data)} sites remaining"
            )

        # Create new instance with filtered data
        filtered_instance = AlternativeIsoform(debug=self.debug)
        filtered_instance.start_sites = filtered_data

        return filtered_instance

    def get_gene_summary(self, gene_name: str) -> Dict:
        """Get detailed summary for a specific gene.

        Args:
            gene_name: Name of the gene

        Returns:
            Dictionary with gene summary information
        """
        gene_sites = self.get_gene_start_sites(gene_name)

        if gene_sites.empty:
            return {"gene_name": gene_name, "error": "No start sites found"}

        # Count by type
        type_counts = gene_sites["start_type"].value_counts().to_dict()

        # Get canonical start info
        canonical_sites = gene_sites[gene_sites["start_type"] == "Annotated"]
        canonical_info = None
        if not canonical_sites.empty:
            best_canonical = canonical_sites.loc[canonical_sites["efficiency"].idxmax()]
            canonical_info = {
                "position": best_canonical["start_position"],
                "efficiency": best_canonical["efficiency"],
                "codon": best_canonical["start_codon"],
            }

        # Get alternative start stats
        alternatives = gene_sites[gene_sites["start_type"] != "Annotated"]
        alt_stats = {}
        if not alternatives.empty:
            alt_stats = {
                "count": len(alternatives),
                "efficiency_range": [
                    alternatives["efficiency"].min(),
                    alternatives["efficiency"].max(),
                ],
                "most_efficient": {
                    "position": alternatives.loc[
                        alternatives["efficiency"].idxmax(), "start_position"
                    ],
                    "efficiency": alternatives["efficiency"].max(),
                    "type": alternatives.loc[
                        alternatives["efficiency"].idxmax(), "start_type"
                    ],
                },
            }

        return {
            "gene_name": gene_name,
            "total_start_sites": len(gene_sites),
            "start_type_counts": type_counts,
            "canonical_start": canonical_info,
            "alternative_starts": alt_stats,
            "chromosome": gene_sites.iloc[0]["chrom"],
            "strand": gene_sites.iloc[0]["strand"],
        }
