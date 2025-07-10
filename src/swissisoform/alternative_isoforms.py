"""Alternative isoform data handling module.

This module provides functionality for loading and processing alternative
isoform data from BED format files, with support for reconstructing
biologically meaningful isoform tracks based on truncation patterns.
"""

import pandas as pd
from typing import Optional, Dict, List, Tuple, Set
from itertools import product


class AlternativeIsoform:
    """Handles alternative isoform data from BED format files.

    The name field in the BED file is expected to contain gene information in the format:
    ENSG00000260916.6_CCPG1_AUG_TruncationToAnno

    This class reconstructs biologically meaningful isoform tracks:
    - Truncations that overlap = same coding region
    - If 1 truncation per region -> merge all (same isoform)
    - If >1 truncation per region -> create combination tracks
    - Shorter truncations imply early termination
    - Strand determines transcription direction
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

    def _regions_overlap(
        self, region1: Tuple[int, int], region2: Tuple[int, int]
    ) -> bool:
        """Check if two genomic regions overlap by any amount.

        Args:
            region1: (start, end) of first region
            region2: (start, end) of second region

        Returns:
            True if regions overlap by at least 1bp
        """
        start1, end1 = region1
        start2, end2 = region2

        # True overlap (no gaps allowed)
        return not (end1 < start2 or end2 < start1)

    def _find_overlapping_groups(
        self, regions: List[Tuple[int, int, int]]
    ) -> List[List[int]]:
        """Find groups of overlapping regions using Union-Find.

        Args:
            regions: List of (start, end, index) tuples

        Returns:
            List of lists, where each inner list contains indices of overlapping regions
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
                region1 = (regions[i][0], regions[i][1])
                region2 = (regions[j][0], regions[j][1])
                if self._regions_overlap(region1, region2):
                    union(i, j)

        # Group indices by their root parent
        components = {}
        for i in range(n):
            root = find(i)
            if root not in components:
                components[root] = []
            components[root].append(i)

        return list(components.values())

    def _get_position_for_sorting(self, start: int, end: int, strand: str) -> int:
        """Get position for sorting truncations based on strand direction.

        Args:
            start: Start position
            end: End position
            strand: Strand (+ or -)

        Returns:
            Position to use for sorting
        """
        if strand == "+":
            return start  # Sort by start position for positive strand
        else:
            return -end  # Sort by end position (reversed) for negative strand

    def _should_terminate_early(
        self,
        truncation_idx: int,
        all_truncations: List[Dict],
        overlapping_groups: List[List[int]],
        strand: str,
    ) -> bool:
        """Determine if a truncation should cause early termination.

        Args:
            truncation_idx: Index of truncation to check
            all_truncations: List of all truncations for this gene
            overlapping_groups: Groups of overlapping truncations
            strand: Gene strand

        Returns:
            True if this truncation should terminate early
        """
        # Find which group this truncation belongs to
        truncation_group = None

        for group in overlapping_groups:
            if truncation_idx in group:
                truncation_group = group
                break

        if truncation_group is None or len(truncation_group) <= 1:
            return False  # No overlapping truncations, no early termination

        # Get truncations in this group sorted by position
        group_truncations = [all_truncations[i] for i in truncation_group]

        # Sort by genomic position based on strand
        if strand == "+":
            group_truncations.sort(key=lambda x: x["start"])
        else:
            group_truncations.sort(key=lambda x: -x["end"])

        # Find our truncation in the sorted list
        our_truncation = all_truncations[truncation_idx]

        # Check if there's a longer truncation in the same group
        for other_trunc in group_truncations:
            if other_trunc != our_truncation:
                our_length = our_truncation["end"] - our_truncation["start"]
                other_length = other_trunc["end"] - other_trunc["start"]
                if other_length > our_length:
                    return True  # There's a longer truncation, so this one terminates early

        return False

    def _generate_isoform_tracks(self, gene_truncations: pd.DataFrame) -> List[Dict]:
        """Generate all possible isoform tracks for a gene.

        Args:
            gene_truncations: DataFrame with truncations for a single gene

        Returns:
            List of isoform tracks, each containing truncations and metadata
        """
        if gene_truncations.empty:
            return []

        # Convert to list of dictionaries for easier manipulation
        truncations = gene_truncations.to_dict("records")

        # Add index to each truncation for tracking
        for i, trunc in enumerate(truncations):
            trunc["original_index"] = i

        # Get basic info
        strand = truncations[0]["strand"]
        gene_name = truncations[0]["gene_name"]

        # Create regions list for overlap detection
        regions = [
            (trunc["start"], trunc["end"], i) for i, trunc in enumerate(truncations)
        ]

        # Find overlapping groups
        overlapping_groups = self._find_overlapping_groups(regions)

        # Sort truncations by genomic position based on strand
        if strand == "+":
            truncations.sort(key=lambda x: x["start"])
        else:
            truncations.sort(key=lambda x: -x["end"])

        # Check if we have a simple case (1 truncation per group)
        if all(len(group) == 1 for group in overlapping_groups):
            # Simple case: merge all truncations into one track
            merged_start = min(trunc["start"] for trunc in truncations)
            merged_end = max(trunc["end"] for trunc in truncations)

            return [
                {
                    "truncations": truncations,
                    "start": merged_start,
                    "end": merged_end,
                    "gene_name": gene_name,
                    "strand": strand,
                    "track_id": f"{gene_name}_merged",
                    "description": f"Merged {len(truncations)} truncations",
                }
            ]

        # Complex case: multiple truncations in some groups
        # Generate all combinations
        isoform_tracks = []

        # Create groups with their truncations
        groups = []
        for group_indices in overlapping_groups:
            group_truncations = [truncations[i] for i in group_indices]
            groups.append(group_truncations)

        # Generate all combinations
        for combination in product(*groups):
            track_truncations = list(combination)

            # Check for early termination
            should_terminate = False
            termination_pos = len(track_truncations)

            for i, trunc in enumerate(track_truncations):
                if self._should_terminate_early(
                    trunc["original_index"], truncations, overlapping_groups, strand
                ):
                    should_terminate = True
                    termination_pos = i + 1
                    break

            # Truncate if needed
            if should_terminate:
                track_truncations = track_truncations[:termination_pos]

            if track_truncations:  # Only add non-empty tracks
                track_start = min(trunc["start"] for trunc in track_truncations)
                track_end = max(trunc["end"] for trunc in track_truncations)

                track_id = f"{gene_name}_track_{len(isoform_tracks) + 1}"
                description = f"Track with {len(track_truncations)} truncations"
                if should_terminate:
                    description += " (early termination)"

                isoform_tracks.append(
                    {
                        "truncations": track_truncations,
                        "start": track_start,
                        "end": track_end,
                        "gene_name": gene_name,
                        "strand": strand,
                        "track_id": track_id,
                        "description": description,
                    }
                )

        return isoform_tracks

    def get_visualization_features(self, gene_name: str) -> pd.DataFrame:
        """Get features formatted for visualization with biologically meaningful isoform tracks.

        Args:
            gene_name: Name of the gene to get features for

        Returns:
            DataFrame: Features formatted for visualization with separate isoform tracks

        Raises:
            ValueError: If no data has been loaded
        """
        if self.isoforms.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        # Get all truncations for this gene
        gene_truncations = self.isoforms[self.isoforms["gene_name"] == gene_name].copy()

        if gene_truncations.empty:
            return pd.DataFrame()

        # Generate isoform tracks
        isoform_tracks = self._generate_isoform_tracks(gene_truncations)

        if not isoform_tracks:
            return pd.DataFrame()

        # Convert tracks to visualization format
        viz_features = []

        for track in isoform_tracks:
            # Create one feature per track
            feature = {
                "chromosome": track["truncations"][0]["chrom"],
                "source": "truncation",
                "feature_type": "alternative_start",
                "start": track["start"],
                "end": track["end"],
                "score": track["truncations"][0]["score"],
                "strand": track["strand"],
                "frame": ".",
                "gene_id": track["truncations"][0]["gene_id"],
                "transcript_id": track["track_id"],
                "gene_name": track["gene_name"],
                "start_codon": track["truncations"][0]["start_codon"],
                "name": track["track_id"],
                "track_description": track["description"],
                "num_truncations": len(track["truncations"]),
            }
            viz_features.append(feature)

        return pd.DataFrame(viz_features)

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
            "total_truncations": len(self.isoforms),
            "unique_genes": len(self.get_gene_list()),
            "chromosomes": sorted(self.isoforms["chrom"].unique().tolist()),
            "start_codons": sorted(self.isoforms["start_codon"].unique().tolist()),
        }

        return stats

    def get_alternative_starts(self, gene_name: Optional[str] = None) -> pd.DataFrame:
        """Get alternative start sites as biologically meaningful isoform tracks.

        Args:
            gene_name: Gene name to filter by (if None, returns all genes)

        Returns:
            DataFrame: Alternative start sites with merged/separated tracks in original format

        Raises:
            ValueError: If no data has been loaded
        """
        if self.isoforms.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        if gene_name is not None:
            gene_truncations = self.isoforms[
                self.isoforms["gene_name"] == gene_name
            ].copy()
            tracks = self._generate_isoform_tracks(gene_truncations)
        else:
            # Process all genes
            tracks = []
            for gene in self.get_gene_list():
                gene_truncations = self.isoforms[
                    self.isoforms["gene_name"] == gene
                ].copy()
                gene_tracks = self._generate_isoform_tracks(gene_truncations)
                tracks.extend(gene_tracks)

        if not tracks:
            return pd.DataFrame()

        # Convert tracks back to the original DataFrame format
        result_rows = []

        for track in tracks:
            # Create one row per track (representing the merged/combined truncations)
            first_trunc = track["truncations"][0]  # Use first truncation as template

            # Create the row in original format
            row = {
                "chrom": first_trunc["chrom"],
                "start": track["start"],  # Use track start/end (merged coordinates)
                "end": track["end"],
                "name": track["track_id"],  # Use track ID as name
                "score": first_trunc["score"],
                "strand": track["strand"],
                "gene_id": first_trunc["gene_id"],
                "gene_name": track["gene_name"],
                "start_codon": first_trunc["start_codon"],
                "isoform_type": first_trunc["isoform_type"],
                "feature_type": "alternative_start",
                "source": "truncation",
                # Add track-specific metadata
                "track_description": track["description"],
                "num_truncations": len(track["truncations"]),
                "original_truncations": track[
                    "truncations"
                ],  # Keep reference to original data
            }
            result_rows.append(row)

        return pd.DataFrame(result_rows)

    def get_isoform_tracks(self, gene_name: Optional[str] = None) -> List[Dict]:
        """Get alternative start sites as biologically meaningful isoform track objects.

        Args:
            gene_name: Gene name to filter by (if None, returns all genes)

        Returns:
            List of isoform track dictionaries with detailed information

        Raises:
            ValueError: If no data has been loaded
        """
        if self.isoforms.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        if gene_name is not None:
            gene_truncations = self.isoforms[
                self.isoforms["gene_name"] == gene_name
            ].copy()
            return self._generate_isoform_tracks(gene_truncations)
        else:
            # Process all genes
            all_tracks = []
            for gene in self.get_gene_list():
                gene_truncations = self.isoforms[
                    self.isoforms["gene_name"] == gene
                ].copy()
                gene_tracks = self._generate_isoform_tracks(gene_truncations)
                all_tracks.extend(gene_tracks)
            return all_tracks

    def get_detailed_track_info(self, gene_name: str) -> Dict:
        """Get detailed information about isoform tracks for a gene.

        Args:
            gene_name: Name of the gene

        Returns:
            Dictionary with detailed track information
        """
        tracks = self.get_isoform_tracks(gene_name)

        if not tracks:
            return {"gene_name": gene_name, "tracks": [], "summary": "No tracks found"}

        track_info = []
        for i, track in enumerate(tracks):
            truncation_details = []
            for trunc in track["truncations"]:
                truncation_details.append(
                    {
                        "position": f"{trunc['start']}-{trunc['end']}",
                        "length": trunc["end"] - trunc["start"],
                        "start_codon": trunc["start_codon"],
                    }
                )

            track_info.append(
                {
                    "track_id": track["track_id"],
                    "description": track["description"],
                    "total_span": f"{track['start']}-{track['end']}",
                    "strand": track["strand"],
                    "truncations": truncation_details,
                }
            )

        return {
            "gene_name": gene_name,
            "num_tracks": len(tracks),
            "tracks": track_info,
            "summary": f"Generated {len(tracks)} biologically meaningful isoform tracks",
        }
