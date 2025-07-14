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

    def __init__(self, debug=False):
        """Initialize an empty isoform handler."""
        self.isoforms = pd.DataFrame()
        self.debug = debug

    def _debug_print(self, message: str):
        """Print debug message if debug mode is enabled."""
        if self.debug:
            print(f"DEBUG: {message}")

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

        self._debug_print(f"Checking overlaps for {n} regions:")
        for i, (start, end, idx) in enumerate(regions):
            self._debug_print(f"  Region {i}: {start}-{end} (original index {idx})")

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

                overlap = self._regions_overlap(region1, region2)
                self._debug_print(
                    f"Overlap check: Region {i} {region1} vs Region {j} {region2} = {overlap}"
                )

                if overlap:
                    self._debug_print(f"  -> Unioning regions {i} and {j}")
                    union(i, j)

        # Group indices by their root parent
        components = {}
        for i in range(n):
            root = find(i)
            if root not in components:
                components[root] = []
            components[root].append(i)

        result = list(components.values())
        self._debug_print(f"Final overlapping groups: {result}")
        return result

    def _generate_isoform_tracks(self, gene_truncations: pd.DataFrame) -> List[Dict]:
        """Generate all possible isoform tracks for a gene using improved logic.

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

        # Remove exact duplicates based on coordinates
        unique_truncations = []
        seen_coordinates = set()

        for trunc in truncations:
            coord_key = (trunc["start"], trunc["end"])
            if coord_key not in seen_coordinates:
                unique_truncations.append(trunc)
                seen_coordinates.add(coord_key)
            else:
                self._debug_print(
                    f"Removing duplicate truncation: {trunc['start']}-{trunc['end']}"
                )

        if len(unique_truncations) != len(truncations):
            self._debug_print(
                f"Removed {len(truncations) - len(unique_truncations)} duplicate truncations"
            )
            truncations = unique_truncations
            # Reassign indices after deduplication
            for i, trunc in enumerate(truncations):
                trunc["original_index"] = i

        # Get basic info
        strand = truncations[0]["strand"]
        gene_name = truncations[0]["gene_name"]

        self._debug_print(f"\n=== TRACK GENERATION FOR {gene_name} ===")
        self._debug_print(f"Gene: {gene_name}, Strand: {strand}")
        self._debug_print(
            f"Input truncations: {[(t['start'], t['end'], t['end'] - t['start']) for t in truncations]}"
        )

        # Step 1: Find the rightmost and leftmost positions
        all_starts = [t["start"] for t in truncations]
        all_ends = [t["end"] for t in truncations]
        leftmost_pos = min(all_starts)
        rightmost_pos = max(all_ends)

        self._debug_print(f"Overall span: {leftmost_pos} to {rightmost_pos}")

        # Create regions list for overlap detection
        regions = [
            (trunc["start"], trunc["end"], i) for i, trunc in enumerate(truncations)
        ]

        # Find overlapping groups
        overlapping_groups = self._find_overlapping_groups(regions)

        self._debug_print(f"Overlapping groups: {overlapping_groups}")

        # Step 2: Check if there are overlapping truncations
        has_overlaps = any(len(group) > 1 for group in overlapping_groups)

        if not has_overlaps:
            self._debug_print("SIMPLE CASE: No overlapping truncations found")
            # Sort truncations by genomic position based on strand
            if strand == "+":
                truncations.sort(key=lambda x: x["start"])
            else:
                truncations.sort(key=lambda x: -x["end"])

            # Create single merged track
            track = {
                "truncations": truncations,
                "start": leftmost_pos,
                "end": rightmost_pos,
                "gene_name": gene_name,
                "strand": strand,
                "track_id": f"{gene_name}_merged",
                "description": f"Merged {len(truncations)} non-overlapping truncations",
            }

            self._debug_print(
                f"Created single merged track: {leftmost_pos}-{rightmost_pos}"
            )
            return [track]

        # Step 3: Process overlapping truncations
        self._debug_print("COMPLEX CASE: Overlapping truncations found")

        # Process each overlapping group
        all_tracks = []

        for group_idx, group_indices in enumerate(overlapping_groups):
            if len(group_indices) <= 1:
                continue  # Skip single-truncation groups for now

            self._debug_print(
                f"\n--- Processing overlapping group {group_idx}: {group_indices} ---"
            )

            # Get truncations in this overlapping group
            group_truncations = []
            for orig_idx in group_indices:
                for trunc in truncations:
                    if trunc["original_index"] == orig_idx:
                        group_truncations.append(trunc)
                        break

            self._debug_print(
                f"Group truncations: {[(t['start'], t['end'], t['end'] - t['start']) for t in group_truncations]}"
            )

            # Check if truncations diverge on the correct side
            track_variants = self._process_overlapping_group(
                group_truncations, truncations, strand, gene_name
            )
            all_tracks.extend(track_variants)

        # If we didn't generate any tracks from overlapping groups, fall back to simple case
        if not all_tracks:
            self._debug_print(
                "No valid tracks from overlapping groups, falling back to simple case"
            )
            if strand == "+":
                truncations.sort(key=lambda x: x["start"])
            else:
                truncations.sort(key=lambda x: -x["end"])

            track = {
                "truncations": truncations,
                "start": leftmost_pos,
                "end": rightmost_pos,
                "gene_name": gene_name,
                "strand": strand,
                "track_id": f"{gene_name}_merged",
                "description": f"Merged {len(truncations)} truncations (fallback)",
            }
            all_tracks = [track]

        self._debug_print(f"\nFinal result: {len(all_tracks)} tracks generated")
        return all_tracks

    def _process_overlapping_group(
        self,
        group_truncations: List[Dict],
        all_truncations: List[Dict],
        strand: str,
        gene_name: str,
    ) -> List[Dict]:
        """Process a group of overlapping truncations according to the new logic.

        Args:
            group_truncations: List of overlapping truncations
            all_truncations: List of all truncations for the gene
            strand: Gene strand
            gene_name: Gene name

        Returns:
            List of track variants for this overlapping group
        """
        self._debug_print(
            f"Processing {len(group_truncations)} overlapping truncations"
        )

        # Step 3a: Check if truncations diverge on the correct side
        group_starts = [t["start"] for t in group_truncations]
        group_ends = [t["end"] for t in group_truncations]

        unique_starts = len(set(group_starts))
        unique_ends = len(set(group_ends))

        self._debug_print(
            f"Unique start positions: {unique_starts}, Unique end positions: {unique_ends}"
        )

        if strand == "-":
            # For negative strand, check if they diverge on the left (start positions)
            if unique_starts <= 1:
                self._debug_print(
                    "ERROR: Overlapping truncations don't diverge on left side (negative strand)"
                )
                self._debug_print(
                    "This appears to be an error in the BED file - treating as single truncation"
                )
                # Keep the longest truncation
                longest = max(group_truncations, key=lambda x: x["end"] - x["start"])
                return self._create_simple_track(
                    [longest], all_truncations, strand, gene_name, "error_corrected"
                )

            # Check if one truncation is completely contained within another
            for i, trunc_i in enumerate(group_truncations):
                for j, trunc_j in enumerate(group_truncations):
                    if i != j:
                        # Check if trunc_j is completely contained within trunc_i
                        if (
                            trunc_i["start"] <= trunc_j["start"]
                            and trunc_j["end"] <= trunc_i["end"]
                        ):
                            self._debug_print(
                                f"ERROR: Truncation {trunc_j['start']}-{trunc_j['end']} is completely contained within {trunc_i['start']}-{trunc_i['end']}"
                            )
                            self._debug_print(
                                "This appears to be an error in the BED file - treating as single truncation"
                            )
                            # Keep the longest (containing) truncation
                            longest = max(
                                group_truncations, key=lambda x: x["end"] - x["start"]
                            )
                            return self._create_simple_track(
                                [longest],
                                all_truncations,
                                strand,
                                gene_name,
                                "containment_error_corrected",
                            )

            self._debug_print("Valid divergence on left side for negative strand")
            divergent_positions = sorted(
                set(group_starts), reverse=True
            )  # Sort high to low for neg strand

        else:  # strand == "+"
            # For positive strand, check if they diverge on the right (end positions)
            if unique_ends <= 1:
                self._debug_print(
                    "ERROR: Overlapping truncations don't diverge on right side (positive strand)"
                )
                self._debug_print(
                    "This appears to be an error in the BED file - treating as single truncation"
                )
                # Keep the longest truncation
                longest = max(group_truncations, key=lambda x: x["end"] - x["start"])
                return self._create_simple_track(
                    [longest], all_truncations, strand, gene_name, "error_corrected"
                )

            # Check if one truncation is completely contained within another
            for i, trunc_i in enumerate(group_truncations):
                for j, trunc_j in enumerate(group_truncations):
                    if i != j:
                        # Check if trunc_j is completely contained within trunc_i
                        if (
                            trunc_i["start"] <= trunc_j["start"]
                            and trunc_j["end"] <= trunc_i["end"]
                        ):
                            self._debug_print(
                                f"ERROR: Truncation {trunc_j['start']}-{trunc_j['end']} is completely contained within {trunc_i['start']}-{trunc_i['end']}"
                            )
                            self._debug_print(
                                "This appears to be an error in the BED file - treating as single truncation"
                            )
                            # Keep the longest (containing) truncation
                            longest = max(
                                group_truncations, key=lambda x: x["end"] - x["start"]
                            )
                            return self._create_simple_track(
                                [longest],
                                all_truncations,
                                strand,
                                gene_name,
                                "containment_error_corrected",
                            )

            self._debug_print("Valid divergence on right side for positive strand")
            divergent_positions = sorted(
                set(group_ends)
            )  # Sort low to high for pos strand

        # Step 3b: Check if divergence happens at the extreme end
        all_starts = [t["start"] for t in all_truncations]
        all_ends = [t["end"] for t in all_truncations]

        if strand == "-":
            # Check if divergence is at extreme left (highest start position)
            extreme_position = max(all_starts)
            divergence_at_extreme = max(group_starts) == extreme_position
            self._debug_print(
                f"Extreme left position: {extreme_position}, Group max start: {max(group_starts)}, At extreme: {divergence_at_extreme}"
            )

        else:  # strand == "+"
            # Check if divergence is at extreme right (highest end position)
            extreme_position = max(all_ends)
            divergence_at_extreme = max(group_ends) == extreme_position
            self._debug_print(
                f"Extreme right position: {extreme_position}, Group max end: {max(group_ends)}, At extreme: {divergence_at_extreme}"
            )

        if divergence_at_extreme:
            self._debug_print(
                "CASE: Divergence at extreme end - creating tracks with different end sites only"
            )
            return self._create_different_endpoint_tracks(
                group_truncations, all_truncations, strand, gene_name
            )
        else:
            self._debug_print(
                "CASE: Divergence at intermediate position - creating tracks with early termination"
            )
            return self._create_early_termination_tracks(
                group_truncations, all_truncations, strand, gene_name
            )

    def _create_simple_track(
        self,
        selected_truncations: List[Dict],
        all_truncations: List[Dict],
        strand: str,
        gene_name: str,
        description_suffix: str = "",
    ) -> List[Dict]:
        """Create a simple track with selected truncations plus all non-overlapping ones."""
        # Get all non-overlapping truncations
        regions = [
            (trunc["start"], trunc["end"], i) for i, trunc in enumerate(all_truncations)
        ]
        overlapping_groups = self._find_overlapping_groups(regions)

        track_truncations = selected_truncations.copy()

        # Add all single truncations
        for group in overlapping_groups:
            if len(group) == 1:
                orig_idx = group[0]
                for trunc in all_truncations:
                    if trunc["original_index"] == orig_idx:
                        if trunc not in track_truncations:
                            track_truncations.append(trunc)
                        break

        # Sort by position
        if strand == "+":
            track_truncations.sort(key=lambda x: x["start"])
        else:
            track_truncations.sort(key=lambda x: -x["end"])

        track_start = min(t["start"] for t in track_truncations)
        track_end = max(t["end"] for t in track_truncations)

        track = {
            "truncations": track_truncations,
            "start": track_start,
            "end": track_end,
            "gene_name": gene_name,
            "strand": strand,
            "track_id": f"{gene_name}_track_1",
            "description": f"Track with {len(track_truncations)} truncations"
            + (f" ({description_suffix})" if description_suffix else ""),
        }

        return [track]

    def _create_different_endpoint_tracks(
        self,
        group_truncations: List[Dict],
        all_truncations: List[Dict],
        strand: str,
        gene_name: str,
    ) -> List[Dict]:
        """Create tracks where truncations differ only in their endpoints."""
        self._debug_print("Creating tracks with different endpoints only")

        tracks = []

        for i, selected_trunc in enumerate(group_truncations):
            self._debug_print(
                f"Creating track {i + 1} with truncation: {selected_trunc['start']}-{selected_trunc['end']}"
            )

            track_truncations = [selected_trunc]

            # Add all non-overlapping truncations
            regions = [
                (trunc["start"], trunc["end"], j)
                for j, trunc in enumerate(all_truncations)
            ]
            overlapping_groups = self._find_overlapping_groups(regions)

            for group in overlapping_groups:
                if len(group) == 1:
                    orig_idx = group[0]
                    for trunc in all_truncations:
                        if trunc["original_index"] == orig_idx:
                            if trunc not in track_truncations:
                                track_truncations.append(trunc)
                            break

            # Sort by position
            if strand == "+":
                track_truncations.sort(key=lambda x: x["start"])
            else:
                track_truncations.sort(key=lambda x: -x["end"])

            track_start = min(t["start"] for t in track_truncations)
            track_end = max(t["end"] for t in track_truncations)

            track = {
                "truncations": track_truncations,
                "start": track_start,
                "end": track_end,
                "gene_name": gene_name,
                "strand": strand,
                "track_id": f"{gene_name}_track_{len(tracks) + 1}",
                "description": f"Track with {len(track_truncations)} truncations (different endpoint)",
            }

            tracks.append(track)
            self._debug_print(
                f"Created track: {track['track_id']} ({track_start}-{track_end})"
            )

        return tracks

    def _create_early_termination_tracks(
        self,
        group_truncations: List[Dict],
        all_truncations: List[Dict],
        strand: str,
        gene_name: str,
    ) -> List[Dict]:
        """Create tracks with early termination logic."""
        self._debug_print("Creating tracks with early termination logic")

        tracks = []

        for i, selected_trunc in enumerate(group_truncations):
            self._debug_print(
                f"Creating track {i + 1} with selected truncation: {selected_trunc['start']}-{selected_trunc['end']}"
            )

            # Start with the selected truncation
            track_truncations = [selected_trunc]

            # Add non-overlapping truncations, but apply early termination logic
            regions = [
                (trunc["start"], trunc["end"], j)
                for j, trunc in enumerate(all_truncations)
            ]
            overlapping_groups = self._find_overlapping_groups(regions)

            # Sort all truncations by position first
            if strand == "+":
                sorted_all = sorted(all_truncations, key=lambda x: x["start"])
            else:
                sorted_all = sorted(all_truncations, key=lambda x: -x["end"])

            # Find position of our selected truncation
            selected_pos = None
            for pos, trunc in enumerate(sorted_all):
                if trunc["original_index"] == selected_trunc["original_index"]:
                    selected_pos = pos
                    break

            if selected_pos is not None:
                if strand == "-":
                    # For negative strand, remove truncations to the left (higher indices) of the selected one
                    valid_truncations = sorted_all[: selected_pos + 1]
                    self._debug_print(
                        f"Negative strand: keeping truncations up to position {selected_pos}"
                    )
                else:
                    # For positive strand, remove truncations to the right (higher indices) of the selected one
                    valid_truncations = sorted_all[: selected_pos + 1]
                    self._debug_print(
                        f"Positive strand: keeping truncations up to position {selected_pos}"
                    )

                # Add valid non-overlapping truncations
                for trunc in valid_truncations:
                    if trunc not in track_truncations:
                        # Check if this truncation is in a single group (non-overlapping)
                        for group in overlapping_groups:
                            if len(group) == 1 and trunc["original_index"] == group[0]:
                                track_truncations.append(trunc)
                                break

            # Sort final track truncations
            if strand == "+":
                track_truncations.sort(key=lambda x: x["start"])
            else:
                track_truncations.sort(key=lambda x: -x["end"])

            track_start = min(t["start"] for t in track_truncations)
            track_end = max(t["end"] for t in track_truncations)

            track = {
                "truncations": track_truncations,
                "start": track_start,
                "end": track_end,
                "gene_name": gene_name,
                "strand": strand,
                "track_id": f"{gene_name}_track_{len(tracks) + 1}",
                "description": f"Track with {len(track_truncations)} truncations (early termination)",
            }

            tracks.append(track)
            self._debug_print(
                f"Created track: {track['track_id']} ({track_start}-{track_end}) with {len(track_truncations)} truncations"
            )

        return tracks

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
