"""Alternative start sites data handling module.

This module provides functionality for loading and processing alternative
start sites from BED format files, calculating truncation/extension regions,
and supporting efficiency-based filtering.
"""

import pandas as pd
import numpy as np
import logging
from typing import Optional, Dict, List, Tuple, Set
from itertools import product

logger = logging.getLogger(__name__)


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

    def __init__(self, debug: bool = False):
        """Initialize an empty start sites handler.

        Args:
            debug (bool): If True, enables debug printing.
        """
        self.start_sites = pd.DataFrame()
        self.debug = debug

    def _debug_print(self, message: str) -> None:
        """Print debug message if debug mode is enabled.

        Args:
            message (str): Message to print.
        """
        if self.debug:
            logger.debug(message)

    def load_bed(self, file_path: str) -> None:
        """Load data from BED format file with flexible column support.

        Supports standard 6-column, enhanced 7-column, and dual-transcript 8-column formats:
        - 6 columns: chrom, start, end, name, score, strand
        - 7 columns: chrom, start, end, name, score, strand, transcript_id
        - 8 columns: chrom, start, end, name, score, strand, ensembl_transcript_id, refseq_transcript_id

        Args:
            file_path (str): Path to the BED file.

        Raises:
            ValueError: If the BED file is empty or has an unsupported format.
        """
        # Read BED format - auto-detect number of columns, skipping comment lines
        with open(file_path, "r") as f:
            first_line = None
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    first_line = line
                    break

            if not first_line:
                raise ValueError("Empty BED file or no data lines found")

            num_columns = len(first_line.split("\t"))

        if num_columns == 6:
            # Standard 6-column format
            column_names = ["chrom", "start", "end", "name", "score", "strand"]
            self._debug_print("Loading standard 6-column BED format")

        elif num_columns == 7:
            # Enhanced 7-column format with transcript_id
            column_names = [
                "chrom",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "transcript_id",
            ]
            self._debug_print(
                "Loading enhanced 7-column BED format with transcript IDs"
            )

        elif num_columns == 8:
            # Dual-transcript 8-column format with both Ensembl and RefSeq IDs
            column_names = [
                "chrom",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "ensembl_transcript_id",
                "refseq_transcript_id",
            ]
            self._debug_print(
                "Loading dual-transcript 8-column BED format with Ensembl and RefSeq IDs"
            )

        else:
            raise ValueError(
                f"Unsupported BED format: {num_columns} columns. Expected 6, 7, or 8."
            )

        # Read the full file, skipping comment lines
        self.start_sites = pd.read_csv(
            file_path,
            sep="\t",
            names=column_names,
            dtype=str,  # Read all as strings first to avoid parsing issues
            comment="#",  # Skip lines starting with #
        )

        # Parse the name field which contains gene information
        def parse_name(name: str) -> Dict[str, str]:
            # Handle cases where name might be numeric (shouldn't happen but safety check)
            if not isinstance(name, str):
                name = str(name)

            parts = name.split("_")
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

        # Calculate start codon positions correctly based on strand
        def get_start_position(row):
            if row["strand"] == "+":
                return row["start"]  # Start codon at 'start' coordinate for + strand
            else:
                return row["end"]  # Start codon at 'end' coordinate for - strand

        self.start_sites["start_position"] = self.start_sites.apply(
            get_start_position, axis=1
        )

        self._debug_print(f"Loaded {len(self.start_sites)} start sites")
        self._debug_print(f"Columns: {list(self.start_sites.columns)}")
        self._debug_print(
            f"Start types: {self.start_sites['start_type'].value_counts().to_dict()}"
        )

        # Report transcript information if available
        if "transcript_id" in self.start_sites.columns:
            transcript_assignments = (self.start_sites["transcript_id"] != "NA").sum()
            self._debug_print(
                f"Transcript assignments: {transcript_assignments}/{len(self.start_sites)}"
            )
        elif "ensembl_transcript_id" in self.start_sites.columns:
            # For 8-column format, report both transcript types
            ensembl_assignments = (
                self.start_sites["ensembl_transcript_id"] != "NA"
            ).sum()
            refseq_assignments = (
                self.start_sites["refseq_transcript_id"] != "NA"
            ).sum()

            self._debug_print(
                f"Ensembl transcript assignments: {ensembl_assignments}/{len(self.start_sites)}"
            )
            self._debug_print(
                f"RefSeq transcript assignments: {refseq_assignments}/{len(self.start_sites)}"
            )

            # For backward compatibility, create a transcript_id column using Ensembl ID
            self.start_sites["transcript_id"] = self.start_sites[
                "ensembl_transcript_id"
            ]

    def check_data_quality(self) -> Dict:
        """Check data quality for loaded start sites.

        Returns:
            Dict: Dictionary with data quality statistics, including counts and lists of genes with multiple or missing annotated starts, distant annotated starts, and start type distribution.

        Raises:
            ValueError: If no data is loaded.
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        # Get all unique genes
        all_genes = self.start_sites["gene_name"].unique()
        total_genes = len(all_genes)

        # Track genes with issues
        multiple_annotated = []
        missing_annotated = []

        logger.info("Checking data quality...")
        logger.info("=" * 50)

        # Check each gene
        for gene_name in all_genes:
            gene_sites = self.start_sites[self.start_sites["gene_name"] == gene_name]
            annotated_sites = gene_sites[gene_sites["start_type"] == "Annotated"]

            if len(annotated_sites) > 1:
                multiple_annotated.append(
                    {
                        "gene": gene_name,
                        "count": len(annotated_sites),
                        "positions": annotated_sites["start_position"].tolist(),
                        "efficiencies": annotated_sites["efficiency"].tolist(),
                    }
                )
            elif len(annotated_sites) == 0:
                missing_annotated.append(gene_name)

        # Calculate percentages
        multiple_annotated_count = len(multiple_annotated)
        missing_annotated_count = len(missing_annotated)

        multiple_percent = (multiple_annotated_count / total_genes) * 100
        missing_percent = (missing_annotated_count / total_genes) * 100

        # Print summary
        logger.info(f"Total genes: {total_genes}")
        logger.info(
            f"Genes with multiple Annotated starts: {multiple_annotated_count} ({multiple_percent:.1f}%)"
        )
        logger.info(
            f"Genes missing Annotated start: {missing_annotated_count} ({missing_percent:.1f}%)"
        )
        logger.info("")

        # Show examples of problematic genes
        if multiple_annotated:
            logger.info("GENES WITH MULTIPLE ANNOTATED STARTS:")
            logger.info("-" * 40)
            for gene_info in multiple_annotated[:10]:  # Show first 10
                gene = gene_info["gene"]
                positions = gene_info["positions"]
                efficiencies = gene_info["efficiencies"]
                logger.info(f"{gene}: {len(positions)} Annotated starts")
                for pos, eff in zip(positions, efficiencies):
                    logger.info(f"  Position {pos}, efficiency {eff:.3f}")

            if len(multiple_annotated) > 10:
                logger.info(f"... and {len(multiple_annotated) - 10} more genes")
            logger.info("")

        if missing_annotated:
            logger.info("GENES MISSING ANNOTATED START:")
            logger.info("-" * 32)
            for gene in missing_annotated[:10]:  # Show first 10
                gene_sites = self.start_sites[self.start_sites["gene_name"] == gene]
                available_types = gene_sites["start_type"].unique()
                logger.info(f"{gene}: has {list(available_types)}")

            if len(missing_annotated) > 10:
                logger.info(f"... and {len(missing_annotated) - 10} more genes")
            logger.info("")

        # Additional statistics
        start_type_counts = self.start_sites["start_type"].value_counts()
        logger.info("START TYPE DISTRIBUTION:")
        logger.info("-" * 23)
        for start_type, count in start_type_counts.items():
            percentage = (count / len(self.start_sites)) * 100
            logger.info(f"{start_type}: {count} ({percentage:.1f}%)")
        logger.info("")

        # Check for very distant Annotated starts within genes
        distant_annotated = []
        for gene_info in multiple_annotated:
            positions = gene_info["positions"]
            if len(positions) >= 2:
                min_pos = min(positions)
                max_pos = max(positions)
                distance = max_pos - min_pos
                if distance > 10000:  # More than 10kb apart
                    distant_annotated.append(
                        {
                            "gene": gene_info["gene"],
                            "distance": distance,
                            "positions": positions,
                        }
                    )

        if distant_annotated:
            logger.info("GENES WITH VERY DISTANT ANNOTATED STARTS (>10kb):")
            logger.info("-" * 45)
            for gene_info in distant_annotated:
                gene = gene_info["gene"]
                distance = gene_info["distance"]
                positions = gene_info["positions"]
                logger.info(f"{gene}: {distance:,} bp apart (positions: {positions})")
            logger.info("")

        # Return structured data
        return {
            "total_genes": total_genes,
            "multiple_annotated": {
                "count": multiple_annotated_count,
                "percentage": multiple_percent,
                "genes": multiple_annotated,
            },
            "missing_annotated": {
                "count": missing_annotated_count,
                "percentage": missing_percent,
                "genes": missing_annotated,
            },
            "distant_annotated": {
                "count": len(distant_annotated),
                "genes": distant_annotated,
            },
            "start_type_distribution": start_type_counts.to_dict(),
        }

    def get_gene_start_sites(self, gene_name: str) -> pd.DataFrame:
        """Get all start sites for a specific gene.

        Args:
            gene_name (str): Name of the gene.

        Returns:
            pd.DataFrame: DataFrame of start sites for the gene, sorted by start position.

        Raises:
            ValueError: If no data is loaded.
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        gene_sites = self.start_sites[self.start_sites["gene_name"] == gene_name].copy()
        return gene_sites.sort_values("start_position")

    def calculate_alternative_start_regions(
        self, gene_name: str, top_n_per_type_per_transcript: Optional[int] = None
    ) -> pd.DataFrame:
        """Calculate truncation and extension regions grouped by transcript.

        Groups by transcript_id so each group contains:
        - 1 annotated start (canonical for that transcript)
        - N alternative starts (truncated/extended for that transcript)

        Args:
            gene_name (str): Name of the gene.
            top_n_per_type_per_transcript (Optional[int]): If specified, keep top N most efficient per TYPE per transcript.

        Returns:
            pd.DataFrame: DataFrame containing regions with transcript-specific canonical starts.
        """
        gene_sites = self.get_gene_start_sites(gene_name)
        if gene_sites.empty:
            return pd.DataFrame()

        self._debug_print(
            f"Processing {gene_name} with {len(gene_sites)} total start sites"
        )

        # Group by transcript_id (skip "NA" assignments)
        valid_transcripts = gene_sites[gene_sites["transcript_id"] != "NA"]

        if valid_transcripts.empty:
            self._debug_print(f"No valid transcript assignments for {gene_name}")
            return pd.DataFrame()

        all_regions = []

        # Process each transcript group separately
        for transcript_id, transcript_sites in valid_transcripts.groupby(
            "transcript_id"
        ):
            self._debug_print(f"\nProcessing transcript {transcript_id}:")
            self._debug_print(f"  Sites in group: {len(transcript_sites)}")

            # Find annotated start for this transcript
            annotated_sites = transcript_sites[
                transcript_sites["start_type"] == "Annotated"
            ]
            alternative_sites = transcript_sites[
                transcript_sites["start_type"].isin(["Truncated", "Extended"])
            ]

            self._debug_print(f"  Annotated starts: {len(annotated_sites)}")
            self._debug_print(f"  Alternative starts: {len(alternative_sites)}")

            # Skip if no annotated start for this transcript
            if annotated_sites.empty:
                self._debug_print(
                    f"  WARNING: No annotated start found for transcript {transcript_id}"
                )
                continue

            # Use the annotated start for this transcript (should be only one)
            if len(annotated_sites) > 1:
                self._debug_print(
                    f"  WARNING: Multiple annotated starts for transcript {transcript_id}, using highest efficiency"
                )
                canonical_start = annotated_sites.loc[
                    annotated_sites["efficiency"].idxmax()
                ]
            else:
                canonical_start = annotated_sites.iloc[0]

            canonical_pos = canonical_start["start_position"]
            strand = canonical_start["strand"]

            self._debug_print(
                f"  Canonical start: {canonical_pos} (efficiency: {canonical_start['efficiency']:.3f})"
            )

            # Skip if no alternatives for this transcript
            if alternative_sites.empty:
                self._debug_print(
                    f"  No alternative starts for transcript {transcript_id}"
                )
                continue

            # Apply efficiency filtering PER TYPE within this transcript
            if top_n_per_type_per_transcript is not None:
                filtered_alternatives = []

                for start_type in ["Truncated", "Extended"]:
                    type_sites = alternative_sites[
                        alternative_sites["start_type"] == start_type
                    ]
                    if not type_sites.empty:
                        top_sites = type_sites.nlargest(
                            top_n_per_type_per_transcript, "efficiency"
                        )
                        filtered_alternatives.append(top_sites)
                        self._debug_print(
                            f"  Filtered to top {min(top_n_per_type_per_transcript, len(type_sites))} {start_type} starts"
                        )

                if filtered_alternatives:
                    alternative_sites = pd.concat(
                        filtered_alternatives, ignore_index=True
                    )
                else:
                    alternative_sites = pd.DataFrame()

            if alternative_sites.empty:
                self._debug_print(f"  No alternative starts remaining after filtering")
                continue

            # Calculate regions for this transcript group
            for _, alt_site in alternative_sites.iterrows():
                alt_pos = alt_site["start_position"]
                alt_type = alt_site["start_type"]

                self._debug_print(
                    f"    Alternative: {alt_pos} ({alt_type}, efficiency: {alt_site['efficiency']:.3f})"
                )

                # Calculate region based on strand and positions
                if strand == "+":
                    if alt_pos > canonical_pos:
                        # Downstream = truncation
                        region_start = canonical_pos
                        region_end = alt_pos - 1
                        region_type = "truncation"
                        biological_effect = "removes_N_terminal_sequence"
                    else:
                        # Upstream = extension
                        region_start = alt_pos
                        region_end = canonical_pos - 1
                        region_type = "extension"
                        biological_effect = "adds_N_terminal_sequence"
                else:  # strand == "-"
                    if alt_pos < canonical_pos:
                        # Lower genomic coord = downstream in transcript = truncation
                        region_start = alt_pos + 1
                        region_end = canonical_pos
                        region_type = "truncation"
                        biological_effect = "removes_N_terminal_sequence"
                    else:
                        # Higher genomic coord = upstream in transcript = extension
                        region_start = canonical_pos + 1
                        region_end = alt_pos
                        region_type = "extension"
                        biological_effect = "adds_N_terminal_sequence"

                region_length = region_end - region_start + 1

                # Skip invalid regions
                if region_length <= 0:
                    self._debug_print(
                        f"    WARNING: Invalid region length {region_length}"
                    )
                    continue

                # Check if BED annotation matches biological prediction
                annotation_matches = (
                    region_type == "truncation" and alt_type == "Truncated"
                ) or (region_type == "extension" and alt_type == "Extended")

                if not annotation_matches:
                    self._debug_print(
                        f"    WARNING: BED annotation '{alt_type}' doesn't match predicted '{region_type}'"
                    )

                self._debug_print(
                    f"    -> {region_type}: {region_start}-{region_end} (length: {region_length})"
                )

                # Create region record
                region = {
                    "chromosome": alt_site["chrom"],
                    "region_start": region_start,
                    "region_end": region_end,
                    "region_type": region_type,
                    "region_length": region_length,
                    "biological_effect": biological_effect,
                    "bed_annotation": alt_type,
                    "annotation_matches_biology": annotation_matches,
                    "strand": strand,
                    "gene_name": gene_name,
                    "gene_id": alt_site["gene_id"],
                    "transcript_id": transcript_id,  # Add transcript ID to region
                    "canonical_start_pos": canonical_pos,
                    "alternative_start_pos": alt_pos,
                    "alternative_start_type": alt_site["start_type"],
                    "alternative_start_codon": alt_site["start_codon"],
                    "efficiency": alt_site["efficiency"],
                    "canonical_efficiency": canonical_start["efficiency"],
                    "region_id": f"{gene_name}_{transcript_id}_{region_type}_{alt_site['start_codon']}_{alt_pos}",
                    "track_id": f"{gene_name}_{transcript_id}_{region_type}_{alt_pos}",
                }

                all_regions.append(region)

        result_df = pd.DataFrame(all_regions)

        if not result_df.empty:
            self._debug_print(f"\nGenerated {len(result_df)} regions for {gene_name}:")
            for transcript_id in result_df["transcript_id"].unique():
                transcript_regions = result_df[
                    result_df["transcript_id"] == transcript_id
                ]
                truncations = len(
                    transcript_regions[
                        transcript_regions["region_type"] == "truncation"
                    ]
                )
                extensions = len(
                    transcript_regions[transcript_regions["region_type"] == "extension"]
                )
                self._debug_print(
                    f"  {transcript_id}: {truncations} truncations, {extensions} extensions"
                )
        else:
            self._debug_print(f"No valid regions generated for {gene_name}")

        return result_df

    def get_transcript_groups(self, gene_name: str) -> Dict[str, Dict]:
        """Get summary of transcript groups for a gene.

        Args:
            gene_name (str): Name of the gene.

        Returns:
            Dict[str, Dict]: Dictionary mapping transcript_id to group info, including counts and efficiency ranges.
        """
        gene_sites = self.get_gene_start_sites(gene_name)

        if gene_sites.empty:
            return {}

        groups = {}
        valid_transcripts = gene_sites[gene_sites["transcript_id"] != "NA"]

        for transcript_id, transcript_sites in valid_transcripts.groupby(
            "transcript_id"
        ):
            annotated = transcript_sites[transcript_sites["start_type"] == "Annotated"]
            alternatives = transcript_sites[
                transcript_sites["start_type"].isin(["Truncated", "Extended"])
            ]

            groups[transcript_id] = {
                "total_sites": len(transcript_sites),
                "annotated_count": len(annotated),
                "alternative_count": len(alternatives),
                "alternative_types": alternatives["start_type"]
                .value_counts()
                .to_dict(),
                "efficiency_range": [
                    alternatives["efficiency"].min()
                    if not alternatives.empty
                    else None,
                    alternatives["efficiency"].max()
                    if not alternatives.empty
                    else None,
                ],
                "has_complete_data": len(annotated) > 0 and len(alternatives) > 0,
            }

        return groups

    def get_mutation_features(
        self,
        gene_name: str,
        top_n_per_type_per_transcript: Optional[int] = None,
        include_extensions: bool = True,
        include_truncations: bool = True,
    ) -> pd.DataFrame:
        """Get features formatted for mutations, including both truncations and extensions.

        Args:
            gene_name (str): Name of the gene to get features for.
            top_n_per_type_per_transcript (Optional[int]): If specified, keep top N most efficient per TYPE (Truncated/Extended separately).
            include_extensions (bool): Whether to include extension regions.
            include_truncations (bool): Whether to include truncation regions.

        Returns:
            pd.DataFrame: Features formatted for visualization.

        Raises:
            ValueError: If no data is loaded.
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        # Use the updated function with per-type filtering
        regions = self.calculate_alternative_start_regions(
            gene_name, top_n_per_type_per_transcript
        )

        if regions.empty:
            return pd.DataFrame()

        # Filter by region type if requested
        if not include_extensions:
            regions = regions[regions["region_type"] != "extension"]
        if not include_truncations:
            regions = regions[regions["region_type"] != "truncation"]

        if regions.empty:
            return pd.DataFrame()

        # Convert to visualization format
        viz_features = []

        for _, region in regions.iterrows():
            feature = {
                "chromosome": region["chromosome"],
                "source": "ribosome_profiling",
                "feature_type": f"alternative_start_{region['region_type']}",
                "start": region["region_start"],
                "end": region["region_end"],
                "mutation_start": region["region_start"] - 3
                if region["region_type"] == "truncation" and region["strand"] == "-"
                else region["region_start"],
                "mutation_end": region["region_end"] + 3
                if region["region_type"] == "truncation" and region["strand"] == "+"
                else region["region_end"],
                "alternative_start_pos": region["alternative_start_pos"],
                "canonical_start_pos": region["canonical_start_pos"],
                "score": region["efficiency"],
                "strand": region["strand"],
                "frame": ".",
                "gene_id": region["gene_id"],
                "transcript_id": region["transcript_id"],
                "refseq_transcript_id": region.get("refseq_transcript_id", "NA"),
                "gene_name": region["gene_name"],
                "start_codon": region["alternative_start_codon"],
                "name": region["region_id"],
                "region_type": region["region_type"],
                "region_length": region["region_length"],
                "biological_effect": region["biological_effect"],
                "bed_annotation": region["bed_annotation"],
                "annotation_matches_biology": region["annotation_matches_biology"],
                "efficiency": region["efficiency"],
                "canonical_efficiency": region["canonical_efficiency"],
                "efficiency_ratio": region["efficiency"]
                / region["canonical_efficiency"]
                if region["canonical_efficiency"] > 0
                else 0,
            }
            viz_features.append(feature)

        return pd.DataFrame(viz_features)

    def get_translation_features(
        self,
        gene_name: str,
        top_n_per_type_per_transcript: Optional[int] = None,
        include_extensions: bool = True,
        include_truncations: bool = True,
    ) -> pd.DataFrame:
        """Get features formatted for translation, including both truncations and extensions.

        Args:
            gene_name (str): Name of the gene to get features for.
            top_n_per_type_per_transcript (Optional[int]): If specified, keep top N most efficient per TYPE.
            include_extensions (bool): Whether to include extension regions.
            include_truncations (bool): Whether to include truncation regions.

        Returns:
            pd.DataFrame: Features formatted for translation with original coordinates.

        Raises:
            ValueError: If no data is loaded.
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        # Use the updated function with per-type filtering
        regions = self.calculate_alternative_start_regions(
            gene_name, top_n_per_type_per_transcript
        )

        if regions.empty:
            return pd.DataFrame()

        # Filter by region type if requested
        if not include_extensions:
            regions = regions[regions["region_type"] != "extension"]
        if not include_truncations:
            regions = regions[regions["region_type"] != "truncation"]

        if regions.empty:
            return pd.DataFrame()

        # Convert to translation format - NO MUTATION COORDINATE EXTENSIONS
        translation_features = []

        for _, region in regions.iterrows():
            feature = {
                "chromosome": region["chromosome"],
                "source": "ribosome_profiling",
                "feature_type": f"alternative_start_{region['region_type']}",
                "start": region["region_start"],  # Original coordinates only
                "end": region["region_end"],  # Original coordinates only
                "alternative_start_pos": region["alternative_start_pos"],
                "canonical_start_pos": region["canonical_start_pos"],
                "score": region["efficiency"],
                "strand": region["strand"],
                "frame": ".",
                "gene_id": region["gene_id"],
                "transcript_id": region["transcript_id"],
                "refseq_transcript_id": region.get("refseq_transcript_id", "NA"),
                "gene_name": region["gene_name"],
                "start_codon": region["alternative_start_codon"],
                "name": region["region_id"],
                "region_type": region["region_type"],
                "region_length": region["region_length"],
                "biological_effect": region["biological_effect"],
                "bed_annotation": region["bed_annotation"],
                "annotation_matches_biology": region["annotation_matches_biology"],
                "efficiency": region["efficiency"],
                "canonical_efficiency": region["canonical_efficiency"],
                "efficiency_ratio": region["efficiency"]
                / region["canonical_efficiency"]
                if region["canonical_efficiency"] > 0
                else 0,
            }
            translation_features.append(feature)

        return pd.DataFrame(translation_features)

    def get_gene_list(self) -> List[str]:
        """Get list of all genes in the dataset.

        Returns:
            List[str]: List of gene names.

        Raises:
            ValueError: If no data is loaded.
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        return sorted(self.start_sites["gene_name"].unique().tolist())

    def get_stats(self) -> Dict:
        """Get basic statistics about the loaded start sites.

        Returns:
            Dict: Dictionary containing various statistics such as total start sites, unique genes, chromosomes, start types, start codons, and efficiency stats.

        Raises:
            ValueError: If no data is loaded.
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
        top_n_per_type: Optional[int] = None,
    ) -> "AlternativeIsoform":
        """Filter start sites with gene-level efficiency logic.

        Args:
            min_efficiency (Optional[float]): Minimum efficiency threshold to filter start sites.
            top_n_per_gene (Optional[int]): Keep top N alternatives per gene (all annotated retained).
            top_n_per_type (Optional[int]): Keep top N truncated and top N extended per gene.

        Returns:
            AlternativeIsoform: New instance containing filtered start sites.
        """
        filtered_data = self.start_sites.copy()

        if min_efficiency is not None:
            filtered_data = filtered_data[filtered_data["efficiency"] >= min_efficiency]

        if top_n_per_gene is not None:
            # Keep all annotated, get top N alternatives PER GENE
            annotated = filtered_data[filtered_data["start_type"] == "Annotated"]
            alternatives = filtered_data[filtered_data["start_type"] != "Annotated"]

            top_alternatives = (
                alternatives.groupby("gene_name")
                .apply(lambda x: x.nlargest(top_n_per_gene, "efficiency"))
                .reset_index(drop=True)
            )

            filtered_data = pd.concat([annotated, top_alternatives], ignore_index=True)

        if top_n_per_type is not None:
            # Keep top N truncated + top N extended PER GENE
            annotated = filtered_data[filtered_data["start_type"] == "Annotated"]
            alternatives = filtered_data[filtered_data["start_type"] != "Annotated"]

            top_by_type = []
            for gene_name in alternatives["gene_name"].unique():
                gene_alts = alternatives[alternatives["gene_name"] == gene_name]

                for start_type in ["Truncated", "Extended"]:
                    type_sites = gene_alts[gene_alts["start_type"] == start_type]
                    if not type_sites.empty:
                        top_sites = type_sites.nlargest(top_n_per_type, "efficiency")
                        top_by_type.append(top_sites)

            if top_by_type:
                filtered_alternatives = pd.concat(top_by_type, ignore_index=True)
                filtered_data = pd.concat(
                    [annotated, filtered_alternatives], ignore_index=True
                )

        # Return filtered instance
        filtered_instance = AlternativeIsoform(debug=self.debug)
        filtered_instance.start_sites = filtered_data
        return filtered_instance
