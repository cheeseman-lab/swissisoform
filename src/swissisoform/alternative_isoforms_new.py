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
        """Load data from BED format file with flexible column support.

        Supports both standard 6-column and enhanced 7-column formats:
        - 6 columns: chrom, start, end, name, score, strand
        - 7 columns: chrom, start, end, name, score, strand, transcript_id

        Args:
            file_path: Path to the BED file
        """
        # Read BED format - auto-detect number of columns
        with open(file_path, "r") as f:
            first_line = f.readline().strip()
            if not first_line:
                raise ValueError("Empty BED file")

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
        else:
            raise ValueError(
                f"Unsupported BED format: {num_columns} columns. Expected 6 or 7."
            )

        # Read the full file
        self.start_sites = pd.read_csv(
            file_path,
            sep="\t",
            names=column_names,
            dtype=str,  # Read all as strings first to avoid parsing issues
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

    def check_data_quality(self) -> Dict:
        """Check data quality for loaded start sites.

        Returns:
            Dictionary with data quality statistics
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        # Get all unique genes
        all_genes = self.start_sites["gene_name"].unique()
        total_genes = len(all_genes)

        # Track genes with issues
        multiple_annotated = []
        missing_annotated = []

        print("Checking data quality...")
        print("=" * 50)

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
        print(f"Total genes: {total_genes}")
        print(
            f"Genes with multiple Annotated starts: {multiple_annotated_count} ({multiple_percent:.1f}%)"
        )
        print(
            f"Genes missing Annotated start: {missing_annotated_count} ({missing_percent:.1f}%)"
        )
        print()

        # Show examples of problematic genes
        if multiple_annotated:
            print("GENES WITH MULTIPLE ANNOTATED STARTS:")
            print("-" * 40)
            for gene_info in multiple_annotated[:10]:  # Show first 10
                gene = gene_info["gene"]
                positions = gene_info["positions"]
                efficiencies = gene_info["efficiencies"]
                print(f"{gene}: {len(positions)} Annotated starts")
                for pos, eff in zip(positions, efficiencies):
                    print(f"  Position {pos}, efficiency {eff:.3f}")

            if len(multiple_annotated) > 10:
                print(f"... and {len(multiple_annotated) - 10} more genes")
            print()

        if missing_annotated:
            print("GENES MISSING ANNOTATED START:")
            print("-" * 32)
            for gene in missing_annotated[:10]:  # Show first 10
                gene_sites = self.start_sites[self.start_sites["gene_name"] == gene]
                available_types = gene_sites["start_type"].unique()
                print(f"{gene}: has {list(available_types)}")

            if len(missing_annotated) > 10:
                print(f"... and {len(missing_annotated) - 10} more genes")
            print()

        # Additional statistics
        start_type_counts = self.start_sites["start_type"].value_counts()
        print("START TYPE DISTRIBUTION:")
        print("-" * 23)
        for start_type, count in start_type_counts.items():
            percentage = (count / len(self.start_sites)) * 100
            print(f"{start_type}: {count} ({percentage:.1f}%)")
        print()

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
            print("GENES WITH VERY DISTANT ANNOTATED STARTS (>10kb):")
            print("-" * 45)
            for gene_info in distant_annotated:
                gene = gene_info["gene"]
                distance = gene_info["distance"]
                positions = gene_info["positions"]
                print(f"{gene}: {distance:,} bp apart (positions: {positions})")
            print()

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
            gene_name: Name of the gene

        Returns:
            DataFrame of start sites for the gene
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        gene_sites = self.start_sites[self.start_sites["gene_name"] == gene_name].copy()
        return gene_sites.sort_values("start_position")

    def calculate_alternative_start_regions(
        self, gene_name: str, top_n_efficiency: Optional[int] = None
    ) -> pd.DataFrame:
        """Calculate truncation and extension regions for a gene.

        This function handles both:
        - Truncations: Regions removed when using downstream alternative starts
        - Extensions: Regions added when using upstream alternative starts

        Args:
            gene_name: Name of the gene
            top_n_efficiency: If specified, keep top N most efficient per TYPE (Truncated/Extended separately)

        Returns:
            DataFrame containing regions with metadata for both truncations and extensions
        """
        gene_sites = self.get_gene_start_sites(gene_name)
        if gene_sites.empty:
            return pd.DataFrame()

        # Find canonical (Annotated) start site
        canonical_sites = gene_sites[gene_sites["start_type"] == "Annotated"]
        if canonical_sites.empty:
            self._debug_print(
                f"WARNING: No Annotated start site found for {gene_name}. Skipping gene."
            )
            return pd.DataFrame()  # Return empty DataFrame if no canonical start

        # Always use the highest efficiency Annotated start by default
        canonical_start = canonical_sites.loc[canonical_sites["efficiency"].idxmax()]
        strand = canonical_start["strand"]
        canonical_pos = canonical_start["start_position"]

        self._debug_print(
            f"Gene {gene_name}: canonical start at {canonical_pos} (strand {strand}) "
            f"[BED: {canonical_start['start']}-{canonical_start['end']}] "
            f"efficiency: {canonical_start['efficiency']:.3f}"
        )

        # Get alternative start sites (Truncated and Extended)
        alternative_sites = gene_sites[
            gene_sites["start_type"].isin(["Truncated", "Extended"])
        ]

        if alternative_sites.empty:
            self._debug_print(f"No alternative start sites found for {gene_name}")
            return pd.DataFrame()

        # Apply efficiency filtering PER TYPE if requested
        if top_n_efficiency is not None:
            filtered_alternatives = []

            for start_type in ["Truncated", "Extended"]:
                type_sites = alternative_sites[
                    alternative_sites["start_type"] == start_type
                ]
                if not type_sites.empty:
                    # Get top N for this specific type
                    top_sites = type_sites.nlargest(top_n_efficiency, "efficiency")
                    filtered_alternatives.append(top_sites)
                    self._debug_print(
                        f"Filtered to top {min(top_n_efficiency, len(type_sites))} "
                        f"{start_type} starts (from {len(type_sites)} total)"
                    )

            if filtered_alternatives:
                alternative_sites = pd.concat(filtered_alternatives, ignore_index=True)
            else:
                alternative_sites = pd.DataFrame()  # No alternatives after filtering

            if alternative_sites.empty:
                self._debug_print(
                    f"No alternative start sites remaining after efficiency filtering"
                )
                return pd.DataFrame()

        # Calculate regions for both truncations and extensions
        regions = []

        for _, alt_site in alternative_sites.iterrows():
            bed_annotation = alt_site["start_type"]  # "Truncated" or "Extended"
            alt_pos = alt_site["start_position"]

            self._debug_print(
                f"  Alternative start: {alt_pos} ({bed_annotation}) "
                f"[BED: {alt_site['start']}-{alt_site['end']}] "
                f"efficiency: {alt_site['efficiency']:.3f}"
            )

            # Determine biological effect based on relative positions and strand
            if strand == "+":
                if alt_pos > canonical_pos:
                    # Downstream alternative start = N-terminal truncation
                    region_start = canonical_pos
                    region_end = alt_pos - 1
                    region_type = "truncation"
                    biological_effect = "removes_N_terminal_sequence"
                else:
                    # Upstream alternative start = N-terminal extension
                    region_start = alt_pos
                    region_end = canonical_pos - 1
                    region_type = "extension"
                    biological_effect = "adds_N_terminal_sequence"

            else:  # strand == "-"
                if alt_pos < canonical_pos:
                    # Lower genomic coordinate = downstream in transcript = truncation
                    region_start = alt_pos
                    region_end = canonical_start[
                        "start"
                    ]  # Use BED start coordinate for the end
                    region_type = "truncation"
                    biological_effect = "removes_N_terminal_sequence"
                else:
                    # Higher genomic coordinate = upstream in transcript = extension
                    region_start = canonical_start[
                        "start"
                    ]  # Use BED start coordinate for the start
                    region_end = alt_pos - 1
                    region_type = "extension"
                    biological_effect = "adds_N_terminal_sequence"

            region_length = region_end - region_start + 1

            # Validate that BED annotation matches biological prediction
            annotation_matches_biology = (
                region_type == "truncation" and bed_annotation == "Truncated"
            ) or (region_type == "extension" and bed_annotation == "Extended")

            if not annotation_matches_biology:
                self._debug_print(
                    f"WARNING: BED annotation '{bed_annotation}' doesn't match "
                    f"predicted biological effect '{region_type}' for {gene_name} "
                    f"alt start at {alt_pos} (canonical at {canonical_pos}, strand {strand})"
                )

            # Skip invalid regions
            if region_length <= 0:
                self._debug_print(
                    f"Warning: Invalid region length {region_length} for {alt_site['name']}"
                )
                continue

            self._debug_print(
                f"  -> {region_type}: {region_start}-{region_end} (length: {region_length})"
            )

            regions.append(
                {
                    "chromosome": alt_site["chrom"],
                    "region_start": region_start,
                    "region_end": region_end,
                    "region_type": region_type,  # "truncation" or "extension"
                    "region_length": region_length,
                    "biological_effect": biological_effect,
                    "bed_annotation": bed_annotation,  # "Truncated" or "Extended" from BED
                    "annotation_matches_biology": annotation_matches_biology,
                    "strand": strand,
                    "gene_name": gene_name,
                    "gene_id": alt_site["gene_id"],
                    "canonical_start_pos": canonical_pos,
                    "alternative_start_pos": alt_pos,
                    "alternative_start_type": alt_site["start_type"],
                    "alternative_start_codon": alt_site["start_codon"],
                    "efficiency": alt_site["efficiency"],
                    "canonical_efficiency": canonical_start["efficiency"],
                    # Create unique identifier
                    "region_id": f"{gene_name}_{region_type}_{alt_site['start_codon']}_{alt_pos}",
                    "track_id": f"{gene_name}_{region_type}_{alt_pos}",
                }
            )

        result_df = pd.DataFrame(regions)

        # Debug output
        if not result_df.empty:
            self._debug_print(f"Generated regions for {gene_name} (strand {strand}):")
            for region_type in ["truncation", "extension"]:
                type_regions = result_df[result_df["region_type"] == region_type]
                if not type_regions.empty:
                    self._debug_print(
                        f"  {region_type.upper()}S ({len(type_regions)}):"
                    )
                    for _, region in type_regions.iterrows():
                        self._debug_print(
                            f"    {region['region_start']}-{region['region_end']} "
                            f"(length: {region['region_length']}, "
                            f"alt_start: {region['alternative_start_pos']}, "
                            f"efficiency: {region['efficiency']:.3f}, "
                            f"BED: {region['bed_annotation']}, "
                            f"matches: {region['annotation_matches_biology']})"
                        )
        else:
            self._debug_print(f"No valid regions generated for {gene_name}")

        return result_df

    def get_visualization_features(
        self,
        gene_name: str,
        top_n_per_type: Optional[int] = None,
        include_extensions: bool = True,
        include_truncations: bool = True,
    ) -> pd.DataFrame:
        """Get features formatted for visualization, including both truncations and extensions.

        Args:
            gene_name: Name of the gene to get features for
            top_n_per_type: If specified, keep top N most efficient per TYPE (Truncated/Extended separately)
            include_extensions: Whether to include extension regions
            include_truncations: Whether to include truncation regions

        Returns:
            DataFrame: Features formatted for visualization
        """
        if self.start_sites.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        # Use the updated function with per-type filtering
        regions = self.calculate_alternative_start_regions(gene_name, top_n_per_type)

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
                "score": region["efficiency"],
                "strand": region["strand"],
                "frame": ".",
                "gene_id": region["gene_id"],
                "transcript_id": region["track_id"],
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
