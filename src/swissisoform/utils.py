"""BED cleanup utilities for ribosome profiling data.

This module provides functions for processing ribosome profiling BED files
with pre-assigned transcript IDs and generating protein sequences.
"""

import pandas as pd
from typing import List, Dict, Optional, Union, Tuple, Set
import logging
from pathlib import Path
from collections import defaultdict
import re

logger = logging.getLogger(__name__)


# NOTE: TranscriptDatabase, BedProcessor, and related classes removed
# These were used for the old position-based transcript matching workflow.
# Current workflow uses pre-assigned transcript IDs from ribosome profiling data.


def find_refseq_gtf_file(genome_data_dir: Union[str, Path]) -> Optional[Path]:
    """Find available RefSeq GTF file in the genome data directory.

    Args:
        genome_data_dir (Union[str, Path]): Directory to search for RefSeq GTF files.

    Returns:
        Optional[Path]: Path to RefSeq GTF file if found, None otherwise.
    """
    genome_data_dir = Path(genome_data_dir)

    possible_files = [
        "GRCh38_latest_genomic.gtf",
        "ncbiRefSeq.gtf",
        "refseq_genomic.gtf",
        "GRCh38_refseq.gtf",
    ]

    for filename in possible_files:
        gtf_path = genome_data_dir / filename
        if gtf_path.exists():
            return gtf_path

    return None


def print_translation_summary(
    dataset: pd.DataFrame,
    successful_genes: int,
    total_genes: int,
    failed_genes: List[str],
    mutations_mode: bool,
    output_dir: str,
) -> None:
    """Print comprehensive summary of protein sequence generation results.

    Args:
        dataset: Generated protein sequence dataset
        successful_genes: Number of successfully processed genes
        total_genes: Total number of genes attempted
        failed_genes: List of gene names that failed processing
        mutations_mode: Whether mutations were included
        output_dir: Output directory path
    """
    logger.info(f"\nProtein Sequence Generation Summary:")
    logger.info(f"  Total genes processed: {total_genes}")

    # Status breakdown
    logger.info(f"\n  Status breakdown:")
    logger.info(f"  Success: {successful_genes}")
    if failed_genes:
        logger.info(f"  Failed: {len(failed_genes)}")
    else:
        logger.info(f"  Failed: 0")

    # Dataset statistics
    if not dataset.empty:
        logger.info(f"\n  Sequence Generation:")
        logger.info(f"  Total sequences generated: {len(dataset):,}")

        # Calculate transcript-isoform pairs
        genes_with_data = dataset["gene"].nunique()
        if mutations_mode and "variant_type" in dataset.columns:
            # Count unique transcript-isoform pairs (canonical + alternative base pairs)
            base_sequences = dataset[
                dataset["variant_type"].isin(
                    ["canonical", "alternative"]
                )  # Updated terminology
            ]
            unique_pairs = (
                len(base_sequences) // 2
                if len(base_sequences) % 2 == 0
                else (len(base_sequences) + 1) // 2
            )
        else:
            # For pairs mode, total sequences / 2 = pairs
            unique_pairs = (
                len(dataset) // 2 if len(dataset) % 2 == 0 else (len(dataset) + 1) // 2
            )

        logger.info(
            f"  Transcript-isoform pairs: {unique_pairs}"
        )  # Updated terminology
        logger.info(
            f"  Average sequences per gene: {len(dataset) / genes_with_data:.1f}"
        )

        # Mode-specific breakdown
        if mutations_mode and "variant_type" in dataset.columns:
            logger.info(f"\n  Sequence breakdown:")
            type_counts = dataset["variant_type"].value_counts()
            for variant_type, count in type_counts.items():
                percentage = (count / len(dataset)) * 100
                logger.info(f"  {variant_type}: {count:,} ({percentage:.1f}%)")

            # Mutation-specific statistics
            if "canonical_mutated" in type_counts:
                mutation_data = dataset[dataset["variant_type"] == "canonical_mutated"]
                if not mutation_data.empty:
                    logger.info(f"\n  Mutation Analysis:")
                    logger.info(f"  Total mutations integrated: {len(mutation_data)}")
                    logger.info(
                        f"  Unique mutation positions: {mutation_data['mutation_position'].nunique()}"
                    )

                    if "aa_change" in mutation_data.columns:
                        aa_changes = mutation_data["aa_change"].dropna()
                        silent_mutations = len(mutation_data) - len(aa_changes)
                        logger.info(f"  Mutations with AA changes: {len(aa_changes)}")
                        if silent_mutations > 0:
                            logger.info(f"  Silent mutations: {silent_mutations}")

                    # Breakdown by impact type if available
                    if "mutation_impact" in mutation_data.columns:
                        impact_counts = mutation_data["mutation_impact"].value_counts()
                        logger.info(f"")
                        logger.info(f"  Breakdown by impact type:")
                        for impact_type, count in impact_counts.items():
                            percentage = (count / len(mutation_data)) * 100
                            logger.info(f"  {impact_type}: {count} ({percentage:.1f}%)")

                    logger.info(
                        f"  Average mutations per gene: {len(mutation_data) / genes_with_data:.1f}"
                    )
        else:
            # Pairs mode breakdown
            if "is_alternative" in dataset.columns:  # Updated column name
                canonical_count = len(dataset[dataset["is_alternative"] == 0])
                alternative_count = len(dataset[dataset["is_alternative"] == 1])
                logger.info(f"\n  Sequence breakdown:")
                logger.info(f"  Canonical: {canonical_count:,}")
                logger.info(
                    f"  Alternative isoform: {alternative_count:,}"
                )  # Updated terminology

        # Length statistics
        logger.info(f"\n  Length statistics:")
        logger.info(f"  Average: {dataset['length'].mean():.1f} amino acids")
        logger.info(f"  Range: {dataset['length'].min()}-{dataset['length'].max()}")
        logger.info(f"  Median: {dataset['length'].median():.1f}")

    else:
        logger.info(f"\n  No sequences generated")

    # Failed genes details
    if failed_genes:
        logger.info(f"\n  Genes with errors:")
        for gene in failed_genes[:5]:  # Show first 5
            logger.info(
                f"  {gene}: No transcript-isoform pairs found"
            )  # Updated terminology
        if len(failed_genes) > 5:
            logger.info(f"  ... and {len(failed_genes) - 5} more")

    # Output files
    output_path = Path(output_dir)
    logger.info(f"\n  Results saved to: {output_dir}")

    if mutations_mode:
        fasta_file = output_path / "protein_sequences_with_mutations.fasta"
        csv_file = output_path / "protein_sequences_with_mutations.csv"
        logger.info(f"  Protein sequences with mutations saved to: {csv_file.name}")
        if fasta_file.exists():
            logger.info(f"  FASTA format saved to: {fasta_file.name}")
    else:
        fasta_file = output_path / "protein_sequences_pairs.fasta"
        csv_file = output_path / "protein_sequences_pairs.csv"
        logger.info(f"  Protein sequence pairs saved to: {csv_file.name}")
        if fasta_file.exists():
            logger.info(f"  FASTA format saved to: {fasta_file.name}")


def update_gencode_gene_names(
    input_gtf_path: Union[str, Path],
    output_gtf_path: Union[str, Path],
    reference_gtf_path: Union[str, Path],
    verbose: bool = True,
) -> Dict:
    """Update gene names in GENCODE GTF using an updated GENCODE reference GTF.

    Args:
        input_gtf_path: Path to input GENCODE GTF file to be updated
        output_gtf_path: Path to save updated GTF file
        reference_gtf_path: Path to reference GTF file with preferred gene names
        verbose: Print detailed update summary

    Returns:
        Dictionary containing update statistics
    """
    input_gtf_path, output_gtf_path = Path(input_gtf_path), Path(output_gtf_path)

    # Extract Ensembl IDs to gene name mappings from both files
    if verbose:
        logger.info(
            f"Creating gene ID to name mappings from reference GTF: {reference_gtf_path}"
        )

    # Dictionary to store updated gene names: key=Ensembl ID (without version), value=updated gene name
    gene_name_updates = {}

    # First, get the current gene names from GENCODE
    gencode_gene_names = {}
    with open(input_gtf_path, "r") as f:
        for line in f:
            if line.startswith("#") or len(line.strip()) == 0:
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            # Extract Ensembl ID and gene name
            attributes = fields[8]
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)

            if gene_id_match and gene_name_match:
                full_gene_id = gene_id_match.group(1)
                gene_name = gene_name_match.group(1)

                # Remove version from Ensembl ID (e.g., ENSG00000264520.1 -> ENSG00000264520)
                base_gene_id = full_gene_id.split(".")[0]
                gencode_gene_names[base_gene_id] = gene_name

    if verbose:
        logger.info(f"Extracted {len(gencode_gene_names)} gene names from GENCODE GTF")

    # Next, get the gene names from the reference GTF
    reference_gene_names = {}
    with open(reference_gtf_path, "r") as f:
        for line in f:
            if line.startswith("#") or len(line.strip()) == 0:
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            # Extract Ensembl ID and gene name
            attributes = fields[8]

            # Look for Ensembl ID pattern (with or without version)
            ensembl_id_match = re.search(r'gene_id "((ENSG\d+)(?:\.\d+)?)"', attributes)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)

            if ensembl_id_match and gene_name_match:
                base_gene_id = ensembl_id_match.group(
                    2
                )  # This gets the ID without version
                gene_name = gene_name_match.group(1)

                reference_gene_names[base_gene_id] = gene_name

    if verbose:
        logger.info(
            f"Extracted {len(reference_gene_names)} gene names from reference GTF"
        )

    # Create update mapping
    for ensembl_id, gencode_name in gencode_gene_names.items():
        if ensembl_id in reference_gene_names:
            reference_name = reference_gene_names[ensembl_id]

            # Only update if names differ
            if gencode_name != reference_name:
                gene_name_updates[ensembl_id] = reference_name

    if verbose:
        logger.info(f"Created {len(gene_name_updates)} gene name updates")

    # Process GENCODE GTF file to update gene names
    stats = {
        "total_lines": 0,
        "updated_lines": 0,
        "genes_processed": 0,
        "genes_updated": 0,
    }

    # Create output directory if it doesn't exist
    output_dir = output_gtf_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(input_gtf_path, "r") as f_in, open(output_gtf_path, "w") as f_out:
        for line in f_in:
            stats["total_lines"] += 1

            # Pass through comment lines
            if line.startswith("#"):
                f_out.write(line)
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                f_out.write(line)
                continue

            # Parse attributes
            attributes = fields[8]
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)

            if gene_id_match:
                full_gene_id = gene_id_match.group(1)
                base_gene_id = full_gene_id.split(".")[0]

                if fields[2] == "gene":
                    stats["genes_processed"] += 1

                # Check if we have an updated name for this gene
                if base_gene_id in gene_name_updates:
                    new_name = gene_name_updates[base_gene_id]

                    # Replace the gene_name attribute
                    new_attributes = re.sub(
                        r'gene_name "([^"]+)"', f'gene_name "{new_name}"', attributes
                    )

                    if new_attributes != attributes:
                        fields[8] = new_attributes
                        stats["updated_lines"] += 1

                        if fields[2] == "gene":
                            stats["genes_updated"] += 1

                        line = "\t".join(fields) + "\n"

            f_out.write(line)

    if verbose:
        logger.info(f"\nGTF Update Summary:")
        logger.info(f"  Total lines processed: {stats['total_lines']}")
        logger.info(f"  Genes processed: {stats['genes_processed']}")
        logger.info(f"  Genes with updated names: {stats['genes_updated']}")
        logger.info(f"  Total lines updated: {stats['updated_lines']}")
        logger.info(f"  Output saved to: {output_gtf_path}")

    return stats


def parse_gene_list(gene_list_path: Union[str, Path]) -> List[str]:
    """Parse a file containing a list of gene names.

    Args:
        gene_list_path (Union[str, Path]): Path to file containing gene names.

    Returns:
        List[str]: List of gene names.

    Raises:
        FileNotFoundError: If the gene list file does not exist.
    """
    gene_list_path = Path(gene_list_path)
    if not gene_list_path.exists():
        raise FileNotFoundError(f"Gene list file not found: {gene_list_path}")
    with open(gene_list_path, "r") as f:
        gene_names = [line.strip() for line in f if line.strip()]
    logger.info(f"Read {len(gene_names)} genes from {gene_list_path}")
    return gene_names


# NOTE: subset_gene_list() function removed - no longer needed for current workflow


def save_gene_level_results(
    results: List[Dict],
    output_dir: Union[str, Path],
    filename: str = "gene_level_results.csv",
) -> None:
    """Save gene-level analysis results to CSV file.

    Args:
        results: List of result dictionaries to save
        output_dir: Directory to save the results file
        filename: Name of the results file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Extract gene-level results, removing pair-specific data
    gene_results = []
    for result in results:
        # Skip pair-specific details to keep the summary concise
        if "pair_results" in result:
            result = {k: v for k, v in result.items() if k != "pair_results"}
        gene_results.append(result)
    # Save to CSV
    gene_df = pd.DataFrame(gene_results)
    output_path = output_dir / filename
    gene_df.to_csv(output_path, index=False)
    logger.info(f"Gene-level results saved to {output_path}")


def save_isoform_level_results(
    results: List[Dict],
    output_dir: Union[str, Path],
    filename: str = "isoform_level_results.csv",
) -> None:
    """Save detailed transcript-isoform pair analysis results to CSV file.

    Args:
        results: List of result dictionaries containing pair_results
        output_dir: Directory to save the results file
        filename: Name of the results file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract all pair results
    all_pairs = []
    for result in results:
        if "pair_results" in result and result["pair_results"]:
            gene_name = result["gene_name"]
            # Add gene name to each pair
            for pair in result["pair_results"]:
                pair_with_gene = {"gene_name": gene_name, **pair}
                all_pairs.append(pair_with_gene)

    # If no pairs found, return early
    if not all_pairs:
        logger.info(f"No transcript-isoform pairs found to save")
        return

    # Convert to DataFrame
    pairs_df = pd.DataFrame(all_pairs)

    # Rename columns for clarity
    column_renames = {
        "mutation_count_total": "total_mutations",
        "mutations_in_alt_start_site": "mutations_in_alt_start",
        "mutation_sources": "source_databases",
        "variant_ids": "all_variant_ids",
        "alt_start_site_variant_ids": "alt_start_variant_ids",
    }

    # Rename mutation count columns for clarity
    impact_type_renames = {}
    for col in pairs_df.columns:
        if col.startswith("mutations_"):
            # mutations_missense_variant -> count_missense
            impact_type = col.replace("mutations_", "")
            new_name = f"count_{impact_type}"
            impact_type_renames[col] = new_name
        elif col.startswith("variant_ids_"):
            # variant_ids_missense_variant -> ids_missense
            impact_type = col.replace("variant_ids_", "")
            new_name = f"ids_{impact_type}"
            impact_type_renames[col] = new_name

    column_renames.update(impact_type_renames)
    pairs_df = pairs_df.rename(columns=column_renames)

    # Organize columns in 5-section structure with count+IDs pairing
    column_order = []

    # SECTION 1: Feature Identification (8 columns)
    section1 = [
        "gene_name",
        "transcript_id",
        "feature_id",
        "feature_type",
        "feature_start",
        "feature_end",
        "total_mutations",
        "source_databases",
    ]
    column_order.extend(section1)

    # SECTION 2: Summary Counts ONLY (no IDs - those are redundant)
    # Order: missense, nonsense, frameshift, inframe_deletion, inframe_insertion, synonymous, mutations_in_alt_start
    section2 = [
        "count_missense_variant",
        "count_nonsense_variant",
        "count_frameshift_variant",
        "count_inframe_deletion",
        "count_inframe_insertion",
        "count_synonymous_variant",
        "mutations_in_alt_start",
    ]
    column_order.extend(section2)

    # SECTION 3: Per-Source Totals (3 columns)
    section3 = [
        "count_clinvar",
        "count_gnomad",
        "count_cosmic",
    ]
    column_order.extend(section3)

    # SECTION 4: Source×Impact Matrix PAIRED (count immediately followed by IDs)
    # Order by source (clinvar, gnomad, cosmic), then by impact (missense, nonsense, frameshift, inframe_deletion, inframe_insertion, synonymous)
    sources = ["clinvar", "gnomad", "cosmic"]
    impacts = [
        "missense_variant",
        "nonsense_variant",
        "frameshift_variant",
        "inframe_deletion",
        "inframe_insertion",
        "synonymous_variant",
    ]

    section4 = []
    for source in sources:
        for impact in impacts:
            count_col = f"count_{source}_{impact}"
            ids_col = f"ids_{source}_{impact}"
            # Add both count and IDs in sequence
            section4.extend([count_col, ids_col])
    column_order.extend(section4)

    # SECTION 5: Special Columns (2 columns)
    section5 = [
        "alt_start_variant_ids",
        "all_variant_ids",
    ]
    column_order.extend(section5)

    # Reorder columns that exist in the dataframe
    available_columns = [col for col in column_order if col in pairs_df.columns]

    # Add any remaining columns that weren't explicitly ordered
    remaining_columns = [
        col for col in pairs_df.columns if col not in available_columns
    ]
    final_column_order = available_columns + remaining_columns

    # Apply the column order and save
    if final_column_order:  # Only reorder if we have columns to reorder
        pairs_df = pairs_df[final_column_order]

    # Save to CSV
    output_path = output_dir / filename
    pairs_df.to_csv(output_path, index=False)
    logger.info(f"Transcript-isoform pair analysis saved to {output_path}")


def print_mutation_summary(results_df, output_dir):
    """Print summary statistics from the analysis results.

    Args:
        results_df: DataFrame containing analysis results
        output_dir: Directory where output files are saved
    """
    logger.info("\nAnalysis Summary:")
    logger.info(f"  Total genes processed: {len(results_df)}")
    logger.info("\n  Status breakdown:")
    for status, count in results_df["status"].value_counts().items():
        logger.info(f"  {status}: {count}")

    # Transcript-isoform statistics (updated for new structure)
    successful_genes = results_df[results_df["status"] == "success"]
    if not successful_genes.empty:
        total_transcripts = successful_genes["total_transcripts"].sum()
        total_features = successful_genes["alternative_features"].sum()
        total_pairs = successful_genes["transcript_feature_pairs"].sum()

        # Calculate average pairs per gene
        avg_pairs_per_gene = (
            total_pairs / len(successful_genes) if len(successful_genes) > 0 else 0
        )

        logger.info(f"\n  Transcript-Isoform Analysis:")
        logger.info(f"  Total transcripts across all genes: {total_transcripts}")
        logger.info(f"  Total alternative isoform features: {total_features}")
        logger.info(f"  Total transcript-isoform pairs: {total_pairs}")
        logger.info(f"  Average pairs per gene: {avg_pairs_per_gene:.2f}")

        # Mutation statistics if available
        if "mutations_filtered" in successful_genes.columns:
            total_mutations = successful_genes["mutations_filtered"].sum()
            logger.info(f"\n  Mutation Analysis:")
            logger.info(
                f"  Total mutations in alternative isoform regions: {total_mutations}"
            )

            # Try to load mutation analysis results to report statistics
            mutation_analysis_path = (
                Path(output_dir) / "isoform_level_results.csv"
            )  # Updated filename
            if mutation_analysis_path.exists():
                try:
                    pairs_df = pd.read_csv(mutation_analysis_path)

                    # Get total mutations across all pairs
                    total_pair_mutations = (
                        pairs_df["total_mutations"].sum()
                        if "total_mutations" in pairs_df.columns
                        else 0
                    )
                    avg_mutations_per_pair = (
                        total_pair_mutations / len(pairs_df) if len(pairs_df) > 0 else 0
                    )
                    logger.info(
                        f"  Total mutations across all pairs: {total_pair_mutations}"
                    )
                    logger.info(
                        f"  Average mutations per transcript-isoform pair: {avg_mutations_per_pair:.2f}"
                    )

                    # Show breakdown by isoform type if available
                    if "feature_type" in pairs_df.columns:
                        logger.info(f"")
                        logger.info(f"  Breakdown by isoform type:")

                        type_mutations = pairs_df.groupby("feature_type")[
                            "total_mutations"
                        ].sum()
                        for isoform_type, count in type_mutations.items():
                            type_pairs = len(
                                pairs_df[pairs_df["feature_type"] == isoform_type]
                            )
                            avg_per_type = count / type_pairs if type_pairs > 0 else 0
                            logger.info(
                                f"  {isoform_type.capitalize()}: {count} total ({avg_per_type:.1f} avg/pair)"
                            )

                    # Print statistics for each mutation category
                    mutation_categories = [
                        col for col in pairs_df.columns if col.startswith("mutations_")
                    ]
                    if mutation_categories:
                        logger.info(f"")
                        logger.info(f"  Breakdown by mutation category:")

                        for category in mutation_categories:
                            # Skip if the column doesn't exist in the dataframe
                            if category not in pairs_df.columns:
                                continue

                            # Convert category name from mutations_missense_variant to "Missense Variant"
                            category_name = category.replace("mutations_", "").replace(
                                "_", " "
                            )
                            category_name = category_name.title()

                            # Calculate statistics for this category
                            category_total = pairs_df[category].sum()
                            category_percent = (
                                (category_total / total_pair_mutations * 100)
                                if total_pair_mutations > 0
                                else 0
                            )

                            logger.info(
                                f"  {category_name}: {category_total} ({category_percent:.1f}%)"
                            )

                    logger.info(
                        f"  Detailed results available in isoform_level_results.csv"  # Updated filename
                    )
                except Exception as e:
                    logger.error(f"Error reading mutation analysis: {str(e)}")
                    logger.info(f"  Error reading detailed mutation analysis: {str(e)}")

    # Genes with errors
    error_genes = results_df[results_df["status"] == "error"]
    if not error_genes.empty:
        logger.info("\n  Genes with errors:")
        for _, row in error_genes.iterrows():
            logger.info(f"  {row['gene_name']}: {row['error']}")

    logger.info(f"\n  Results saved to: {output_dir}")
    logger.info(
        f"  Gene-level results saved to: {Path(output_dir) / 'gene_level_results.csv'}"
    )
    logger.info(
        f"  Detailed mutation analysis by pair saved to: {Path(output_dir) / 'isoform_level_results.csv'}"  # Updated filename
    )


def print_translation_summary(
    dataset: pd.DataFrame,
    successful_genes: int,
    total_genes: int,
    failed_genes: List[str],
    mutations_mode: bool,
    output_dir: str,
) -> None:
    """Print comprehensive summary of protein sequence generation results.

    Args:
        dataset: Generated protein sequence dataset
        successful_genes: Number of successfully processed genes
        total_genes: Total number of genes attempted
        failed_genes: List of gene names that failed processing
        mutations_mode: Whether mutations were included
        output_dir: Output directory path
    """
    logger.info(f"\nProtein Sequence Generation Summary:")
    logger.info(f"  Total genes processed: {total_genes}")

    # Status breakdown
    logger.info(f"\n  Status breakdown:")
    logger.info(f"  Success: {successful_genes}")
    if failed_genes:
        logger.info(f"  Failed: {len(failed_genes)}")
    else:
        logger.info(f"  Failed: 0")

    # Dataset statistics
    if not dataset.empty:
        logger.info(f"\n  Sequence Generation:")
        logger.info(f"  Total sequences generated: {len(dataset):,}")

        # Calculate transcript-isoform pairs
        genes_with_data = dataset["gene"].nunique()
        if mutations_mode and "variant_type" in dataset.columns:
            # Count unique transcript-isoform pairs (canonical + alternative base pairs)
            base_sequences = dataset[
                dataset["variant_type"].isin(
                    ["canonical", "alternative"]
                )  # Updated terminology
            ]
            unique_pairs = (
                len(base_sequences) // 2
                if len(base_sequences) % 2 == 0
                else (len(base_sequences) + 1) // 2
            )
        else:
            # For pairs mode, total sequences / 2 = pairs
            unique_pairs = (
                len(dataset) // 2 if len(dataset) % 2 == 0 else (len(dataset) + 1) // 2
            )

        logger.info(
            f"  Transcript-isoform pairs: {unique_pairs}"
        )  # Updated terminology
        logger.info(
            f"  Average sequences per gene: {len(dataset) / genes_with_data:.1f}"
        )

        # Mode-specific breakdown
        if mutations_mode and "variant_type" in dataset.columns:
            logger.info(f"\n  Sequence breakdown:")
            type_counts = dataset["variant_type"].value_counts()
            for variant_type, count in type_counts.items():
                percentage = (count / len(dataset)) * 100
                logger.info(f"  {variant_type}: {count:,} ({percentage:.1f}%)")

            # Mutation-specific statistics
            if "canonical_mutated" in type_counts:
                mutation_data = dataset[dataset["variant_type"] == "canonical_mutated"]
                if not mutation_data.empty:
                    logger.info(f"\n  Mutation Analysis:")
                    logger.info(f"  Total mutations integrated: {len(mutation_data)}")
                    logger.info(
                        f"  Unique mutation positions: {mutation_data['mutation_position'].nunique()}"
                    )

                    if "aa_change" in mutation_data.columns:
                        aa_changes = mutation_data["aa_change"].dropna()
                        silent_mutations = len(mutation_data) - len(aa_changes)
                        logger.info(f"  Mutations with AA changes: {len(aa_changes)}")
                        if silent_mutations > 0:
                            logger.info(f"  Silent mutations: {silent_mutations}")

                    # Breakdown by impact type if available
                    if "mutation_impact" in mutation_data.columns:
                        impact_counts = mutation_data["mutation_impact"].value_counts()
                        logger.info(f"")
                        logger.info(f"  Breakdown by impact type:")
                        for impact_type, count in impact_counts.items():
                            percentage = (count / len(mutation_data)) * 100
                            logger.info(f"  {impact_type}: {count} ({percentage:.1f}%)")

                    logger.info(
                        f"  Average mutations per gene: {len(mutation_data) / genes_with_data:.1f}"
                    )
        else:
            # Pairs mode breakdown
            if "is_alternative" in dataset.columns:  # Updated column name
                canonical_count = len(dataset[dataset["is_alternative"] == 0])
                alternative_count = len(dataset[dataset["is_alternative"] == 1])
                logger.info(f"\n  Sequence breakdown:")
                logger.info(f"  Canonical: {canonical_count:,}")
                logger.info(
                    f"  Alternative isoform: {alternative_count:,}"
                )  # Updated terminology

        # Length statistics
        logger.info(f"\n  Length statistics:")
        logger.info(f"  Average: {dataset['length'].mean():.1f} amino acids")
        logger.info(f"  Range: {dataset['length'].min()}-{dataset['length'].max()}")
        logger.info(f"  Median: {dataset['length'].median():.1f}")

    else:
        logger.info(f"\n  No sequences generated")

    # Failed genes details
    if failed_genes:
        logger.info(f"\n  Genes with errors:")
        for gene in failed_genes[:5]:  # Show first 5
            logger.info(
                f"  {gene}: No transcript-isoform pairs found"
            )  # Updated terminology
        if len(failed_genes) > 5:
            logger.info(f"  ... and {len(failed_genes) - 5} more")

    # Output files
    output_path = Path(output_dir)
    logger.info(f"\n  Results saved to: {output_dir}")

    if mutations_mode:
        fasta_file = output_path / "protein_sequences_with_mutations.fasta"
        csv_file = output_path / "protein_sequences_with_mutations.csv"
        logger.info(f"  Protein sequences with mutations saved to: {csv_file.name}")
        if fasta_file.exists():
            logger.info(f"  FASTA format saved to: {fasta_file.name}")
    else:
        fasta_file = output_path / "protein_sequences_pairs.fasta"
        csv_file = output_path / "protein_sequences_pairs.csv"
        logger.info(f"  Protein sequence pairs saved to: {csv_file.name}")
        if fasta_file.exists():
            logger.info(f"  FASTA format saved to: {fasta_file.name}")


def load_pre_validated_variants(mutations_file: str) -> Dict[str, Set[str]]:
    """Load pre-validated variant IDs from step 2 results.

    Args:
        mutations_file (str): Path to isoform_level_results.csv

    Returns:
        Dict[str, Set[str]]: Dictionary mapping gene_name -> set of variant IDs
    """
    logger.info(f"Loading pre-validated variant IDs from {mutations_file}")

    # Read the mutation results
    mutations_df = pd.read_csv(mutations_file)

    if mutations_df.empty:
        logger.warning("No mutation records found in results file")
        return {}

    logger.info(f"Found {len(mutations_df)} mutation records")

    # Build list of source×impact ID columns dynamically
    # After Task 2, we use source-specific columns instead of summary columns
    variant_id_columns = []
    sources = ["clinvar", "gnomad", "cosmic"]
    impacts = [
        "missense_variant",
        "nonsense_variant",
        "frameshift_variant",
        "inframe_deletion",
        "inframe_insertion",
        "synonymous_variant",
    ]

    for source in sources:
        for impact in impacts:
            variant_id_columns.append(f"ids_{source}_{impact}")

    # Extract variant IDs by gene
    pre_validated_variants = {}
    total_variant_ids = 0

    for _, row in mutations_df.iterrows():
        gene_name = row.get("gene_name", "")
        if not gene_name:
            continue

        if gene_name not in pre_validated_variants:
            pre_validated_variants[gene_name] = set()

        # Collect variant IDs from all impact type columns
        for col in variant_id_columns:
            if col in row and pd.notna(row[col]) and str(row[col]).strip():
                # Parse comma-separated variant IDs
                variant_ids = [
                    vid.strip() for vid in str(row[col]).split(",") if vid.strip()
                ]
                pre_validated_variants[gene_name].update(variant_ids)
                total_variant_ids += len(variant_ids)

    genes_with_variants = len([g for g in pre_validated_variants.values() if g])
    logger.info(
        f"Loaded {total_variant_ids} pre-validated variant IDs across {genes_with_variants} genes"
    )

    # Print summary by gene (top 10)
    gene_counts = [
        (gene, len(variants))
        for gene, variants in pre_validated_variants.items()
        if variants
    ]
    gene_counts.sort(key=lambda x: x[1], reverse=True)

    if gene_counts:
        logger.info("Top genes by variant count:")
        for gene, count in gene_counts[:10]:
            logger.info(f"  {gene}: {count} variants")
        if len(gene_counts) > 10:
            logger.info(f"  ... and {len(gene_counts) - 10} more genes")

    return pre_validated_variants


def simple_transcript_based_cleanup(
    input_bed_path: Union[str, Path],
    gtf_path: Union[str, Path],
    output_bed_path: Union[str, Path],
    refseq_gtf_path: Optional[Union[str, Path]] = None,
    categories_to_keep: Set[str] = {"Annotated", "Extended", "Truncated", "uORF"},
    verbose: bool = True,
) -> Dict:
    """Simplified cleanup for BED files that already have transcript IDs assigned.

    This function is designed for ribosome profiling data where transcript IDs are already
    in the BED file name field (format: CODON_CATEGORY_GENE_TRANSCRIPT).

    Steps:
    1. Parse BED entries and extract transcript ID from name field
    2. Load GTF for transcript → gene mapping
    3. Load RefSeq mapping
    4. Validate entries and map to RefSeq IDs
    5. Add canonical starts where missing (except for uORF-only transcripts)
    6. Remove Extended/Truncated entries for transcripts without annotated starts
    7. Validate and deduplicate Annotated starts (max 1 per transcript)
    8. Generate statistics including multi-TIS per transcript/category

    Args:
        input_bed_path: Path to input BED file with transcript IDs in name field
        gtf_path: Path to GTF annotation file
        output_bed_path: Path to save enhanced 8-column BED file
        refseq_gtf_path: Path to RefSeq GTF for mapping (optional)
        categories_to_keep: Set of categories to include in output
        verbose: Print detailed progress information

    Returns:
        Dictionary with comprehensive processing statistics including:
        - total_entries: Total input entries
        - filtered_by_category: Entries filtered out
        - valid_entries: Entries kept after filtering
        - ensembl_mapped: Entries with Ensembl transcript IDs
        - refseq_mapped: Entries with RefSeq IDs
        - canonical_starts_added: Number of canonical starts added from GTF
        - extended_truncated_removed: Number of Extended/Truncated entries removed (no annotated start)
        - transcripts_with_removed_entries: Number of transcripts affected by removal
        - duplicate_annotated_removed: Number of duplicate Annotated entries removed
        - multi_tis_transcripts: Dict of {category: count} for transcripts with multiple TIS
        - total_transcripts_with_multi_tis: Number of transcripts with any multi-TIS
        - final_gene_count: Number of unique genes
    """
    if verbose:
        logger.info(f"\n{'=' * 60}")
        logger.info("Simple Transcript-Based BED Cleanup")
        logger.info(f"{'=' * 60}")
        logger.info(f"  Input BED: {input_bed_path}")
        logger.info(f"  GTF: {gtf_path}")
        logger.info(f"  Output BED: {output_bed_path}")

    # Step 1: Load GTF transcript information
    if verbose:
        logger.info("\nStep 1: Loading GTF transcript information...")

    transcript_info = {}  # transcript_id → {gene_id, gene_name, chrom, strand, start_pos}
    canonical_starts = {}  # transcript_id → {chrom, start, end, strand}

    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feature = fields[2]
            attrs = fields[8]

            # Extract transcript info
            transcript_match = re.search(r'transcript_id "([^"]+)"', attrs)
            if not transcript_match:
                continue

            transcript_id = transcript_match.group(1)

            # Get transcript-level info
            if feature == "transcript":
                gene_match = re.search(r'gene_id "([^"]+)"', attrs)
                gene_name_match = re.search(r'gene_name "([^"]+)"', attrs)

                if gene_match and gene_name_match:
                    transcript_info[transcript_id] = {
                        "gene_id": gene_match.group(1),
                        "gene_name": gene_name_match.group(1),
                        "chrom": fields[0],
                        "strand": fields[6],
                        "start_pos": int(fields[3])
                        if fields[6] == "+"
                        else int(fields[4]),
                    }

            # Get canonical start codon positions
            elif feature == "start_codon" and transcript_id not in canonical_starts:
                canonical_starts[transcript_id] = {
                    "chrom": fields[0],
                    "start": int(fields[3]),  # Keep 1-based to match input BED format
                    "end": int(fields[4]),
                    "strand": fields[6],
                }

    if verbose:
        logger.info(f"    ✓ Loaded {len(transcript_info)} transcripts")
        logger.info(f"    ✓ Loaded {len(canonical_starts)} canonical start sites")

    # Step 2: Load RefSeq mapping (with gene name fallback)
    if verbose:
        logger.info("\nStep 2: Loading RefSeq mapping...")

    refseq_mapping = {}  # ensembl_base_id → refseq_id
    refseq_gene_mapping = {}  # gene_name → [refseq_ids]  # For fallback lookups

    if refseq_gtf_path and Path(refseq_gtf_path).exists():
        with open(refseq_gtf_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 9 or fields[2] != "transcript":
                    continue

                attrs = fields[8]

                # Look for Ensembl cross-reference and RefSeq ID
                ensembl_match = re.search(r'db_xref "Ensembl:([^"]+)"', attrs)
                refseq_match = re.search(r'transcript_id "([^"]+)"', attrs)
                gene_name_match = re.search(r'gene "([^"]+)"', attrs)

                refseq_id = refseq_match.group(1) if refseq_match else None

                # Only process RefSeq transcript IDs
                if refseq_id and refseq_id.startswith(("NM_", "NR_", "XM_", "XR_")):
                    # Primary mapping: Ensembl ID → RefSeq ID
                    if ensembl_match:
                        ensembl_id = ensembl_match.group(1)
                        base_ensembl = ensembl_id.split(".")[0]
                        refseq_mapping[base_ensembl] = refseq_id
                        refseq_mapping[ensembl_id] = refseq_id  # Also store versioned

                    # Secondary mapping: Gene name → RefSeq IDs (for fallback)
                    if gene_name_match:
                        gene_name = gene_name_match.group(1)
                        if gene_name not in refseq_gene_mapping:
                            refseq_gene_mapping[gene_name] = []
                        if refseq_id not in refseq_gene_mapping[gene_name]:
                            refseq_gene_mapping[gene_name].append(refseq_id)

        if verbose:
            logger.info(f"    ✓ Loaded {len(refseq_mapping)} Ensembl→RefSeq mappings")
            logger.info(
                f"    ✓ Loaded {len(refseq_gene_mapping)} gene name→RefSeq mappings (for fallback)"
            )
    else:
        if verbose:
            logger.info(f"    ⚠ RefSeq GTF not found, RefSeq mapping will be limited")

    # Step 3: Parse BED file
    if verbose:
        logger.info("\nStep 3: Parsing BED file...")

    stats = defaultdict(int)
    entries = []
    transcript_category_counts = defaultdict(
        lambda: defaultdict(int)
    )  # transcript → category → count

    with open(input_bed_path, "r") as f:
        for line in f:
            if not line.strip():
                continue

            stats["total_entries"] += 1
            fields = line.strip().split("\t")

            if len(fields) < 5:
                stats["invalid_format"] += 1
                continue

            chrom, start, end, name = fields[:4]
            # Input BED has 7 columns: chr, start, end, name, strand, dna_sequence, protein_sequence:strand
            # Input BED has no score column, so always use '0'
            score = "0"
            strand = (
                fields[6].split(":")[-1]
                if len(fields) > 6 and ":" in fields[6]
                else (fields[4] if len(fields) > 4 else ".")
            )

            # Extract sequences from input
            dna_seq = fields[5] if len(fields) > 5 else ""
            # Protein sequence is in column 7, format: "PROTEIN_SEQUENCE:STRAND"
            protein_seq = (
                fields[6].rsplit(":", 1)[0]
                if len(fields) > 6 and ":" in fields[6]
                else ""
            )

            # Parse name: CODON_CATEGORY_GENE_TRANSCRIPT
            parts = name.split("_")
            if len(parts) < 4:
                stats["invalid_name_format"] += 1
                continue

            codon = parts[0]
            category = parts[1]
            gene = parts[2]  # This will be replaced with GTF gene name later
            transcript = parts[3]

            # Filter by category
            if category not in categories_to_keep:
                stats["filtered_by_category"] += 1
                continue

            # Validate transcript exists in GTF
            if transcript not in transcript_info:
                stats["transcript_not_in_gtf"] += 1
                continue

            stats["valid_entries"] += 1

            # Track transcript-category counts for multi-TIS detection
            transcript_category_counts[transcript][category] += 1

            entries.append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "name": name,
                    "score": score,
                    "strand": strand,
                    "transcript": transcript,
                    "gene": gene,  # Original gene from input; will use GTF name in output
                    "category": category,
                    "codon": codon,  # Store codon for output formatting
                    "dna_sequence": dna_seq,
                    "protein_sequence": protein_seq,
                }
            )

    if verbose:
        logger.info(f"    ✓ Total entries: {stats['total_entries']}")
        logger.info(f"    ✓ Valid entries: {stats['valid_entries']}")
        logger.info(f"    ✓ Filtered by category: {stats['filtered_by_category']}")

    # Step 4: Add missing canonical starts
    if verbose:
        logger.info("\nStep 4: Adding missing canonical starts...")

    # Find transcripts with alternatives but no annotated start
    transcripts_with_annotated = {
        e["transcript"] for e in entries if e["category"] == "Annotated"
    }
    transcripts_without_annotated = set()

    for entry in entries:
        if entry["transcript"] not in transcripts_with_annotated:
            transcripts_without_annotated.add(entry["transcript"])

    added_canonical = 0
    added_canonical_transcripts = []  # Track which transcripts had canonical starts added

    for transcript_id in transcripts_without_annotated:
        if transcript_id not in canonical_starts:
            continue

        # Check if this transcript only has uORFs - skip adding canonical start for uORF-only transcripts
        transcript_entries = [e for e in entries if e["transcript"] == transcript_id]
        transcript_categories = {e["category"] for e in transcript_entries}
        if transcript_categories == {"uORF"}:
            continue

        canonical = canonical_starts[transcript_id]
        trans_info = transcript_info[transcript_id]

        # Check if this position already exists
        existing = any(
            e["start"] == canonical["start"] and e["transcript"] == transcript_id
            for e in entries
        )

        if not existing:
            entries.append(
                {
                    "chrom": canonical["chrom"],
                    "start": canonical["start"],
                    "end": canonical["end"],
                    "name": f"{trans_info['gene_name']}_{trans_info['gene_id']}_Annotated_ATG",
                    "score": "0",
                    "strand": canonical["strand"],
                    "transcript": transcript_id,
                    "gene": trans_info["gene_name"],
                    "category": "Annotated",
                    "codon": "ATG",  # Canonical starts are always ATG
                    "dna_sequence": "",  # No sequence for added canonical starts
                    "protein_sequence": "",
                }
            )
            added_canonical += 1
            added_canonical_transcripts.append(transcript_id)
            transcript_category_counts[transcript_id]["Annotated"] += 1

    stats["canonical_starts_added"] = added_canonical
    stats["canonical_starts_added_transcripts"] = added_canonical_transcripts[
        :10
    ]  # Store up to 10 examples

    if verbose:
        logger.info(
            f"    ✓ Added {added_canonical} canonical starts to transcripts with alternatives"
        )
        if added_canonical > 0 and added_canonical_transcripts:
            logger.info(f"    Examples: {', '.join(added_canonical_transcripts[:3])}")

    # Step 4b: Remove Extended/Truncated entries without annotated starts
    if verbose:
        logger.info(
            "\nStep 4b: Removing Extended/Truncated entries without annotated starts..."
        )

    # Re-find transcripts with annotated starts (now including added ones)
    transcripts_with_annotated = {
        e["transcript"] for e in entries if e["category"] == "Annotated"
    }

    # Filter out Extended/Truncated entries for transcripts without annotated starts
    filtered_entries = []
    removed_extended_truncated = 0
    removed_transcripts = set()

    for entry in entries:
        # Keep all entries that are not Extended or Truncated
        if entry["category"] not in {"Extended", "Truncated"}:
            filtered_entries.append(entry)
        # For Extended/Truncated, only keep if transcript has an annotated start
        elif entry["transcript"] in transcripts_with_annotated:
            filtered_entries.append(entry)
        else:
            # Remove this Extended/Truncated entry
            removed_extended_truncated += 1
            removed_transcripts.add(entry["transcript"])
            # Update transcript_category_counts
            if entry["transcript"] in transcript_category_counts:
                if entry["category"] in transcript_category_counts[entry["transcript"]]:
                    transcript_category_counts[entry["transcript"]][
                        entry["category"]
                    ] -= 1

    entries = filtered_entries

    stats["extended_truncated_removed"] = removed_extended_truncated
    stats["transcripts_with_removed_entries"] = len(removed_transcripts)
    stats["removed_transcript_examples"] = list(removed_transcripts)[
        :10
    ]  # Store up to 10 examples

    if verbose:
        logger.info(
            f"    ✓ Removed {removed_extended_truncated} Extended/Truncated entries from {len(removed_transcripts)} transcripts without annotated starts"
        )
        if removed_extended_truncated > 0 and removed_transcripts:
            logger.info(f"    Examples: {', '.join(list(removed_transcripts)[:3])}")

    # Step 4c: Validate and deduplicate Annotated starts per transcript
    if verbose:
        logger.info("\nStep 4c: Validating Annotated starts (max 1 per transcript)...")

    # Group entries by transcript and category
    transcript_annotated_map = defaultdict(list)
    for i, entry in enumerate(entries):
        if entry["category"] == "Annotated":
            transcript_annotated_map[entry["transcript"]].append((i, entry))

    # Find transcripts with multiple Annotated starts and keep only the valid one
    entries_to_remove = set()
    duplicate_annotated_removed = 0
    transcripts_with_duplicates = []  # Track affected transcripts

    for transcript_id, annotated_entries in transcript_annotated_map.items():
        if len(annotated_entries) > 1:
            transcripts_with_duplicates.append(transcript_id)

            # Keep the one that matches the canonical start position from GTF
            if transcript_id in canonical_starts:
                canonical_pos = canonical_starts[transcript_id]["start"]
                valid_entry_idx = None
                kept_entry_start = None

                for idx, entry in annotated_entries:
                    if entry["start"] == canonical_pos:
                        valid_entry_idx = idx
                        kept_entry_start = canonical_pos
                        break

                # Only mark entries for removal if we found a valid match
                if valid_entry_idx is not None:
                    # Mark all others for removal
                    for idx, entry in annotated_entries:
                        if idx != valid_entry_idx:
                            entries_to_remove.add(idx)
                            duplicate_annotated_removed += 1
                else:
                    # No entry matched canonical position - keep all and log warning
                    if verbose:
                        logger.warning(
                            f"    ⚠ Transcript {transcript_id}: No Annotated entry matches GTF canonical position {canonical_pos}"
                        )
                        logger.warning(
                            f"       Found positions: {[e[1]['start'] for e in annotated_entries]}"
                        )
                        logger.warning(
                            f"       Keeping all {len(annotated_entries)} duplicate Annotated entries"
                        )
            else:
                # No canonical start in GTF - cannot validate, log warning
                if verbose:
                    logger.warning(
                        f"    ⚠ Transcript {transcript_id}: No GTF canonical start available for validation"
                    )
                    logger.warning(
                        f"       Keeping all {len(annotated_entries)} duplicate Annotated entries"
                    )

    # Remove marked entries and update transcript_category_counts
    if entries_to_remove:
        entries = [e for i, e in enumerate(entries) if i not in entries_to_remove]
        stats["duplicate_annotated_removed"] = duplicate_annotated_removed
        stats["transcripts_with_duplicate_annotated"] = len(transcripts_with_duplicates)
        stats["duplicate_annotated_examples"] = transcripts_with_duplicates[
            :10
        ]  # Store up to 10 examples

        # Update transcript_category_counts to reflect removed entries
        # This ensures multi-TIS statistics are calculated from FINAL data, not intermediate
        for transcript_id, annotated_entries in transcript_annotated_map.items():
            if len(annotated_entries) > 1:
                # Decrement count by number of removed entries for this transcript
                removed_for_transcript = sum(
                    1 for idx, _ in annotated_entries if idx in entries_to_remove
                )
                if removed_for_transcript > 0:
                    transcript_category_counts[transcript_id]["Annotated"] -= (
                        removed_for_transcript
                    )

        if verbose:
            logger.info(
                f"    ✓ Removed {duplicate_annotated_removed} duplicate Annotated entries from {len(transcripts_with_duplicates)} transcripts"
            )
            if transcripts_with_duplicates:
                logger.info(
                    f"    Examples: {', '.join(transcripts_with_duplicates[:3])}"
                )
    else:
        stats["duplicate_annotated_removed"] = 0
        stats["transcripts_with_duplicate_annotated"] = 0
        stats["duplicate_annotated_examples"] = []
        if verbose:
            logger.info(f"    ✓ No duplicate Annotated entries found")

    # Step 5: Map to Ensembl and RefSeq IDs
    if verbose:
        logger.info("\nStep 5: Mapping to Ensembl and RefSeq IDs...")

    output_entries = []
    refseq_gene_name_matches = 0  # Track fallback successes

    for entry in entries:
        transcript_id = entry["transcript"]
        base_transcript = transcript_id.split(".")[0]

        # Ensembl ID is already in the entry
        ensembl_id = transcript_id
        stats["ensembl_mapped"] += 1

        # Get gene name from GTF (authoritative source, with v47 names)
        gene_name_from_gtf = (
            transcript_info[transcript_id]["gene_name"]
            if transcript_id in transcript_info
            else entry["gene"]
        )

        # Primary lookup: Direct Ensembl ID → RefSeq ID mapping
        refseq_id = refseq_mapping.get(
            base_transcript, refseq_mapping.get(transcript_id, None)
        )

        # Fallback lookup: Try gene name if direct lookup failed
        if refseq_id is None and gene_name_from_gtf in refseq_gene_mapping:
            # Use first RefSeq ID for this gene name
            refseq_candidates = refseq_gene_mapping[gene_name_from_gtf]
            if refseq_candidates:
                refseq_id = refseq_candidates[0]  # Use first match
                refseq_gene_name_matches += 1

        if refseq_id is None:
            refseq_id = "NA"

        if refseq_id != "NA":
            stats["refseq_mapped"] += 1

        output_entries.append(
            {
                **entry,
                "ensembl_id": ensembl_id,
                "refseq_id": refseq_id,
                "gene_name_gtf": gene_name_from_gtf,  # Add GTF gene name for output formatting
            }
        )

    stats["refseq_gene_name_matches"] = refseq_gene_name_matches

    if verbose:
        logger.info(
            f"    ✓ Ensembl mapped: {stats['ensembl_mapped']}/{len(output_entries)}"
        )
        logger.info(
            f"    ✓ RefSeq mapped: {stats['refseq_mapped']}/{len(output_entries)}"
        )
        if refseq_gene_name_matches > 0:
            logger.info(
                f"    ✓ RefSeq via gene name fallback: {refseq_gene_name_matches}"
            )

        # Show some examples of RefSeq mappings for debugging
        if stats["refseq_mapped"] > 0:
            refseq_examples = [e for e in output_entries if e["refseq_id"] != "NA"][:3]
            if refseq_examples:
                logger.info(f"    RefSeq mapping examples:")
                for ex in refseq_examples:
                    logger.info(f"      {ex['transcript']} → {ex['refseq_id']}")

    # Step 6: Calculate multi-TIS statistics and collect examples
    multi_tis_transcripts = defaultdict(int)
    multi_tis_examples = defaultdict(list)  # Store example transcript IDs per category
    total_multi_tis_transcripts = set()

    for transcript, category_counts in transcript_category_counts.items():
        for category, count in category_counts.items():
            if count > 1:
                multi_tis_transcripts[category] += 1
                total_multi_tis_transcripts.add(transcript)
                # Collect up to 3 examples per category
                if len(multi_tis_examples[category]) < 3:
                    multi_tis_examples[category].append(transcript)

    stats["multi_tis_transcripts"] = dict(multi_tis_transcripts)
    stats["multi_tis_examples"] = {
        k: v for k, v in multi_tis_examples.items()
    }  # Convert to regular dict
    stats["total_transcripts_with_multi_tis"] = len(total_multi_tis_transcripts)

    # Count unique genes
    unique_genes = {
        transcript_info[e["transcript"]]["gene_name"]
        for e in output_entries
        if e["transcript"] in transcript_info
    }
    stats["final_gene_count"] = len(unique_genes)
    stats["final_entries"] = len(
        output_entries
    )  # Total output entries including added canonical starts

    # Step 7: Write output
    if verbose:
        logger.info("\nStep 6: Writing output...")

    output_path = Path(output_bed_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        for entry in sorted(output_entries, key=lambda e: (e["chrom"], e["start"])):
            # Get gene info from GTF (authoritative source)
            gene_name_from_gtf = entry["gene_name_gtf"]
            gene_id = (
                transcript_info[entry["transcript"]]["gene_id"]
                if entry["transcript"] in transcript_info
                else "UNKNOWN"
            )

            # Reformat name to match AlternativeIsoform expected format: GENENAME_GENEID_CATEGORY_CODON
            # Use GTF gene name, not input BED gene name
            # Note: Score removed from name (input BED has no score column)
            # Use stored category and codon from entry dict instead of parsing name
            new_name = (
                f"{gene_name_from_gtf}_{gene_id}_{entry['category']}_{entry['codon']}"
            )

            # Write output in BED6+ format: standard BED6 columns first, then extra columns
            # BED6: chr, start, end, name, score, strand
            # Extra: ensembl_id, refseq_id, dna_sequence, protein_sequence
            # Use '.' for missing/empty sequences (BED standard)
            dna_seq = entry["dna_sequence"] if entry["dna_sequence"] else "."
            protein_seq = (
                entry["protein_sequence"] if entry["protein_sequence"] else "."
            )

            f.write(f"{entry['chrom']}\t{entry['start']}\t{entry['end']}\t")
            f.write(f"{new_name}\t{entry['score']}\t{entry['strand']}\t")
            f.write(f"{entry['ensembl_id']}\t{entry['refseq_id']}\t")
            f.write(f"{dna_seq}\t{protein_seq}\n")

    if verbose:
        logger.info(f"    ✓ Wrote {len(output_entries)} entries to {output_path}")

    # Summary
    if verbose:
        logger.info(f"\n{'=' * 60}")
        logger.info("Summary")
        logger.info(f"{'=' * 60}")
        logger.info(f"  Total entries processed: {stats['total_entries']}")
        logger.info(f"  Valid input entries: {stats['valid_entries']}")
        logger.info(f"  Canonical starts added: {stats['canonical_starts_added']}")
        logger.info(
            f"  Extended/Truncated entries removed (no annotated start): {stats['extended_truncated_removed']}"
        )
        logger.info(
            f"  Duplicate Annotated entries removed: {stats['duplicate_annotated_removed']}"
        )
        logger.info(f"  Final output entries: {stats['final_entries']}")
        logger.info(f"  Unique genes: {stats['final_gene_count']}")

        if stats["final_entries"] > 0:
            logger.info(
                f"  RefSeq mapping rate: {stats['refseq_mapped']}/{stats['final_entries']} ({stats['refseq_mapped'] / stats['final_entries'] * 100:.1f}%)"
            )
        else:
            logger.info(f"  RefSeq mapping rate: N/A (no output entries)")

        logger.info(f"\n  Multi-TIS Statistics:")
        if stats["multi_tis_transcripts"]:
            logger.info(f"  Transcripts with multiple TIS in same category:")
            for category in sorted(stats["multi_tis_transcripts"].keys()):
                count = stats["multi_tis_transcripts"][category]
                examples = stats.get("multi_tis_examples", {}).get(category, [])
                if examples:
                    examples_str = ", ".join(examples[:3])
                    logger.info(
                        f"    {category}: {count} transcripts (e.g., {examples_str})"
                    )
                else:
                    logger.info(f"    {category}: {count} transcripts")
            logger.info(
                f"  Total transcripts with any multi-TIS: {stats['total_transcripts_with_multi_tis']}"
            )
        else:
            logger.info(f"  No transcripts with multiple TIS in same category")

        logger.info(f"\n  Duplicate Annotated Entry Statistics:")
        if stats.get("transcripts_with_duplicate_annotated", 0) > 0:
            logger.info(
                f"  Transcripts with duplicate Annotated entries: {stats['transcripts_with_duplicate_annotated']}"
            )
            logger.info(f"  Entries removed: {stats['duplicate_annotated_removed']}")
            examples = stats.get("duplicate_annotated_examples", [])
            if examples:
                examples_str = ", ".join(examples[:3])
                logger.info(f"  Examples: {examples_str}")
            logger.info(
                f"  Detailed log: {Path(str(output_bed_path) + '.duplicates_removed.tsv')}"
            )
        else:
            logger.info(f"  No duplicate Annotated entries detected")

    return dict(stats)
