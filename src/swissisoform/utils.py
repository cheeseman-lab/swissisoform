"""Utility functions for the SwissIsoform package.

This module provides various helper functions for mutation analysis,
file handling, and data processing used throughout the package.
"""

import pandas as pd
from typing import List, Dict, Optional, Union
import logging
from pathlib import Path
from collections import defaultdict
import re
import json


# Configure logger
logger = logging.getLogger(__name__)


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
    print(f"Gene-level results saved to {output_path}")


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
        print(f"No transcript-isoform pairs found to save")
        return

    # Convert to DataFrame
    pairs_df = pd.DataFrame(all_pairs)

    # Organize columns in preferred order
    column_order = [
        "gene_name",  # Gene name first
        "transcript_id",  # Transcript ID second
        "feature_id",  # Alternative feature ID third
        "feature_type",  # truncation or extension
        "feature_start",
        "feature_end",
        "mutation_count_total",  # Total mutation count
    ]

    # Add mutation category columns to the order (they start with "mutations_")
    mutation_category_columns = [
        col for col in pairs_df.columns if col.startswith("mutations_")
    ]
    column_order.extend(mutation_category_columns)

    # Add variant IDs columns
    variant_id_columns = [
        col for col in pairs_df.columns if col.startswith("variant_ids_")
    ]
    column_order.extend(variant_id_columns)

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
    print(f"Transcript-isoform pair analysis saved to {output_path}")


def print_mutation_summary(results_df, output_dir):
    """Print summary statistics from the analysis results.

    Args:
        results_df: DataFrame containing analysis results
        output_dir: Directory where output files are saved
    """
    print("\nAnalysis Summary:")
    print(f"  ├─ Total genes processed: {len(results_df)}")
    print("\n  ├─ Status breakdown:")
    for status, count in results_df["status"].value_counts().items():
        print(f"  │  ├─ {status}: {count}")

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

        print(f"\n  ├─ Transcript-Isoform Analysis:")
        print(f"  │  ├─ Total transcripts across all genes: {total_transcripts}")
        print(f"  │  ├─ Total alternative isoform features: {total_features}")
        print(f"  │  ├─ Total transcript-isoform pairs: {total_pairs}")
        print(f"  │  └─ Average pairs per gene: {avg_pairs_per_gene:.2f}")

        # Mutation statistics if available
        if "mutations_filtered" in successful_genes.columns:
            total_mutations = successful_genes["mutations_filtered"].sum()
            print(f"\n  ├─ Mutation Analysis:")
            print(
                f"  │  ├─ Total mutations in alternative isoform regions: {total_mutations}"
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
                        pairs_df["mutation_count_total"].sum()
                        if "mutation_count_total" in pairs_df.columns
                        else 0
                    )
                    avg_mutations_per_pair = (
                        total_pair_mutations / len(pairs_df) if len(pairs_df) > 0 else 0
                    )
                    print(
                        f"  │  ├─ Total mutations across all pairs: {total_pair_mutations}"
                    )
                    print(
                        f"  │  ├─ Average mutations per transcript-isoform pair: {avg_mutations_per_pair:.2f}"
                    )

                    # Show breakdown by isoform type if available
                    if "feature_type" in pairs_df.columns:
                        print(f"  │  │")
                        print(f"  │  ├─ Breakdown by isoform type:")

                        type_mutations = pairs_df.groupby("feature_type")[
                            "mutation_count_total"
                        ].sum()
                        for isoform_type, count in type_mutations.items():
                            type_pairs = len(
                                pairs_df[pairs_df["feature_type"] == isoform_type]
                            )
                            avg_per_type = count / type_pairs if type_pairs > 0 else 0
                            print(
                                f"  │  │  ├─ {isoform_type.capitalize()}: {count} total ({avg_per_type:.1f} avg/pair)"
                            )

                    # Print statistics for each mutation category
                    mutation_categories = [
                        col for col in pairs_df.columns if col.startswith("mutations_")
                    ]
                    if mutation_categories:
                        print(f"  │  │")
                        print(f"  │  ├─ Breakdown by mutation category:")

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

                            print(
                                f"  │  │  ├─ {category_name}: {category_total} ({category_percent:.1f}%)"
                            )

                    print(
                        f"  │  └─ Detailed results available in isoform_level_results.csv"  # Updated filename
                    )
                except Exception as e:
                    logger.error(f"Error reading mutation analysis: {str(e)}")
                    print(f"  │  └─ Error reading detailed mutation analysis: {str(e)}")

    # Genes with errors
    error_genes = results_df[results_df["status"] == "error"]
    if not error_genes.empty:
        print("\n  ├─ Genes with errors:")
        for _, row in error_genes.iterrows():
            print(f"  │  ├─ {row['gene_name']}: {row['error']}")

    print(f"\n  ├─ Results saved to: {output_dir}")
    print(
        f"  ├─ Gene-level results saved to: {Path(output_dir) / 'gene_level_results.csv'}"
    )
    print(
        f"  ├─ Detailed mutation analysis by pair saved to: {Path(output_dir) / 'isoform_level_results.csv'}"  # Updated filename
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
    print(f"\nProtein Sequence Generation Summary:")
    print(f"  ├─ Total genes processed: {total_genes}")

    # Status breakdown
    print(f"\n  ├─ Status breakdown:")
    print(f"  │  ├─ Success: {successful_genes}")
    if failed_genes:
        print(f"  │  └─ Failed: {len(failed_genes)}")
    else:
        print(f"  │  └─ Failed: 0")

    # Dataset statistics
    if not dataset.empty:
        print(f"\n  ├─ Sequence Generation:")
        print(f"  │  ├─ Total sequences generated: {len(dataset):,}")

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

        print(
            f"  │  ├─ Transcript-isoform pairs: {unique_pairs}"
        )  # Updated terminology
        print(
            f"  │  └─ Average sequences per gene: {len(dataset) / genes_with_data:.1f}"
        )

        # Mode-specific breakdown
        if mutations_mode and "variant_type" in dataset.columns:
            print(f"\n  ├─ Sequence breakdown:")
            type_counts = dataset["variant_type"].value_counts()
            for variant_type, count in type_counts.items():
                percentage = (count / len(dataset)) * 100
                print(f"  │  ├─ {variant_type}: {count:,} ({percentage:.1f}%)")

            # Mutation-specific statistics
            if "canonical_mutated" in type_counts:
                mutation_data = dataset[dataset["variant_type"] == "canonical_mutated"]
                if not mutation_data.empty:
                    print(f"\n  ├─ Mutation Analysis:")
                    print(f"  │  ├─ Total mutations integrated: {len(mutation_data)}")
                    print(
                        f"  │  ├─ Unique mutation positions: {mutation_data['mutation_position'].nunique()}"
                    )

                    if "aa_change" in mutation_data.columns:
                        aa_changes = mutation_data["aa_change"].dropna()
                        silent_mutations = len(mutation_data) - len(aa_changes)
                        print(f"  │  ├─ Mutations with AA changes: {len(aa_changes)}")
                        if silent_mutations > 0:
                            print(f"  │  ├─ Silent mutations: {silent_mutations}")

                    # Breakdown by impact type if available
                    if "mutation_impact" in mutation_data.columns:
                        impact_counts = mutation_data["mutation_impact"].value_counts()
                        print(f"  │  │")
                        print(f"  │  ├─ Breakdown by impact type:")
                        for impact_type, count in impact_counts.items():
                            percentage = (count / len(mutation_data)) * 100
                            print(
                                f"  │  │  ├─ {impact_type}: {count} ({percentage:.1f}%)"
                            )

                    print(
                        f"  │  └─ Average mutations per gene: {len(mutation_data) / genes_with_data:.1f}"
                    )
        else:
            # Pairs mode breakdown
            if "is_alternative" in dataset.columns:  # Updated column name
                canonical_count = len(dataset[dataset["is_alternative"] == 0])
                alternative_count = len(dataset[dataset["is_alternative"] == 1])
                print(f"\n  ├─ Sequence breakdown:")
                print(f"  │  ├─ Canonical: {canonical_count:,}")
                print(
                    f"  │  └─ Alternative isoform: {alternative_count:,}"
                )  # Updated terminology

        # Length statistics
        print(f"\n  ├─ Length statistics:")
        print(f"  │  ├─ Average: {dataset['length'].mean():.1f} amino acids")
        print(f"  │  ├─ Range: {dataset['length'].min()}-{dataset['length'].max()}")
        print(f"  │  └─ Median: {dataset['length'].median():.1f}")

    else:
        print(f"\n  ├─ No sequences generated")

    # Failed genes details
    if failed_genes:
        print(f"\n  ├─ Genes with errors:")
        for gene in failed_genes[:5]:  # Show first 5
            print(
                f"  │  ├─ {gene}: No transcript-isoform pairs found"
            )  # Updated terminology
        if len(failed_genes) > 5:
            print(f"  │  └─ ... and {len(failed_genes) - 5} more")

    # Output files
    output_path = Path(output_dir)
    print(f"\n  ├─ Results saved to: {output_dir}")

    if mutations_mode:
        fasta_file = output_path / "protein_sequences_with_mutations.fasta"
        csv_file = output_path / "protein_sequences_with_mutations.csv"
        print(f"  ├─ Protein sequences with mutations saved to: {csv_file.name}")
        if fasta_file.exists():
            print(f"  ├─ FASTA format saved to: {fasta_file.name}")
    else:
        fasta_file = output_path / "protein_sequences_pairs.fasta"
        csv_file = output_path / "protein_sequences_pairs.csv"
        print(f"  ├─ Protein sequence pairs saved to: {csv_file.name}")
        if fasta_file.exists():
            print(f"  ├─ FASTA format saved to: {fasta_file.name}")


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
        print(
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
        print(f"Extracted {len(gencode_gene_names)} gene names from GENCODE GTF")

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
        print(f"Extracted {len(reference_gene_names)} gene names from reference GTF")

    # Create update mapping
    for ensembl_id, gencode_name in gencode_gene_names.items():
        if ensembl_id in reference_gene_names:
            reference_name = reference_gene_names[ensembl_id]

            # Only update if names differ
            if gencode_name != reference_name:
                gene_name_updates[ensembl_id] = reference_name

    if verbose:
        print(f"Created {len(gene_name_updates)} gene name updates")

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
        print(f"\nGTF Update Summary:")
        print(f"  Total lines processed: {stats['total_lines']}")
        print(f"  Genes processed: {stats['genes_processed']}")
        print(f"  Genes with updated names: {stats['genes_updated']}")
        print(f"  Total lines updated: {stats['updated_lines']}")
        print(f"  Output saved to: {output_gtf_path}")

    return stats


def parse_bed_line(line):
    """Parse a single BED file line to extract ENSEMBL ID and gene name.

    Handles format: GENE_NAME_ENSEMBL_ID_TYPE_CODON_EFFICIENCY
    Example: TM4SF1_ENSG00000169908.11_Annotated_AUG_13.62202605

    Args:
        line: BED file line to parse
    Returns:
        Dictionary containing parsed data or None if invalid
    """
    fields = line.strip().split("\t")
    if len(fields) < 4:
        return None

    name_parts = fields[3].split("_")
    if len(name_parts) < 4:  # Need at least GENE_ENSEMBL_TYPE_CODON
        return None

    # Updated parsing for your format: GENE_NAME_ENSEMBL_ID_TYPE_CODON_EFFICIENCY
    gene_name = name_parts[0]
    full_ensembl_id = name_parts[1]
    ensembl_id = full_ensembl_id.split(".")[0]  # Remove version number

    # Validate it's an Ensembl ID
    if not ensembl_id.startswith("ENSG"):
        return None

    return {
        "ensembl_id": ensembl_id,
        "full_ensembl_id": full_ensembl_id,
        "gene_name": gene_name,
        "fields": fields,
        "line": line.strip(),
    }


def parse_alternative_start_name(name: str) -> Dict[str, str]:
    """Parse alternative start site name field.

    Handles format: GENE_NAME_ENSEMBL_ID_TYPE_CODON_EFFICIENCY
    Example: TM4SF1_ENSG00000169908.11_Annotated_AUG_13.62202605

    Args:
        name: Name field from BED file

    Returns:
        Dictionary with parsed components
    """
    parts = name.split("_")

    if len(parts) >= 5:
        return {
            "gene_name": parts[0],
            "gene_id": parts[1],
            "start_type": parts[2],  # Annotated, Truncated, Extended, uORF
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


def extract_gene_mapping_from_gtf(
    gtf_path: Union[str, Path], is_ucsc: bool = False
) -> Dict[str, str]:
    """Extract gene ID to gene name mappings from a GTF file.

    Args:
        gtf_path: Path to the GTF annotation file
        is_ucsc: Whether the GTF is from UCSC (different format than GENCODE)

    Returns:
        Dictionary mapping gene IDs to gene names
    """
    gtf_path = Path(gtf_path)
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")

    reference_data = {}
    gene_count = 0

    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#"):  # Skip comment lines
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:  # Ensure we have attribute field
                continue

            # Extract attributes to dictionary
            attributes = fields[8]
            attr_dict = {}
            for attr in attributes.split(";"):
                attr = attr.strip()
                if not attr:
                    continue

                parts = attr.split(" ", 1)
                if len(parts) == 2:
                    key, value = parts
                    # Remove quotes from value
                    value = value.strip('"')
                    attr_dict[key] = value

            if is_ucsc:
                # UCSC format - gene_id and gene_name directly available
                if "gene_id" in attr_dict and "gene_name" in attr_dict:
                    gene_id = attr_dict["gene_id"]
                    gene_name = attr_dict["gene_name"]

                    # Store mapping with gene_id as key
                    reference_data[gene_id] = gene_name
                    gene_count += 1
            else:
                # GENCODE format
                if (
                    fields[2] == "gene"
                    and "gene_id" in attr_dict
                    and "gene_name" in attr_dict
                ):
                    full_gene_id = attr_dict["gene_id"]
                    gene_name = attr_dict["gene_name"]

                    # Store both with and without version number
                    reference_data[full_gene_id] = gene_name

                    # Also store base gene ID without version
                    if "." in full_gene_id:
                        base_gene_id = full_gene_id.split(".")[0]
                        reference_data[base_gene_id] = gene_name

                    gene_count += 1

    print(f"Extracted {gene_count} unique gene ID to name mappings from GTF")
    return reference_data


def comprehensive_cleanup_bed_with_transcripts(
    input_bed_path: Union[str, Path],
    gtf_path: Union[str, Path],
    preferred_transcripts_path: Union[str, Path],
    output_bed_path: Union[str, Path],
    verbose: bool = True,
) -> Dict:
    """Comprehensive BED cleanup with transcript mapping for each entry.

    This function:
    1. Parses and validates all BED entries
    2. Maps each entry to its best matching transcript
    3. Adds missing annotated starts from GTF
    4. Uses biologically-relevant annotated start selection
    5. Outputs enhanced BED with transcript_id column
    6. Filters to complete genes only

    Args:
        input_bed_path: Path to input BED file
        gtf_path: Path to GTF annotation file
        preferred_transcripts_path: Path to preferred transcripts file
        output_bed_path: Path to save enhanced BED file
        verbose: Print detailed progress

    Returns:
        Dictionary with processing statistics
    """
    import re
    from collections import defaultdict

    input_bed_path = Path(input_bed_path)
    gtf_path = Path(gtf_path)
    output_bed_path = Path(output_bed_path)

    if verbose:
        print("Starting comprehensive BED cleanup with transcript mapping...")
        print(f"  ├─ Input BED: {input_bed_path}")
        print(f"  ├─ GTF: {gtf_path}")
        print(f"  └─ Output: {output_bed_path}")

    # Load gene name mappings and preferred transcripts
    if verbose:
        print("\nLoading reference data...")

    gene_name_mapping = extract_gene_mapping_from_gtf(gtf_path)
    preferred_transcripts = load_preferred_transcripts(preferred_transcripts_path)

    if verbose:
        print(f"  ├─ Gene mappings: {len(gene_name_mapping)}")
        print(f"  └─ Preferred transcripts: {len(preferred_transcripts)}")

    # Build transcript position index from GTF
    if verbose:
        print("\nBuilding transcript index...")

    transcript_info = {}  # transcript_id -> (gene_id, chrom, start_pos, end_pos, strand)
    gene_transcripts = defaultdict(list)  # gene_id -> [transcript_info]

    # Build transcript position index from GTF start_codon features
    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#") or "\tstart_codon\t" not in line:
                continue

            fields = line.strip().split("\t")
            chrom, start, end, strand = (
                fields[0],
                int(fields[3]),
                int(fields[4]),
                fields[6],
            )

            gene_match = re.search(r'gene_id "([^"]+)"', fields[8])
            transcript_match = re.search(r'transcript_id "([^"]+)"', fields[8])

            if gene_match and transcript_match:
                gene_id = gene_match.group(1).split(".")[0]
                transcript_id = transcript_match.group(1)
                start_pos = start if strand == "+" else end

                transcript_info[transcript_id] = (
                    gene_id,
                    chrom,
                    start_pos,
                    end,
                    strand,
                )
                gene_transcripts[gene_id].append(
                    {
                        "transcript_id": transcript_id,
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "start_pos": start_pos,
                        "strand": strand,
                    }
                )

    if verbose:
        print(
            f"  └─ Indexed {len(transcript_info)} transcripts for {len(gene_transcripts)} genes"
        )

    # Parse and validate BED entries
    if verbose:
        print("\nParsing BED entries...")

    valid_entries = []
    gene_entries = defaultdict(list)  # gene_id -> [entry_info]

    stats = {
        "total_entries": 0,
        "invalid_format": 0,
        "invalid_gene": 0,
        "uorfs_filtered": 0,
        "valid_entries": 0,
        "updated_gene_names": 0,
    }

    with open(input_bed_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            if not line.strip():
                continue

            stats["total_entries"] += 1
            fields = line.strip().split("\t")

            if len(fields) < 6:
                stats["invalid_format"] += 1
                continue

            # Parse name field
            try:
                name_info = parse_alternative_start_name(fields[3])
                gene_id = name_info["gene_id"].split(".")[0]

                # Skip problematic entries
                if '"' in name_info["gene_id"] or gene_id not in gene_name_mapping:
                    stats["invalid_gene"] += 1
                    continue

                # Filter uORFs
                if name_info["start_type"] == "uORF":
                    stats["uorfs_filtered"] += 1
                    continue

                # Only keep target types
                if name_info["start_type"] not in [
                    "Annotated",
                    "Truncated",
                    "Extended",
                ]:
                    continue

                # Calculate start position based on strand
                chrom, start, end, strand = (
                    fields[0],
                    int(fields[1]),
                    int(fields[2]),
                    fields[5],
                )
                start_pos = start if strand == "+" else end

                # Update gene name if needed
                ref_gene_name = gene_name_mapping[gene_id]
                updated_name = fields[3]
                if name_info["gene_name"].upper() != ref_gene_name.upper():
                    # Reconstruct the name with updated gene name
                    updated_name = f"{ref_gene_name}_{name_info['gene_id']}_{name_info['start_type']}_{name_info['start_codon']}_{name_info['efficiency']}"
                    stats["updated_gene_names"] = stats.get("updated_gene_names", 0) + 1

                entry_info = {
                    "line_num": line_num,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "start_pos": start_pos,
                    "gene_id": gene_id,
                    "gene_name": ref_gene_name,  # Use the reference gene name
                    "start_type": name_info["start_type"],
                    "start_codon": name_info["start_codon"],
                    "efficiency": name_info["efficiency"],
                    "original_name": updated_name,  # Use the updated name
                    "score": fields[4],
                }

                valid_entries.append(entry_info)
                gene_entries[gene_id].append(entry_info)
                stats["valid_entries"] += 1

            except (IndexError, ValueError, KeyError):
                stats["invalid_format"] += 1
                continue

    if verbose:
        print(f"  ├─ Total entries: {stats['total_entries']}")
        print(f"  ├─ Valid entries: {stats['valid_entries']}")
        print(f"  ├─ Invalid format: {stats['invalid_format']}")
        print(f"  ├─ Invalid genes: {stats['invalid_gene']}")
        print(f"  ├─ Gene names updated: {stats['updated_gene_names']}")
        print(f"  └─ uORFs filtered: {stats['uorfs_filtered']}")

    # Add missing annotated starts
    if verbose:
        print("\nAdding missing annotated starts...")

    genes_needing_annotated = set()
    for gene_id, entries in gene_entries.items():
        has_annotated = any(e["start_type"] == "Annotated" for e in entries)
        has_alternatives = any(
            e["start_type"] in ["Truncated", "Extended"] for e in entries
        )

        if has_alternatives and not has_annotated:
            genes_needing_annotated.add(gene_id)

    added_annotated = 0
    for gene_id in genes_needing_annotated:
        if gene_id not in gene_transcripts:
            continue

        # Select best transcript (prefer from preferred list)
        gene_transcript_list = gene_transcripts[gene_id]
        selected_transcript = None

        for transcript_info_dict in gene_transcript_list:
            transcript_id = transcript_info_dict["transcript_id"]
            if transcript_id in preferred_transcripts:
                selected_transcript = transcript_info_dict
                break

        if not selected_transcript:
            selected_transcript = gene_transcript_list[0]  # Fallback to first

        # Create annotated entry
        chrom = selected_transcript["chrom"]
        strand = selected_transcript["strand"]

        if strand == "+":
            bed_start = selected_transcript["start"]
            bed_end = selected_transcript["start"] + 2
        else:
            bed_start = selected_transcript["end"] - 2
            bed_end = selected_transcript["end"]

        gene_name = gene_name_mapping[gene_id]

        new_entry = {
            "line_num": -1,  # Mark as added
            "chrom": chrom,
            "start": bed_start,
            "end": bed_end,
            "strand": strand,
            "start_pos": selected_transcript["start_pos"],
            "gene_id": gene_id,
            "gene_name": gene_name,
            "start_type": "Annotated",
            "start_codon": "AUG",
            "efficiency": 0.0,
            "original_name": f"{gene_name}_{gene_id}_Annotated_AUG_0.0",
            "score": "0",
        }

        gene_entries[gene_id].append(new_entry)
        valid_entries.append(new_entry)
        added_annotated += 1

    if verbose:
        print(f"  └─ Added {added_annotated} annotated starts")

    # Map entries to transcripts with biological relevance
    if verbose:
        print("\nMapping entries to transcripts...")

    def select_relevant_annotated_start(
        gene_id, alternative_pos, alternative_type, strand
    ):
        """Select biologically relevant annotated start."""
        annotated_entries = [
            e for e in gene_entries[gene_id] if e["start_type"] == "Annotated"
        ]

        if len(annotated_entries) <= 1:
            return annotated_entries[0] if annotated_entries else None

        relevant_starts = []

        for ann_entry in annotated_entries:
            ann_pos = ann_entry["start_pos"]

            if alternative_type == "Extended":
                if strand == "+":
                    # Extension upstream, canonical downstream
                    if ann_pos > alternative_pos:
                        relevant_starts.append((ann_entry, ann_pos - alternative_pos))
                else:
                    # Extension upstream in transcript (lower genomic coord)
                    if ann_pos < alternative_pos:
                        relevant_starts.append((ann_entry, alternative_pos - ann_pos))

            elif alternative_type == "Truncated":
                if strand == "+":
                    # Truncation downstream, canonical upstream
                    if ann_pos < alternative_pos:
                        relevant_starts.append((ann_entry, alternative_pos - ann_pos))
                else:
                    # Truncation downstream in transcript (higher genomic coord)
                    if ann_pos > alternative_pos:
                        relevant_starts.append((ann_entry, ann_pos - alternative_pos))

        if relevant_starts:
            return min(relevant_starts, key=lambda x: x[1])[0]
        else:
            # Fallback to highest efficiency
            return max(annotated_entries, key=lambda x: x["efficiency"])

    def find_best_transcript(entry, reference_annotated=None):
        """Find best matching transcript for an entry."""
        gene_id = entry["gene_id"]

        if gene_id not in gene_transcripts:
            return "NA"

        # For annotated starts, find closest transcript
        if entry["start_type"] == "Annotated":
            best_transcript = None
            best_distance = float("inf")

            for transcript_dict in gene_transcripts[gene_id]:
                distance = abs(transcript_dict["start_pos"] - entry["start_pos"])
                if distance < best_distance:
                    best_distance = distance
                    best_transcript = transcript_dict["transcript_id"]

            return best_transcript or "NA"

        # For alternatives, use the transcript from relevant annotated start
        else:
            if reference_annotated:
                ref_transcript = find_best_transcript(reference_annotated)
                return ref_transcript
            else:
                return "NA"

    enhanced_entries = []
    transcript_assignments = 0

    for gene_id, entries in gene_entries.items():
        # Filter to complete genes (have both annotated and alternatives)
        has_annotated = any(e["start_type"] == "Annotated" for e in entries)
        has_alternatives = any(
            e["start_type"] in ["Truncated", "Extended"] for e in entries
        )

        if not (has_annotated and has_alternatives):
            continue

        for entry in entries:
            if entry["start_type"] == "Annotated":
                transcript_id = find_best_transcript(entry)
            else:
                # Find relevant annotated start
                relevant_annotated = select_relevant_annotated_start(
                    gene_id, entry["start_pos"], entry["start_type"], entry["strand"]
                )
                transcript_id = find_best_transcript(entry, relevant_annotated)

            # Create enhanced BED line
            enhanced_line = (
                f"{entry['chrom']}\t{entry['start']}\t{entry['end']}\t"
                f"{entry['original_name']}\t{entry['score']}\t{entry['strand']}\t{transcript_id}"
            )

            enhanced_entries.append(enhanced_line)
            if transcript_id != "NA":
                transcript_assignments += 1

    if verbose:
        print(f"  ├─ Enhanced entries: {len(enhanced_entries)}")
        print(f"  └─ Transcript assignments: {transcript_assignments}")

    # Write output
    if verbose:
        print(f"\nWriting enhanced BED file...")

    output_bed_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_bed_path, "w") as f:
        f.write("\n".join(enhanced_entries))

    # Final statistics
    final_genes = len(
        set(
            re.search(r"_([^_]+)_", line).group(1).split(".")[0]
            for line in enhanced_entries
            if re.search(r"_([^_]+)_", line)
        )
    )

    final_stats = {
        **stats,
        "added_annotated": added_annotated,
        "enhanced_entries": len(enhanced_entries),
        "transcript_assignments": transcript_assignments,
        "final_genes": final_genes,
    }

    if verbose:
        print(f"\nComprehensive cleanup complete:")
        print(f"  ├─ Enhanced entries: {len(enhanced_entries)}")
        print(f"  ├─ Final genes: {final_genes}")
        print(f"  ├─ Added annotated starts: {added_annotated}")
        print(f"  ├─ Transcript assignments: {transcript_assignments}")
        print(f"  └─ Output: {output_bed_path}")

    return final_stats


def extract_ensembl_refseq_mapping_from_refseq_gtf_direct(
    refseq_gtf_path: Union[str, Path], verbose: bool = True
) -> tuple[Dict[str, str], Dict[str, List[str]], Dict[str, List[str]]]:
    """Extract direct Ensembl to RefSeq mapping from RefSeq GTF cross-references.

    Args:
        refseq_gtf_path: Path to RefSeq GTF file
        verbose: Print extraction progress

    Returns:
        Tuple of (ensembl_to_refseq_dict, refseq_to_ensembl_dict, gene_to_refseq_with_ensembl_dict)
    """
    import re
    from collections import defaultdict

    ensembl_to_refseq = {}
    refseq_to_ensembl = defaultdict(list)
    gene_to_refseq_with_ensembl = defaultdict(list)
    mapped_transcripts = 0

    with open(refseq_gtf_path, "r") as f:
        for line in f:
            # Skip comments and non-transcript lines
            if line.startswith("#") or "\ttranscript\t" not in line:
                continue

            attributes = line.strip().split("\t")[8]

            # Extract gene name, RefSeq transcript ID, and Ensembl cross-reference
            gene_match = re.search(r'gene "([^"]+)"', attributes) or re.search(
                r'gene_id "([^"]+)"', attributes
            )
            refseq_match = re.search(r'transcript_id "([^"]+)"', attributes)
            ensembl_match = re.search(r'db_xref "Ensembl:([^"]+)"', attributes)

            if refseq_match:
                refseq_id = refseq_match.group(1)
                gene_name = gene_match.group(1) if gene_match else "Unknown"

                # Only keep RefSeq transcript IDs (NM_, NR_, XM_, XR_)
                if refseq_id.startswith(("NM_", "NR_", "XM_", "XR_")):
                    # If this RefSeq transcript has an Ensembl cross-reference
                    if ensembl_match:
                        ensembl_id = ensembl_match.group(1)

                        # Store direct mappings
                        ensembl_to_refseq[ensembl_id] = refseq_id
                        base_ensembl = ensembl_id.split(".")[0]
                        ensembl_to_refseq[base_ensembl] = refseq_id

                        # Store reverse mapping
                        refseq_to_ensembl[refseq_id].append(ensembl_id)

                        # Store gene-level mapping
                        gene_to_refseq_with_ensembl[gene_name].append(refseq_id)

                        mapped_transcripts += 1

    if verbose:
        print(
            f"  ├─ Found {mapped_transcripts:,} RefSeq transcripts with Ensembl cross-references"
        )
        print(f"  └─ Covering {len(gene_to_refseq_with_ensembl):,} genes")

    return ensembl_to_refseq, refseq_to_ensembl, gene_to_refseq_with_ensembl


def find_refseq_gtf_file(genome_data_dir: Union[str, Path]) -> Optional[Path]:
    """Find available RefSeq GTF file in the genome data directory.

    Args:
        genome_data_dir: Directory containing genome files

    Returns:
        Path to RefSeq GTF file or None if not found
    """
    genome_data_dir = Path(genome_data_dir)

    # Possible RefSeq GTF filenames (in order of preference)
    possible_files = [
        "GRCh38_latest_genomic.gtf",  # NCBI RefSeq
        "ncbiRefSeq.gtf",  # UCSC RefSeq
        "refseq_genomic.gtf",  # Alternative naming
        "GRCh38_refseq.gtf",  # Alternative naming
    ]

    for filename in possible_files:
        gtf_path = genome_data_dir / filename
        if gtf_path.exists():
            return gtf_path

    return None


def comprehensive_cleanup_bed_with_refseq_gtf_direct(
    input_bed_path: Union[str, Path],
    gtf_path: Union[str, Path],
    preferred_transcripts_path: Union[str, Path],
    output_bed_path: Union[str, Path],
    refseq_gtf_path: Optional[Union[str, Path]] = None,
    verbose: bool = True,
) -> Dict:
    """Comprehensive BED cleanup using direct RefSeq GTF cross-references with fallback logic.

    Args:
        input_bed_path: Path to input BED file
        gtf_path: Path to GENCODE GTF annotation file
        preferred_transcripts_path: Path to preferred transcripts file
        output_bed_path: Path to save enhanced BED file
        refseq_gtf_path: Path to RefSeq GTF file (auto-detected if None)
        verbose: Print detailed progress

    Returns:
        Dictionary with processing statistics
    """
    from pathlib import Path

    # Step 1: Find RefSeq GTF if not provided
    if refseq_gtf_path is None:
        genome_data_dir = Path(gtf_path).parent
        refseq_gtf_path = find_refseq_gtf_file(genome_data_dir)

        if refseq_gtf_path is None:
            if verbose:
                print("  ❌ No RefSeq GTF found!")
                print("  💡 Please run: bash 0_download_genome.sh")
            raise FileNotFoundError(
                f"RefSeq GTF file not found in {genome_data_dir}. "
                f"Please run 0_download_genome.sh to download RefSeq annotations."
            )

    if verbose:
        print(f"  ├─ Using RefSeq GTF: {refseq_gtf_path.name}")

    # Step 2: Run standard cleanup to get Ensembl transcript mappings
    temp_bed_path = str(output_bed_path).replace(".bed", "_temp.bed")

    standard_stats = comprehensive_cleanup_bed_with_transcripts(
        input_bed_path=input_bed_path,
        gtf_path=gtf_path,
        preferred_transcripts_path=preferred_transcripts_path,
        output_bed_path=temp_bed_path,
        verbose=verbose,
    )

    # Step 3: Extract RefSeq mapping from GTF cross-references
    if verbose:
        print("  ├─ Extracting RefSeq mappings...")

    ensembl_to_refseq, refseq_to_ensembl, gene_to_refseq_with_ensembl = (
        extract_ensembl_refseq_mapping_from_refseq_gtf_direct(
            refseq_gtf_path=refseq_gtf_path, verbose=verbose
        )
    )

    # Step 4: Add RefSeq column to BED file with enhanced fallback logic
    if verbose:
        print("  ├─ Mapping transcripts to RefSeq...")

    enhanced_entries = []
    refseq_assignments = 0
    direct_matches = 0
    fallback_matches = 0
    gene_level_matches = 0

    with open(temp_bed_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) >= 7:
                ensembl_transcript_id = parts[6]
                # Extract gene name from BED name field (assumes format: GENE_ENSG_...)
                bed_name_parts = parts[3].split("_")
                gene_name = bed_name_parts[0] if bed_name_parts else "Unknown"

                refseq_id = "NA"

                # Method 1: Direct lookup (versioned)
                if ensembl_transcript_id in ensembl_to_refseq:
                    refseq_id = ensembl_to_refseq[ensembl_transcript_id]
                    direct_matches += 1
                else:
                    # Method 2: Direct lookup (base, without version)
                    base_ensembl = ensembl_transcript_id.split(".")[0]
                    if base_ensembl in ensembl_to_refseq:
                        refseq_id = ensembl_to_refseq[base_ensembl]
                        direct_matches += 1
                    else:
                        # Method 3: Fallback - find RefSeq with only one Ensembl match
                        for (
                            refseq_transcript,
                            ensembl_list,
                        ) in refseq_to_ensembl.items():
                            if len(ensembl_list) == 1:
                                target_ensembl = ensembl_list[0]
                                if (
                                    ensembl_transcript_id == target_ensembl
                                    or base_ensembl == target_ensembl.split(".")[0]
                                ):
                                    refseq_id = refseq_transcript
                                    fallback_matches += 1
                                    break

                        # Method 4: Gene-level fallback
                        if (
                            refseq_id == "NA"
                            and gene_name in gene_to_refseq_with_ensembl
                        ):
                            refseq_candidates = gene_to_refseq_with_ensembl[gene_name]
                            unique_refseq = list(set(refseq_candidates))

                            if len(unique_refseq) == 1:
                                refseq_id = unique_refseq[0]
                                gene_level_matches += 1
                            elif len(unique_refseq) > 1:
                                # Pick the best one (prefer NM_ over XM_, then shortest)
                                best_refseq = min(
                                    unique_refseq,
                                    key=lambda x: (not x.startswith("NM_"), len(x), x),
                                )
                                refseq_id = best_refseq
                                gene_level_matches += 1

                if refseq_id != "NA":
                    refseq_assignments += 1

                # Create enhanced line
                enhanced_line = line.strip() + f"\t{refseq_id}"
                enhanced_entries.append(enhanced_line)

    # Write final enhanced BED file
    output_bed_path = Path(output_bed_path)
    output_bed_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_bed_path, "w") as f:
        f.write("\n".join(enhanced_entries))

    # Clean up temporary file
    Path(temp_bed_path).unlink()

    # Calculate final statistics
    refseq_mapping_rate = (
        (refseq_assignments / len(enhanced_entries) * 100) if enhanced_entries else 0
    )

    final_stats = {
        **standard_stats,
        "refseq_assignments": refseq_assignments,
        "refseq_mapping_rate": refseq_mapping_rate,
        "direct_matches": direct_matches,
        "fallback_matches": fallback_matches,
        "gene_level_matches": gene_level_matches,
        "mapping_source": "refseq_gtf_cross_references_with_gene_fallback",
    }

    if verbose:
        total_fallback = fallback_matches + gene_level_matches
        print(f"  ├─ Mapped: {direct_matches} direct, {total_fallback} fallback")
        print(f"  └─ RefSeq success: {refseq_mapping_rate:.1f}%")

    return final_stats


def parse_gene_list(gene_list_path: Union[str, Path]) -> List[str]:
    """Parse a file containing a list of gene names.

    Args:
        gene_list_path: Path to file containing gene names
    Returns:
        List of gene names
    """
    gene_list_path = Path(gene_list_path)
    if not gene_list_path.exists():
        raise FileNotFoundError(f"Gene list file not found: {gene_list_path}")
    with open(gene_list_path, "r") as f:
        gene_names = [line.strip() for line in f if line.strip()]
    print(f"Read {len(gene_names)} genes from {gene_list_path}")
    return gene_names


def load_preferred_transcripts(file_path: str) -> set:
    """Load a list of preferred transcript IDs from a file.

    Args:
        file_path: Path to the file containing transcript IDs, one per line

    Returns:
        Set of preferred transcript IDs
    """
    with open(file_path, "r") as f:
        # Strip whitespace and ignore empty lines
        transcripts = {line.strip() for line in f if line.strip()}
    return transcripts


def subset_gene_list(
    gene_list_path: Union[str, Path],
    subset_genes: List[str] = [
        # Truncations of interest
        "AKR7A2",
        "ALDH9A1",
        "C15orf40",
        "CHCHD1",
        "CMPK1",
        "ERGIC3",
        "FAAH2",
        "FH",
        "GADD45GIP1",
        "GARS1",
        "GLRX2",
        "GSR",
        "LAGE3",
        "MYG1",
        "NAXE",
        "NTHL1",
        "PCBD2",
        "PNPO",
        "REXO2",
        "TRNT1",
        "UXS1",
        # Extensions of interest
        "PTEN",
        "FTL",
        "MAPK14",
    ],
) -> List[str]:
    """Subset a gene list to only include specified genes.

    Args:
        gene_list_path: Path to the file containing the full gene list
        subset_genes: List of gene names to include in the subset (defaults to a predefined list)

    Returns:
        List of genes from the full list that are in the subset
    """
    full_gene_list = parse_gene_list(gene_list_path)
    subset_set = set(subset_genes)
    subsetted_genes = [gene for gene in full_gene_list if gene in subset_set]
    print(f"Subsetted to {len(subsetted_genes)} genes from {len(full_gene_list)} total")
    return subsetted_genes
