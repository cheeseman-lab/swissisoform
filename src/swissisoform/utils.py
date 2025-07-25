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


def save_truncation_level_results(
    results: List[Dict],
    output_dir: Union[str, Path],
    filename: str = "truncation_level_results.csv",
) -> None:
    """Save detailed transcript-truncation pair analysis results to CSV file.

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
        print(f"No transcript-truncation pairs found to save")
        return
    # Convert to DataFrame
    pairs_df = pd.DataFrame(all_pairs)
    # Organize columns in preferred order
    column_order = [
        "gene_name",  # Gene name first
        "transcript_id",  # Transcript ID second
        "truncation_id",  # Alternative feature ID third
        "truncation_start",
        "truncation_end",
        "mutation_count_total",  # Total mutation count
    ]

    # Add mutation category columns to the order (they start with "mutations_")
    mutation_category_columns = [
        col for col in pairs_df.columns if col.startswith("mutations_")
    ]
    column_order.extend(mutation_category_columns)

    # Add the ClinVar IDs column
    clinvar_columns = [col for col in pairs_df.columns if col.startswith("clinvar")]
    column_order.extend(clinvar_columns)

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
    print(f"Transcript-truncation pair analysis saved to {output_path}")


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

    # Transcript-truncation statistics
    successful_genes = results_df[results_df["status"] == "success"]
    if not successful_genes.empty:
        total_transcripts = successful_genes["total_transcripts"].sum()
        total_truncations = successful_genes["truncation_features"].sum()
        total_pairs = successful_genes["transcript_truncation_pairs"].sum()

        # Calculate average pairs per gene
        avg_pairs_per_gene = (
            total_pairs / len(successful_genes) if len(successful_genes) > 0 else 0
        )

        print(f"\n  ├─ Transcript-Truncation Analysis:")
        print(f"  │  ├─ Total transcripts across all genes: {total_transcripts}")
        print(f"  │  ├─ Total truncation features: {total_truncations}")
        print(f"  │  ├─ Total transcript-truncation pairs: {total_pairs}")
        print(f"  │  └─ Average pairs per gene: {avg_pairs_per_gene:.2f}")

        # Mutation statistics if available
        if "mutations_filtered" in successful_genes.columns:
            total_mutations = successful_genes["mutations_filtered"].sum()
            print(f"\n  ├─ Mutation Analysis:")
            print(f"  │  ├─ Total mutations in truncation regions: {total_mutations}")

            # Try to load mutation analysis results to report statistics
            mutation_analysis_path = Path(output_dir) / "mutation_analysis_by_pair.csv"
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
                        f"  │  ├─ Average mutations per transcript-truncation pair: {avg_mutations_per_pair:.2f}"
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
                        f"  │  └─ Detailed results available in mutation_analysis_by_pair.csv"
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
        f"  ├─ Detailed mutation analysis by pair saved to: {Path(output_dir) / 'mutation_analysis_by_pair.csv'}"
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

        # Calculate transcript-truncation pairs
        genes_with_data = dataset["gene"].nunique()
        if mutations_mode and "variant_type" in dataset.columns:
            # Count unique transcript-truncation pairs (canonical + truncated base pairs)
            base_sequences = dataset[
                dataset["variant_type"].isin(["canonical", "truncated"])
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

        print(f"  │  ├─ Transcript-truncation pairs: {unique_pairs}")
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
            if "is_truncated" in dataset.columns:
                canonical_count = len(dataset[dataset["is_truncated"] == 0])
                truncated_count = len(dataset[dataset["is_truncated"] == 1])
                print(f"\n  ├─ Sequence breakdown:")
                print(f"  │  ├─ Canonical: {canonical_count:,}")
                print(f"  │  └─ Truncated: {truncated_count:,}")

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
            print(f"  │  ├─ {gene}: No transcript-truncation pairs found")
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

    Args:
        line: BED file line to parse
    Returns:
        Dictionary containing parsed data or None if invalid
    """
    fields = line.strip().split("\t")
    if len(fields) < 4:
        return None
    name_parts = fields[3].split("_")
    if len(name_parts) < 3:
        return None
    full_ensembl_id = name_parts[0]
    ensembl_id = full_ensembl_id.split(".")[0]
    gene_name = name_parts[1]
    return {
        "ensembl_id": ensembl_id,
        "full_ensembl_id": full_ensembl_id,
        "gene_name": gene_name,
        "fields": fields,
        "line": line.strip(),
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


def cleanup_bed(
    input_bed_path: Union[str, Path],
    output_bed_path: Union[str, Path],
    gtf_path: Union[str, Path] = "../data/genome_data/gencode.v25.annotation.gtf",
    verbose: bool = True,
) -> Dict:
    """Clean BED file: fix gene names, remove invalid entries and duplicates.

    Args:
        input_bed_path: Path to input BED file
        output_bed_path: Path to save cleaned BED file
        gtf_path: Path to GTF annotation file for gene name mappings
        verbose: Print detailed cleanup summary

    Returns:
        Dictionary containing cleanup statistics
    """
    input_bed_path, output_bed_path = Path(input_bed_path), Path(output_bed_path)

    # Get reference data from GTF instead of Ensembl Biomart
    if verbose:
        print(f"Extracting gene mapping from GTF: {gtf_path}")
    reference_data = extract_gene_mapping_from_gtf(gtf_path)
    if verbose:
        print(f"Retrieved {len(reference_data)} gene name mappings")

    # Process file (rest of function remains the same)
    valid_entries, seen_positions = [], set()
    stats = {
        "total": 0,
        "invalid_format": 0,
        "invalid_ensembl": 0,
        "duplicates": 0,
        "updated": 0,
        "valid": 0,
    }

    with open(input_bed_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            if not line.strip():
                continue
            stats["total"] += 1

            # Validate format and Ensembl ID
            parsed = parse_bed_line(line)
            if not parsed or not parsed["ensembl_id"].startswith("ENSG"):
                stats["invalid_format"] += 1
                continue
            if parsed["ensembl_id"] not in reference_data:
                stats["invalid_ensembl"] += 1
                continue

            # Update gene name if needed
            ref_gene_name = reference_data[parsed["ensembl_id"]]
            fields = parsed["fields"].copy()
            name_parts = fields[3].split("_")

            if name_parts[1].upper() != ref_gene_name.upper():
                name_parts[1] = ref_gene_name
                fields[3] = "_".join(name_parts)
                stats["updated"] += 1

            # Check for duplicates
            position_key = f"{parsed['ensembl_id']}_{fields[0]}_{fields[1]}_{fields[2]}"
            if position_key in seen_positions:
                stats["duplicates"] += 1
                continue

            # Add valid entry
            seen_positions.add(position_key)
            valid_entries.append("\t".join(fields))
            stats["valid"] += 1

    # Write output
    output_bed_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_bed_path, "w") as f_out:
        f_out.write("\n".join(valid_entries))

    if verbose:
        print(f"\nCleanup Summary:")
        print(f"  Total entries: {stats['total']}")
        print(
            f"  Invalid entries removed: {stats['invalid_format'] + stats['invalid_ensembl']}"
        )
        print(f"  Duplicates removed: {stats['duplicates']}")
        print(f"  Gene names updated: {stats['updated']}")
        print(f"  Valid entries in final file: {stats['valid']}")

    return stats


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
        "PPA2",
        "REXO2",
        "TRNT1",
        "UXS1",
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
