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


def cleanup_bed(
    input_bed_path: Union[str, Path],
    output_bed_path: Union[str, Path],
    gtf_path: Union[str, Path],
    preferred_transcripts_path: Optional[Union[str, Path]] = None,
    filter_uorfs: bool = True,
    add_annotated_starts: bool = True,
    verbose: bool = True,
) -> Dict:
    """Enhanced BED file cleanup with corrected parsing for your data format.

    Args:
        input_bed_path: Path to input BED file
        output_bed_path: Path to save cleaned BED file
        gtf_path: Path to GTF annotation file
        preferred_transcripts_path: Path to file with preferred transcript IDs
        filter_uorfs: Whether to remove uORF entries
        add_annotated_starts: Whether to add annotated start sites from GTF
        verbose: Print detailed cleanup summary

    Returns:
        Dictionary containing cleanup statistics
    """
    input_bed_path, output_bed_path = Path(input_bed_path), Path(output_bed_path)

    # Load preferred transcripts if provided
    preferred_transcripts = set()
    if preferred_transcripts_path and Path(preferred_transcripts_path).exists():
        preferred_transcripts = load_preferred_transcripts(preferred_transcripts_path)
        if verbose:
            print(f"Loaded {len(preferred_transcripts)} preferred transcripts")

    # Get gene name mappings from GTF
    if verbose:
        print(f"Extracting gene mapping from GTF: {gtf_path}")
    gene_name_mapping = extract_gene_mapping_from_gtf(gtf_path)

    # Extract transcript start sites from GTF (same as before)
    transcript_start_sites = {}
    gene_transcript_mapping = {}

    if add_annotated_starts:
        if verbose:
            print("Extracting transcript start sites from GTF...")

        with open(gtf_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 9 or fields[2] != "transcript":
                    continue

                # Parse attributes
                attributes = fields[8]
                attr_dict = {}
                for attr in attributes.split(";"):
                    attr = attr.strip()
                    if not attr:
                        continue
                    parts = attr.split(" ", 1)
                    if len(parts) == 2:
                        key, value = parts
                        attr_dict[key] = value.strip('"')

                if "gene_id" in attr_dict and "transcript_id" in attr_dict:
                    gene_id = attr_dict["gene_id"].split(".")[0]  # Remove version
                    transcript_id = attr_dict["transcript_id"]
                    transcript_start = int(fields[3])
                    transcript_end = int(fields[4])
                    strand = fields[6]
                    chromosome = fields[0]

                    transcript_info = {
                        "transcript_id": transcript_id,
                        "start": transcript_start,
                        "end": transcript_end,
                        "strand": strand,
                        "chromosome": chromosome,
                        "gene_id": gene_id,
                    }

                    if gene_id not in gene_transcript_mapping:
                        gene_transcript_mapping[gene_id] = []
                    gene_transcript_mapping[gene_id].append(transcript_info)

        if verbose:
            print(f"Extracted transcript info for {len(gene_transcript_mapping)} genes")

    # Process BED file with corrected parsing
    valid_entries = []
    seen_positions = set()
    added_annotated_starts = []

    stats = {
        "total": 0,
        "invalid_format": 0,
        "invalid_ensembl": 0,
        "duplicates": 0,
        "uorfs_filtered": 0,
        "updated_gene_names": 0,
        "valid_alternatives": 0,
        "annotated_starts_added": 0,
        "genes_with_alternatives": set(),
        "genes_missing_annotated": set(),
        "parsing_errors": [],
    }

    # First pass: collect existing entries
    existing_entries = []
    genes_with_existing_annotated = set()

    with open(input_bed_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            if not line.strip():
                continue
            stats["total"] += 1

            parsed = parse_bed_line(line)
            if not parsed:
                stats["invalid_format"] += 1
                if verbose and len(stats["parsing_errors"]) < 5:  # Show first 5 errors
                    stats["parsing_errors"].append(
                        f"Line {line_num}: {line.strip()[:100]}..."
                    )
                continue

            ensembl_id = parsed["ensembl_id"]
            if not ensembl_id.startswith("ENSG"):
                stats["invalid_ensembl"] += 1
                continue

            if ensembl_id not in gene_name_mapping:
                stats["invalid_ensembl"] += 1
                continue

            # Parse the name field using the corrected parser
            name_info = parse_alternative_start_name(parsed["fields"][3])
            start_type = name_info["start_type"]

            # Filter uORFs if requested
            if filter_uorfs and start_type == "uORF":
                stats["uorfs_filtered"] += 1
                continue

            # Track genes with annotated starts
            if start_type == "Annotated":
                genes_with_existing_annotated.add(ensembl_id)

            # Track genes with alternatives
            if start_type in ["Truncated", "Extended"]:
                stats["genes_with_alternatives"].add(ensembl_id)

            existing_entries.append((parsed, name_info, line.strip()))

    # Second pass: process existing entries and add annotated starts
    processed_genes = set()

    for parsed, name_info, original_line in existing_entries:
        ensembl_id = parsed["ensembl_id"]
        fields = parsed["fields"].copy()

        # Update gene name if needed
        ref_gene_name = gene_name_mapping[ensembl_id]
        if name_info["gene_name"].upper() != ref_gene_name.upper():
            # Reconstruct the name with updated gene name
            new_name = f"{ref_gene_name}_{name_info['gene_id']}_{name_info['start_type']}_{name_info['start_codon']}_{name_info['efficiency']}"
            fields[3] = new_name
            stats["updated_gene_names"] += 1

        # Check for duplicates
        position_key = f"{ensembl_id}_{fields[0]}_{fields[1]}_{fields[2]}"
        if position_key in seen_positions:
            stats["duplicates"] += 1
            continue

        seen_positions.add(position_key)
        valid_entries.append("\t".join(fields))
        stats["valid_alternatives"] += 1

        # Add annotated start for this gene if needed (same logic as before)
        if (
            add_annotated_starts
            and ensembl_id not in genes_with_existing_annotated
            and ensembl_id not in processed_genes
            and ensembl_id in gene_transcript_mapping
        ):
            # [Same annotated start addition logic as before]
            gene_transcripts = gene_transcript_mapping[ensembl_id]
            selected_transcript = None

            if preferred_transcripts:
                for transcript_info in gene_transcripts:
                    transcript_id = transcript_info["transcript_id"]
                    if transcript_id in preferred_transcripts:
                        selected_transcript = transcript_info
                        break
                    base_id = transcript_id.split(".")[0]
                    if any(base_id in pref for pref in preferred_transcripts):
                        selected_transcript = transcript_info
                        break

            if selected_transcript is None and gene_transcripts:
                selected_transcript = gene_transcripts[0]

            if selected_transcript:
                chrom = selected_transcript["chromosome"]
                strand = selected_transcript["strand"]

                if strand == "+":
                    start_pos = selected_transcript["start"]
                    bed_start = start_pos
                    bed_end = start_pos + 2
                else:
                    start_pos = selected_transcript["end"]
                    bed_start = start_pos - 2
                    bed_end = start_pos

                gene_name = gene_name_mapping[ensembl_id]
                annotated_name = f"{gene_name}_{ensembl_id}_Annotated_AUG_0.0"
                annotated_entry = (
                    f"{chrom}\t{bed_start}\t{bed_end}\t{annotated_name}\t0\t{strand}"
                )

                position_key = f"{ensembl_id}_{chrom}_{bed_start}_{bed_end}"
                if position_key not in seen_positions:
                    valid_entries.append(annotated_entry)
                    seen_positions.add(position_key)
                    added_annotated_starts.append(gene_name)
                    stats["annotated_starts_added"] += 1

                processed_genes.add(ensembl_id)

    # Write output
    output_bed_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_bed_path, "w") as f_out:
        f_out.write("\n".join(valid_entries))

    # Print summary
    if verbose:
        print(f"\nEnhanced BED Cleanup Summary:")
        print(f"  ├─ Total entries processed: {stats['total']}")

        if stats["parsing_errors"]:
            print(f"  ├─ Parsing errors (first 5):")
            for error in stats["parsing_errors"][:5]:
                print(f"  │  └─ {error}")

        print(
            f"  ├─ Invalid entries removed: {stats['invalid_format'] + stats['invalid_ensembl']}"
        )
        if filter_uorfs:
            print(f"  ├─ uORF entries filtered: {stats['uorfs_filtered']}")
        print(f"  ├─ Duplicates removed: {stats['duplicates']}")
        print(f"  ├─ Gene names updated: {stats['updated_gene_names']}")
        print(f"  ├─ Valid alternative entries: {stats['valid_alternatives']}")

        if add_annotated_starts:
            print(f"  ├─ Annotated starts added: {stats['annotated_starts_added']}")

        print(f"  ├─ Genes with alternatives: {len(stats['genes_with_alternatives'])}")
        print(f"  └─ Final entries in cleaned file: {len(valid_entries)}")

    # Convert sets to lists for JSON serialization
    stats["genes_with_alternatives"] = list(stats["genes_with_alternatives"])
    stats["genes_missing_annotated"] = list(stats["genes_missing_annotated"])
    stats["added_annotated_starts"] = added_annotated_starts

    return stats


def filter_genes_with_alternatives_only(
    input_bed_path: Union[str, Path],
    output_bed_path: Union[str, Path],
    verbose: bool = True,
) -> Dict:
    """Filter BED file to only include genes that have both Annotated and alternative starts.

    Args:
        input_bed_path: Path to input BED file
        output_bed_path: Path to save filtered BED file
        verbose: Print summary

    Returns:
        Dictionary with filtering statistics
    """
    input_bed_path, output_bed_path = Path(input_bed_path), Path(output_bed_path)

    # First pass: identify genes with each type of start
    genes_with_annotated = set()
    genes_with_alternatives = set()
    all_entries = []

    with open(input_bed_path, "r") as f:
        for line in f:
            if not line.strip():
                continue

            parsed = parse_bed_line(line)
            if not parsed:
                continue

            name_parts = parsed["fields"][3].split("_")
            if len(name_parts) >= 3:
                ensembl_id = parsed["ensembl_id"]
                start_type = name_parts[2]

                if start_type == "Annotated":
                    genes_with_annotated.add(ensembl_id)
                elif start_type in ["Truncated", "Extended"]:
                    genes_with_alternatives.add(ensembl_id)

                all_entries.append((ensembl_id, start_type, line.strip()))

    # Find genes with both annotated and alternative starts
    complete_genes = genes_with_annotated.intersection(genes_with_alternatives)

    # Second pass: keep only entries for complete genes
    filtered_entries = []
    for ensembl_id, start_type, line in all_entries:
        if ensembl_id in complete_genes:
            filtered_entries.append(line)

    # Write filtered output
    output_bed_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_bed_path, "w") as f_out:
        f_out.write("\n".join(filtered_entries))

    stats = {
        "total_genes_input": len(genes_with_annotated.union(genes_with_alternatives)),
        "genes_with_annotated": len(genes_with_annotated),
        "genes_with_alternatives": len(genes_with_alternatives),
        "complete_genes": len(complete_genes),
        "total_entries_input": len(all_entries),
        "total_entries_output": len(filtered_entries),
    }

    if verbose:
        print(f"\nGene Filtering Summary:")
        print(f"  ├─ Total genes in input: {stats['total_genes_input']}")
        print(f"  ├─ Genes with Annotated starts: {stats['genes_with_annotated']}")
        print(f"  ├─ Genes with alternative starts: {stats['genes_with_alternatives']}")
        print(f"  ├─ Complete genes (both types): {stats['complete_genes']}")
        print(f"  ├─ Entries input: {stats['total_entries_input']}")
        print(f"  └─ Entries output: {stats['total_entries_output']}")

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


def generate_empirical_preferred_transcripts(
    ribosome_bed_path: Union[str, Path],
    original_transcripts_path: Union[str, Path],
    gtf_path: Union[str, Path],
    output_path: Union[str, Path],
) -> None:
    """Generate preferred transcripts based on ribosome profiling data.

    For each gene in the ribosome profiling data, finds the GTF transcript that
    best matches the highest efficiency annotated start position. Only includes
    genes present in the ribosome profiling dataset.

    Args:
        ribosome_bed_path: Path to ribosome profiling BED file
        original_transcripts_path: Path to original preferred transcripts
        gtf_path: Path to GTF annotation file
        output_path: Path to save updated preferred transcripts
    """
    # Get genes from ribosome profiling data
    riboprof_genes = set()
    with open(ribosome_bed_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 4:
                name_info = parse_alternative_start_name(fields[3])
                gene_id = name_info["gene_id"].split(".")[0]
                riboprof_genes.add(gene_id)

    # Extract best annotated starts from ribosome profiling (only for genes with annotated starts)
    best_starts = {}  # gene_id -> (chrom, position, strand, efficiency)

    with open(ribosome_bed_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 6:
                continue

            chrom, start, end, name, score, strand = fields[:6]
            name_info = parse_alternative_start_name(name)

            if name_info["start_type"] != "Annotated":
                continue

            gene_id = name_info["gene_id"].split(".")[0]
            efficiency = name_info["efficiency"]
            start_pos = int(start) if strand == "+" else int(end)

            if gene_id not in best_starts or efficiency > best_starts[gene_id][3]:
                best_starts[gene_id] = (chrom, start_pos, strand, efficiency)

    # Load original transcript mappings (transcript -> gene)
    original_transcripts = set()
    with open(original_transcripts_path, "r") as f:
        original_transcripts = {line.strip() for line in f if line.strip()}

    transcript_to_gene = {}  # transcript_id -> gene_id
    gene_to_original_transcript = {}  # gene_id -> transcript_id

    # Map original transcripts to genes
    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#") or "\ttranscript\t" not in line:
                continue

            gene_match = re.search(r'gene_id "([^"]+)"', line)
            transcript_match = re.search(r'transcript_id "([^"]+)"', line)

            if gene_match and transcript_match:
                gene_id = gene_match.group(1).split(".")[0]
                transcript_id = transcript_match.group(1)
                transcript_to_gene[transcript_id] = gene_id

                if transcript_id in original_transcripts:
                    gene_to_original_transcript[gene_id] = transcript_id

    # Find best matching transcripts for empirical starts
    empirical_matches = {}  # gene_id -> transcript_id

    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#") or "\ttranscript\t" not in line:
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

            if not (gene_match and transcript_match):
                continue

            gene_id = gene_match.group(1).split(".")[0]
            transcript_id = transcript_match.group(1)

            if gene_id in best_starts:
                emp_chrom, emp_pos, emp_strand, _ = best_starts[gene_id]
                transcript_pos = start if strand == "+" else end

                if chrom == emp_chrom and strand == emp_strand:
                    distance = abs(transcript_pos - emp_pos)
                    if (
                        gene_id not in empirical_matches
                        or distance < empirical_matches[gene_id][1]
                    ):
                        empirical_matches[gene_id] = (transcript_id, distance)

    # Build final transcript list: one transcript per gene in ribosome profiling data
    final_transcripts = set()
    genes_changed = 0
    genes_kept_original = 0
    genes_no_transcript = 0

    for gene_id in riboprof_genes:
        chosen_transcript = None

        # Priority 1: Use empirical match if available
        if gene_id in empirical_matches:
            chosen_transcript = empirical_matches[gene_id][0]

            # Check if this differs from original
            if gene_id in gene_to_original_transcript:
                if chosen_transcript != gene_to_original_transcript[gene_id]:
                    genes_changed += 1
                else:
                    genes_kept_original += 1
            else:
                genes_changed += 1  # New transcript for gene not in original list

        # Priority 2: Use original transcript if gene has one
        elif gene_id in gene_to_original_transcript:
            chosen_transcript = gene_to_original_transcript[gene_id]
            genes_kept_original += 1

        # Priority 3: No transcript available
        else:
            genes_no_transcript += 1

        if chosen_transcript:
            final_transcripts.add(chosen_transcript)

    # Write output
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        for transcript_id in sorted(final_transcripts):
            f.write(f"{transcript_id}\n")

    print(
        f"Updated preferred transcripts for {len(riboprof_genes)} genes from ribosome profiling data:"
    )
    print(f"  ├─ Transcripts changed based on empirical data: {genes_changed}")
    print(f"  ├─ Transcripts kept from original list: {genes_kept_original}")
    print(f"  ├─ Genes without available transcripts: {genes_no_transcript}")
    print(f"  └─ Final transcript count: {len(final_transcripts)}")
