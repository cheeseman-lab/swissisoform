#!/usr/bin/env python3
"""Unified protein sequence generator with two modes: pairs and mutations.

This script generates amino acid sequences from alternative transcript variants
with two distinct modes:

Modes:
    1. Pairs mode: Generate canonical + alternative protein sequence pairs.
    2. Mutations mode: Generate canonical + alternative + mutated protein sequences.

Arguments:
    gene_list (str): Path to file containing gene names.
    output_dir (str): Directory to save output files.
    --genome (str): Path to genome FASTA file.
    --annotation (str): Path to genome annotation GTF file.
    --bed (str): Path to alternative isoform BED file.
    --mutations (bool): Enable mutations mode.
    --impact-types (List[str], optional): Mutation impact types to include (only used with --mutations).
    --min-length (int): Minimum protein length to include.
    --max-length (int): Maximum protein length to include.
    --format (str): Output format: fasta, csv, or fasta,csv.
    --top-n-per-type (int): Number of top alternative start sites to keep per type per transcript.

Returns:
    None

Raises:
    Exception: If any gene fails processing, error is printed and gene is skipped.
"""

import asyncio
import pandas as pd
from datetime import datetime
from pathlib import Path
import argparse
import logging
import warnings
from typing import Optional, List

# Suppress pandas FutureWarning
warnings.filterwarnings("ignore", category=FutureWarning)

from swissisoform.genome import GenomeHandler
from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.mutations import MutationHandler
from swissisoform.translation import AlternativeProteinGenerator
from swissisoform.utils import (
    parse_gene_list,
    print_translation_summary,
)

# Configure logger
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


async def main(
    gene_list_path: str,
    output_dir: str,
    genome_path: str,
    annotation_path: str,
    bed_path: str,
    mutations_mode: bool = False,
    impact_types: List[str] = None,
    min_length: int = 10,
    max_length: int = 100000,
    output_format: str = "fasta,csv",
    top_n_per_type: int = 1,
):
    """Main function to process genes for protein sequence generation.

    Args:
        gene_list_path (str): Path to file containing gene names.
        output_dir (str): Directory to save output files.
        genome_path (str): Path to the genome FASTA file.
        annotation_path (str): Path to the genome annotation GTF file.
        bed_path (str): Path to the alternative isoform BED file.
        mutations_mode (bool): If True, generate with mutations; if False, generate pairs only.
        impact_types (Optional[List[str]]): Mutation impact types to include (only used in mutations mode).
        min_length (int): Minimum protein length to include.
        max_length (int): Maximum protein length to include.
        output_format (str): Format to save sequences ('fasta', 'csv', or 'fasta,csv').
        top_n_per_type (int): Number of top alternative start sites to keep per type per transcript.

    Returns:
        None

    Raises:
        Exception: If any gene fails processing, error is printed and gene is skipped.
    """
    start_time = datetime.now()
    print(
        f"Starting protein sequence generation at {start_time.strftime('%Y-%m-%d %H:%M:%S')}"
    )

    # Create output directory
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # Initialize handlers
    print("\nInitializing components:")
    print("  ├─ Loading genome data...")
    genome = GenomeHandler(genome_path, annotation_path)

    print("  ├─ Loading alternative isoform data...")
    alt_isoforms = AlternativeIsoform()
    alt_isoforms.load_bed(bed_path)

    # Initialize mutation handler only if needed
    mutation_handler = None
    if mutations_mode:
        print("  ├─ Initializing mutation handler...")
        mutation_handler = MutationHandler()

    print("  └─ Initializing protein generator...")
    protein_generator = AlternativeProteinGenerator(
        genome_handler=genome,
        alt_isoform_handler=alt_isoforms,
        output_dir=output_dir,
        mutation_handler=mutation_handler,
        debug=False,
    )

    # Set default impact types if not provided and mutations are requested
    if mutations_mode and impact_types is None:
        impact_types = ["missense variant", "5 prime UTR variant"]

    # Read gene list
    print(f"\nReading gene list from {gene_list_path}")
    gene_names = parse_gene_list(gene_list_path)

    total_genes = len(gene_names)
    mode_name = "Mutations" if mutations_mode else "Pairs"
    print(f"\nStarting protein sequence generation for {total_genes} genes")
    print(f"Configuration:")
    print(f"  ├─ Mode: {mode_name}")
    print(f"  ├─ Length range: {min_length}-{max_length} amino acids")
    if mutations_mode:
        print(f"  ├─ Impact types: {impact_types}")
    print(f"  └─ Output format: {output_format}")

    # Process genes with progress reporting and validation
    successful_genes = 0
    failed_genes = []
    validation_stats = {
        "total_pairs": 0,
        "valid_pairs": 0,
        "invalid_extensions": 0,
        "invalid_truncations": 0,
        "length_filtered": 0,
        "identical_sequences": 0,
    }

    for gene_idx, gene_name in enumerate(gene_names, 1):
        print(f"\nProcessing gene {gene_idx}/{total_genes}: {gene_name}")

        try:
            # Test if gene has any data before full processing
            test_pairs = protein_generator.extract_gene_proteins(
                gene_name,
                top_n_per_type,
            )
            if not test_pairs:
                print(f"  └─ ❌ No transcript-truncation pairs found")
                failed_genes.append(gene_name)
                continue

            print(f"  ├─ Found {len(test_pairs)} transcript-alternative pairs")

            # Validate the pairs
            valid_pairs = 0
            for pair in test_pairs:
                validation_stats["total_pairs"] += 1

                canonical_protein = pair["canonical"]["protein"]
                alternative_protein = pair["alternative"]["protein"]
                region_type = pair["region_type"]

                # Check for identical sequences
                if canonical_protein == alternative_protein:
                    validation_stats["identical_sequences"] += 1
                    continue

                # Check length constraints
                if not (min_length <= len(canonical_protein) <= max_length):
                    validation_stats["length_filtered"] += 1
                    continue
                if not (min_length <= len(alternative_protein) <= max_length):
                    validation_stats["length_filtered"] += 1
                    continue

                if protein_generator.validate_protein_pair(
                    canonical_protein,
                    alternative_protein,
                    region_type,
                    gene_name,
                    pair["transcript_id"],
                    verbose=True,  # Set to False to reduce output
                ):
                    valid_pairs += 1
                    validation_stats["valid_pairs"] += 1
                else:
                    # Update specific validation stats based on region type
                    if region_type == "extension":
                        validation_stats["invalid_extensions"] += 1
                    elif region_type == "truncation":
                        validation_stats["invalid_truncations"] += 1

            if valid_pairs == 0:
                print(f"  └─ ❌ No valid pairs after validation")
                failed_genes.append(gene_name)
                continue

            print(f"  ├─ Valid pairs after validation: {valid_pairs}/{len(test_pairs)}")

            # Additional info for mutations mode
            if mutations_mode:
                # Quick check for mutations in alternative regions
                alt_features = alt_isoforms.get_translation_features(gene_name)
                if not alt_features.empty:
                    print(
                        f"  ├─ Found {len(alt_features)} alternative regions for mutation analysis"
                    )

            print(f"  └─ ✅ Gene processed successfully")
            successful_genes += 1

        except Exception as e:
            print(f"  └─ ❌ Error processing gene: {str(e)}")
            failed_genes.append(gene_name)

    # Print validation summary
    print(f"\nValidation Summary:")
    print(f"  ├─ Total transcript-alternative pairs: {validation_stats['total_pairs']}")
    print(f"  ├─ Valid pairs: {validation_stats['valid_pairs']}")
    print(f"  ├─ Invalid extensions: {validation_stats['invalid_extensions']}")
    print(f"  ├─ Invalid truncations: {validation_stats['invalid_truncations']}")
    print(f"  ├─ Length filtered: {validation_stats['length_filtered']}")
    print(f"  └─ Identical sequences: {validation_stats['identical_sequences']}")

    # Generate the full dataset
    print(f"\nGenerating final dataset...")

    if mutations_mode:
        dataset = (
            await protein_generator.create_protein_sequence_dataset_with_mutations(
                gene_list=gene_names,
                include_mutations=True,
                impact_types=impact_types,
                output_format=output_format,
                min_length=min_length,
                max_length=max_length,
            )
        )
    else:
        dataset = protein_generator.create_protein_sequence_dataset_pairs(
            gene_list=gene_names,
            output_format=output_format,
            min_length=min_length,
            max_length=max_length,
        )

    # Final summary
    end_time = datetime.now()
    duration = end_time - start_time

    # Create and print summary using the utils function
    print_translation_summary(
        dataset=dataset,
        successful_genes=successful_genes,
        total_genes=total_genes,
        failed_genes=failed_genes,
        mutations_mode=mutations_mode,
        output_dir=output_dir,
    )

    print(
        f"\n  └─ Analysis completed in {duration} at {end_time.strftime('%Y-%m-%d %H:%M:%S')}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate protein sequences in two modes: pairs or mutations"
    )
    parser.add_argument("gene_list", help="Path to file containing gene names")
    parser.add_argument("output_dir", help="Directory to save output files")
    parser.add_argument(
        "--genome",
        default="../data/genome_data/GRCh38.p7.genome.fa",
        help="Path to genome FASTA file",
    )
    parser.add_argument(
        "--annotation",
        default="../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf",
        help="Path to genome annotation GTF file",
    )
    parser.add_argument(
        "--bed",
        default="../data/ribosome_profiling/isoforms_with_transcripts.bed",
        help="Path to alternative isoform BED file",
    )
    parser.add_argument(
        "--mutations",
        action="store_true",
        help="Enable mutations mode (generate canonical + truncated + mutated sequences)",
    )
    parser.add_argument(
        "--impact-types",
        nargs="+",
        default=None,
        help="Mutation impact types to include (space-separated, only used with --mutations)",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=10,
        help="Minimum protein length to include",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=100000,
        help="Maximum protein length to include",
    )
    parser.add_argument(
        "--format",
        default="fasta,csv",
        help="Output format: fasta, csv, or fasta,csv",
    )
    parser.add_argument(
        "--top-n-per-type",
        type=int,
        default=1,
        help="Number of top alternative start sites to keep per type (Truncated/Extended) per transcript",
    )

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    asyncio.run(
        main(
            gene_list_path=args.gene_list,
            output_dir=args.output_dir,
            genome_path=args.genome,
            annotation_path=args.annotation,
            bed_path=args.bed,
            mutations_mode=args.mutations,
            impact_types=args.impact_types,
            min_length=args.min_length,
            max_length=args.max_length,
            output_format=args.format,
            top_n_per_type=args.top_n_per_type,
        )
    )
