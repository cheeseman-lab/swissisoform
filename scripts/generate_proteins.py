#!/usr/bin/env python3
"""Fast protein sequence generator using pre-validated mutation IDs from step 2.

This script loads pre-validated variant IDs from step 2 results and passes them
to the existing AlternativeProteinGenerator methods to enable fast mode processing.

Arguments:
    gene_list (str): Path to file containing gene names.
    output_dir (str): Directory to save output files.
    --mutations-file (str): Path to isoform_level_results.csv from step 2.
    --genome (str): Path to genome FASTA file.
    --annotation (str): Path to genome annotation GTF file.
    --bed (str): Path to alternative isoform BED file.
    --min-length (int): Minimum protein length to include.
    --max-length (int): Maximum protein length to include.
    --format (str): Output format: fasta, csv, or fasta,csv.
    --fast-mode (bool): Enable fast mode (skip validation).
"""

import asyncio
import pandas as pd
from datetime import datetime
from pathlib import Path
import argparse
import logging
import warnings
from typing import Optional, List, Dict, Set

# Suppress pandas FutureWarning
warnings.filterwarnings("ignore", category=FutureWarning)

from swissisoform.genome import GenomeHandler
from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.mutations import MutationHandler
from swissisoform.translation import AlternativeProteinGenerator
from swissisoform.utils import (
    parse_gene_list,
    load_pre_validated_variants,
)

# Configure logger
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


async def main(
    gene_list_path: str,
    output_dir: str,
    mutations_file: str,
    genome_path: str,
    annotation_path: str,
    bed_path: str,
    sources: List[str] = None,
    impact_types: List[str] = None,
    min_length: int = 10,
    max_length: int = 100000,
    output_format: str = "fasta,csv",
    fast_mode: bool = True,
):
    """Main function for fast protein sequence generation.

    Args:
        gene_list_path (str): Path to file containing gene names
        output_dir (str): Directory to save output files
        mutations_file (str): Path to isoform_level_results.csv from step 2
        genome_path (str): Path to genome FASTA file
        annotation_path (str): Path to genome annotation GTF file
        bed_path (str): Path to alternative isoform BED file
        sources (List[str]): List of mutation sources to query
        impact_types (List[str]): List of impact types to include
        min_length (int): Minimum protein length to include
        max_length (int): Maximum protein length to include
        output_format (str): Output format specification
        fast_mode (bool): Enable fast mode

    Returns:
        None
    """
    # Set defaults
    if sources is None:
        sources = ["clinvar"]
    if impact_types is None:
        impact_types = ["missense variant", "nonsense variant", "frameshift variant"]

    start_time = datetime.now()
    print(
        f"Starting fast protein sequence generation at {start_time.strftime('%Y-%m-%d %H:%M:%S')}"
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

    # Load pre-validated variant IDs
    print(f"\nLoading pre-validated variant IDs from {mutations_file}...")
    pre_validated_variants = load_pre_validated_variants(mutations_file)

    # Read gene list
    print(f"\nReading gene list from {gene_list_path}")
    gene_names = parse_gene_list(gene_list_path)

    total_genes = len(gene_names)
    print(f"\nStarting fast protein sequence generation for {total_genes} genes")
    print(f"Configuration:")
    print(f"  ├─ Fast mode: {fast_mode} (skip validation)")
    print(f"  ├─ Pre-validated variants file: {mutations_file}")
    print(f"  ├─ Sources: {', '.join(sources)}")
    print(f"  ├─ Impact types: {', '.join(impact_types)}")
    print(f"  ├─ Length range: {min_length}-{max_length} amino acids")
    print(f"  └─ Output format: {output_format}")

    # Generate datasets
    print(f"\nGenerating datasets...")

    # Generate pairs dataset (canonical + alternative)
    print(f"\n1. Generating canonical + alternative pairs dataset...")
    pairs_dataset = protein_generator.create_protein_sequence_dataset_pairs(
        gene_list=gene_names,
        output_format=output_format,
        min_length=min_length,
        max_length=max_length,
    )

    # Generate mutations dataset using pre-validated variants
    print(f"\n2. Generating mutations dataset with pre-validated variants...")
    mutations_dataset = await protein_generator.create_protein_sequence_dataset_with_mutations(
        gene_list=gene_names,
        include_mutations=True,
        sources=sources,
        impact_types=impact_types,
        output_format=output_format,
        min_length=min_length,
        max_length=max_length,
        pre_validated_variants=pre_validated_variants,  # Pass pre-validated variants
        skip_validation=fast_mode,  # Enable fast mode
    )

    # Final summary
    end_time = datetime.now()
    duration = end_time - start_time

    print(f"\n🎉 Fast protein sequence generation completed!")
    print(f"  ├─ Duration: {duration}")
    print(f"  ├─ Pairs dataset: {len(pairs_dataset)} sequences")
    print(f"  ├─ Mutations dataset: {len(mutations_dataset)} sequences")
    print(f"  └─ Performance: ⚡ Used pre-validated variants (skipped validation)")

    print(f"\nGeneration completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate protein sequences using pre-validated mutations (fast mode)"
    )
    parser.add_argument("gene_list", help="Path to file containing gene names")
    parser.add_argument("output_dir", help="Directory to save output files")
    parser.add_argument(
        "--mutations-file",
        required=True,
        help="Path to isoform_level_results.csv from step 2",
    )
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
        "--sources",
        nargs="+",
        default=["clinvar"],
        help="Mutation sources to query (space-separated, default: clinvar)",
    )
    parser.add_argument(
        "--impact-types",
        nargs="+",
        default=["missense variant", "nonsense variant", "frameshift variant"],
        help="Mutation impact types to include (space-separated)",
    )
    parser.add_argument(
        "--format",
        default="fasta,csv",
        help="Output format: fasta, csv, or fasta,csv",
    )
    parser.add_argument(
        "--fast-mode",
        action="store_true",
        help="Enable fast mode using pre-validated mutations",
    )

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    asyncio.run(
        main(
            gene_list_path=args.gene_list,
            output_dir=args.output_dir,
            mutations_file=args.mutations_file,
            genome_path=args.genome,
            annotation_path=args.annotation,
            bed_path=args.bed,
            sources=args.sources,
            impact_types=args.impact_types,
            min_length=args.min_length,
            max_length=args.max_length,
            output_format=args.format,
            fast_mode=args.fast_mode,
        )
    )
