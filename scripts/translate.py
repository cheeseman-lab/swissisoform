#!/usr/bin/env python3
"""Unified protein sequence generator with optional mutation integration.

This script generates amino acid sequences from truncated transcript variants
with optional mutation integration, replacing both generate_protein_sequences.py
and generate_protein_sequences_mutations.py.
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
from swissisoform.translation_new import TruncatedProteinGenerator
from swissisoform.utils import (
    parse_gene_list,
    load_preferred_transcripts,
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
    preferred_transcripts_path: Optional[str] = None,
    include_canonical: bool = True,
    include_mutations: bool = False,
    impact_types: List[str] = None,
    pairs_only: bool = False,
    min_length: int = 10,
    max_length: int = 100000,
    output_format: str = "fasta,csv",
    debug: bool = False,
):
    """Main function to process genes for protein sequence generation.

    Args:
        gene_list_path: Path to file containing gene names
        output_dir: Directory to save output files
        genome_path: Path to the genome FASTA file
        annotation_path: Path to the genome annotation GTF file
        bed_path: Path to the alternative isoform BED file
        preferred_transcripts_path: Optional path to file with preferred transcript IDs
        include_canonical: Whether to include canonical sequences
        include_mutations: Whether to include mutation variants
        impact_types: Mutation impact types to include
        pairs_only: Only include canonical/truncated pairs where truncation affects the transcript
        min_length: Minimum protein length to include
        max_length: Maximum protein length to include
        output_format: Format to save sequences ('fasta', 'csv', or 'fasta,csv')
        debug: Enable debug mode for detailed output
    """
    start_time = datetime.now()
    print(f"Starting protein sequence generation at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

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
    if include_mutations:
        print("  ├─ Initializing mutation handler...")
        mutation_handler = MutationHandler()
    
    print("  └─ Initializing protein generator...")
    protein_generator = TruncatedProteinGenerator(
        genome_handler=genome,
        alt_isoform_handler=alt_isoforms,
        output_dir=output_dir,
        mutation_handler=mutation_handler,
        debug=debug
    )

    # Load preferred transcripts if provided
    preferred_transcripts = None
    if preferred_transcripts_path:
        print(f"\nLoading preferred transcripts from {preferred_transcripts_path}")
        preferred_transcripts = load_preferred_transcripts(preferred_transcripts_path)
        print(f"Loaded {len(preferred_transcripts)} preferred transcript IDs")

    # Set default impact types if not provided and mutations are requested
    if include_mutations and impact_types is None:
        impact_types = ["missense variant", "nonsense variant", "frameshift variant"]

    # Read gene list
    print(f"\nReading gene list from {gene_list_path}")
    gene_names = parse_gene_list(gene_list_path)
    total_genes = len(gene_names)

    # Print configuration
    print(f"\nStarting protein sequence generation for {total_genes} genes")
    print(f"Parameters:")
    print(f"  ├─ Length range: {min_length}-{max_length} amino acids")
    print(f"  ├─ Include canonical: {include_canonical}")
    print(f"  ├─ Include mutations: {include_mutations}")
    if include_mutations:
        print(f"  ├─ Impact types: {impact_types}")
    print(f"  ├─ Pairs only: {pairs_only}")
    print(f"  └─ Output format: {output_format}")
    
    if preferred_transcripts:
        print(f"Using {len(preferred_transcripts)} preferred transcripts when available")

    # Generate sequences based on the mode
    if pairs_only and not include_mutations:
        # Legacy pairs-only mode (for compatibility)
        dataset = protein_generator.create_protein_sequence_dataset_pairs(
            gene_list=gene_names,
            preferred_transcripts=preferred_transcripts,
            output_format=output_format,
            min_length=min_length,
            max_length=max_length,
        )
    elif include_mutations:
        # New unified mode with mutation integration
        dataset = await protein_generator.create_protein_sequence_dataset_with_mutations(
            gene_list=gene_names,
            preferred_transcripts=preferred_transcripts,
            include_mutations=include_mutations,
            impact_types=impact_types,
            output_format=output_format,
            min_length=min_length,
            max_length=max_length,
        )
    else:
        # Standard mode without mutations
        dataset = protein_generator.create_protein_sequence_dataset(
            gene_list=gene_names,
            preferred_transcripts=preferred_transcripts,
            output_format=output_format,
            include_canonical=include_canonical,
            min_length=min_length,
            max_length=max_length,
        )

    # Final summary
    end_time = datetime.now()
    duration = end_time - start_time

    print(f"\n  └─ Analysis completed in {duration} at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate amino acid sequences from truncated transcripts with optional mutation integration"
    )
    parser.add_argument("gene_list", help="Path to file containing gene names")
    parser.add_argument("output_dir", help="Directory to save output files")
    parser.add_argument(
        "--genome",
        default="../data/genome_data/GRCh38.p7.genome.fa",
        help="Path to genome FASTA",
    )
    parser.add_argument(
        "--annotation",
        default="../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf",
        help="Path to genome annotation",
    )
    parser.add_argument(
        "--bed",
        default="../data/ribosome_profiling/full_truncations_JL_cleaned.bed",
        help="Path to alternative isoform BED file",
    )
    parser.add_argument(
        "--preferred-transcripts",
        default="../data/genome_data/hela_top_transcript.txt",
        help="Path to file containing preferred transcript IDs",
    )
    parser.add_argument(
        "--include-canonical", 
        action="store_true", 
        default=True,
        help="Include canonical sequences in output"
    )
    parser.add_argument(
        "--include-mutations", 
        action="store_true", 
        help="Include mutation variants in output"
    )
    parser.add_argument(
        "--impact-types",
        nargs="+",
        default=None,
        help="Mutation impact types to include (space-separated)"
    )
    parser.add_argument(
        "--pairs-only",
        action="store_true",
        help="Only include canonical/truncated pairs where truncation affects the transcript"
    )
    parser.add_argument(
        "--min-length", type=int, default=10, help="Minimum protein length"
    )
    parser.add_argument(
        "--max-length", type=int, default=100000, help="Maximum protein length"
    )
    parser.add_argument(
        "--format", default="fasta,csv", help="Output format: fasta, csv, or fasta,csv"
    )
    parser.add_argument(
        "--debug", 
        action="store_true", 
        help="Enable detailed debug output"
    )

    args = parser.parse_args()

    # Validate arguments
    if args.include_mutations and args.pairs_only:
        print("Warning: --pairs-only is ignored when --include-mutations is specified")

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
            preferred_transcripts_path=args.preferred_transcripts,
            include_canonical=args.include_canonical,
            include_mutations=args.include_mutations,
            impact_types=args.impact_types,
            pairs_only=args.pairs_only,
            min_length=args.min_length,
            max_length=args.max_length,
            output_format=args.format,
            debug=args.debug,
        )
    )