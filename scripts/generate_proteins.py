#!/usr/bin/env python3
"""Unified protein sequence generator with two modes: pairs and mutations.

This script generates amino acid sequences from truncated transcript variants
with two distinct modes:
1. Pairs mode: Generate canonical + truncated protein sequence pairs
2. Mutations mode: Generate canonical + truncated + mutated protein sequences
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
from swissisoform.translation import TruncatedProteinGenerator
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
    mutations_mode: bool = False,
    impact_types: List[str] = None,
    min_length: int = 10,
    max_length: int = 100000,
    output_format: str = "fasta,csv",
):
    """Main function to process genes for protein sequence generation.

    Args:
        gene_list_path: Path to file containing gene names
        output_dir: Directory to save output files
        genome_path: Path to the genome FASTA file
        annotation_path: Path to the genome annotation GTF file
        bed_path: Path to the alternative isoform BED file
        preferred_transcripts_path: Optional path to file with preferred transcript IDs
        mutations_mode: If True, generate with mutations; if False, generate pairs only
        impact_types: Mutation impact types to include (only used in mutations mode)
        min_length: Minimum protein length to include
        max_length: Maximum protein length to include
        output_format: Format to save sequences ('fasta', 'csv', or 'fasta,csv')
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
    protein_generator = TruncatedProteinGenerator(
        genome_handler=genome,
        alt_isoform_handler=alt_isoforms,
        output_dir=output_dir,
        mutation_handler=mutation_handler,
        debug=False,
    )

    # Load preferred transcripts if provided
    preferred_transcripts = None
    if preferred_transcripts_path:
        print(f"\nLoading preferred transcripts from {preferred_transcripts_path}")
        preferred_transcripts = load_preferred_transcripts(preferred_transcripts_path)
        print(f"Loaded {len(preferred_transcripts)} preferred transcript IDs")

    # Set default impact types if not provided and mutations are requested
    if mutations_mode and impact_types is None:
        impact_types = ["missense variant"]

    # Read gene list
    print(f"\nReading gene list from {gene_list_path}")
    gene_names = parse_gene_list(gene_list_path)
    total_genes = len(gene_names)

    # Print configuration
    mode_name = "Mutations" if mutations_mode else "Pairs"
    print(f"\nStarting protein sequence generation for {total_genes} genes")
    print(f"Configuration:")
    print(f"  ├─ Mode: {mode_name}")
    print(f"  ├─ Length range: {min_length}-{max_length} amino acids")
    if mutations_mode:
        print(f"  ├─ Impact types: {impact_types}")
    print(f"  └─ Output format: {output_format}")

    if preferred_transcripts:
        print(f"Using {len(preferred_transcripts)} preferred transcripts when available")

    # Process genes with progress reporting
    successful_genes = 0
    failed_genes = []
    
    for gene_idx, gene_name in enumerate(gene_names, 1):
        print(f"\nProcessing gene {gene_idx}/{total_genes}: {gene_name}")
        
        try:
            if mutations_mode:
                # Test if gene has any data before full processing
                test_pairs = protein_generator.extract_gene_proteins(
                    gene_name, preferred_transcripts
                )
                if not test_pairs:
                    print(f"  └─ ❌ No transcript-truncation pairs found")
                    failed_genes.append(gene_name)
                    continue
                
                print(f"  ├─ Found {len(test_pairs)} transcript-truncation pairs")
                print(f"  └─ ✅ Gene processed successfully")
                successful_genes += 1
            else:
                # For pairs mode, just check basic extraction
                gene_pairs = protein_generator.extract_gene_proteins(
                    gene_name, preferred_transcripts
                )
                if not gene_pairs:
                    print(f"  └─ ❌ No transcript-truncation pairs found")
                    failed_genes.append(gene_name)
                    continue
                
                print(f"  ├─ Found {len(gene_pairs)} transcript-truncation pairs")
                print(f"  └─ ✅ Gene processed successfully")
                successful_genes += 1
                
        except Exception as e:
            print(f"  └─ ❌ Error processing gene: {str(e)}")
            failed_genes.append(gene_name)

    # Generate the full dataset
    print(f"\nGenerating final dataset...")
    
    if mutations_mode:
        dataset = await protein_generator.create_protein_sequence_dataset_with_mutations(
            gene_list=gene_names,
            preferred_transcripts=preferred_transcripts,
            include_mutations=True,
            impact_types=impact_types,
            output_format=output_format,
            min_length=min_length,
            max_length=max_length,
        )
    else:
        dataset = protein_generator.create_protein_sequence_dataset_pairs(
            gene_list=gene_names,
            preferred_transcripts=preferred_transcripts,
            output_format=output_format,
            min_length=min_length,
            max_length=max_length,
        )

    # Final summary
    end_time = datetime.now()
    duration = end_time - start_time
    
    print(f"\n🏁 FINAL SUMMARY")
    print("=" * 60)
    print(f"  ├─ Genes processed successfully: {successful_genes}/{total_genes}")
    if failed_genes:
        print(f"  ├─ Failed genes: {len(failed_genes)}")
        if len(failed_genes) <= 5:
            print(f"  │   └─ {', '.join(failed_genes)}")
        else:
            print(f"  │   └─ {', '.join(failed_genes[:5])}, ... and {len(failed_genes)-5} more")
    
    if not dataset.empty:
        print(f"  ├─ Total sequences generated: {len(dataset)}")
        if mutations_mode and 'variant_type' in dataset.columns:
            type_counts = dataset['variant_type'].value_counts()
            for variant_type, count in type_counts.items():
                print(f"  │   ├─ {variant_type}: {count}")
        print(f"  └─ Average sequence length: {dataset['length'].mean():.1f}")
    else:
        print(f"  └─ ❌ No sequences generated")

    print(
        f"\nCompleted in {duration} at {end_time.strftime('%Y-%m-%d %H:%M:%S')}"
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
        help="Path to genome FASTA file"
    )
    parser.add_argument(
        "--annotation",
        default="../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf",
        help="Path to genome annotation GTF file"
    )
    parser.add_argument(
        "--bed",
        default="../data/ribosome_profiling/truncations_cleaned.bed",
        help="Path to alternative isoform BED file"
    )
    parser.add_argument(
        "--preferred-transcripts",
        default="../data/genome_data/hela_top_transcript.txt",
        help="Path to file containing preferred transcript IDs"
    )
    parser.add_argument(
        "--mutations",
        action="store_true",
        help="Enable mutations mode (generate canonical + truncated + mutated sequences)"
    )
    parser.add_argument(
        "--impact-types",
        nargs="+",
        default=None,
        help="Mutation impact types to include (space-separated, only used with --mutations)"
    )
    parser.add_argument(
        "--min-length", 
        type=int, 
        default=10, 
        help="Minimum protein length to include"
    )
    parser.add_argument(
        "--max-length", 
        type=int, 
        default=100000, 
        help="Maximum protein length to include"
    )
    parser.add_argument(
        "--format", 
        default="fasta,csv", 
        help="Output format: fasta, csv, or fasta,csv"
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
            preferred_transcripts_path=args.preferred_transcripts,
            mutations_mode=args.mutations,
            impact_types=args.impact_types,
            min_length=args.min_length,
            max_length=args.max_length,
            output_format=args.format,
        )
    )