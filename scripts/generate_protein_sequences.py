#!/usr/bin/env python3
"""Generate protein sequences from truncated transcripts.

This script generates amino acid sequences from truncated transcript variants
for a list of genes, creating datasets suitable for deep learning
applications.
"""

import asyncio
import pandas as pd
from pathlib import Path
import argparse
import logging
from typing import List, Optional
from swissisoform.genome import GenomeHandler
from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.translation import TruncatedProteinGenerator
from swissisoform.utils import parse_gene_list

# Configure logger
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

async def generate_protein_sequences(
    gene_list_path: str,
    output_dir: str,
    genome_path: str = "../data/genome_data/hg38.fa",
    annotation_path: str = "../data/genome_data/hg38.ncbiRefSeq.gtf",
    bed_path: str = "../data/ribosome_profiling/RiboTISHV6_MD2025_AnnoToTruncation_exonintersect.bed",
    include_canonical: bool = True,
    min_length: int = 50,
    max_length: int = 1000,
    output_format: str = "fasta,csv",
    pairs_only: bool = False,
):
    """Generate amino acid sequences from truncated transcripts for a list of genes.
    
    Args:
        gene_list_path: Path to file containing gene names
        output_dir: Directory to save output files
        genome_path: Path to the genome FASTA file
        annotation_path: Path to the genome annotation GTF file
        bed_path: Path to the alternative isoform BED file
        include_canonical: Whether to include canonical transcripts in the dataset
        min_length: Minimum protein length to include
        max_length: Maximum protein length to include
        output_format: Format to save sequences ('fasta', 'csv', or 'fasta,csv')
        pairs_only: Only include WT and truncated pairs where truncation affects the transcript
    """
    print("Initializing components:")
    # Initialize required handlers
    print("  ├─ Loading genome data...")
    genome = GenomeHandler(genome_path, annotation_path)
    print("  ├─ Loading alternative isoform data...")
    alt_isoforms = AlternativeIsoform()
    alt_isoforms.load_bed(bed_path)
    print("  └─ Initializing protein generator...")
    protein_generator = TruncatedProteinGenerator(
        genome_handler=genome, alt_isoform_handler=alt_isoforms, output_dir=output_dir
    )
    
    # Read gene list
    print(f"\nReading gene list from {gene_list_path}")
    gene_names = parse_gene_list(gene_list_path)
    total_genes = len(gene_names)
    print(f"Found {total_genes} genes to process")
    
    if pairs_only:
        # Generate dataset with only paired WT and truncations that affect transcripts
        print("\nGenerating paired canonical and truncated sequences (truncations must affect canonical)")
        dataset = protein_generator.create_protein_sequence_dataset_pairs(
            gene_list=gene_names,
            output_format=output_format,
            min_length=min_length,
            max_length=max_length,
        )
    else:
        # Generate all sequences using the original method
        print("\nGenerating amino acid sequences for deep learning dataset")
        dataset = protein_generator.create_protein_sequence_dataset(
            gene_list=gene_names,
            output_format=output_format,
            include_canonical=include_canonical,
            min_length=min_length,
            max_length=max_length,
        )
    
    # Print summary
    print("\nDataset summary:")
    print(f"  ├─ Total sequences: {len(dataset)}")
    
    if not dataset.empty:
        truncated_count = dataset[dataset["is_truncated"] == 1].shape[0]
        canonical_count = dataset[dataset["is_truncated"] == 0].shape[0]
        
        print(f"  ├─ Truncated sequences: {truncated_count}")
        print(f"  ├─ Canonical sequences: {canonical_count}")
        print(f"  ├─ Average sequence length: {dataset['length'].mean():.1f}")
        print(f"  ├─ Minimum sequence length: {dataset['length'].min()}")
        print(f"  ├─ Maximum sequence length: {dataset['length'].max()}")
        
        genes_with_data = dataset["gene"].nunique()
        print(f"  └─ Genes with valid sequences: {genes_with_data}/{total_genes}")
        
        if pairs_only:
            # Count the number of transcript-truncation pairs
            transcript_pairs = dataset.groupby("transcript_id").filter(
                lambda x: (x["is_truncated"] == 0).any() and (x["is_truncated"] == 1).any()
            )
            unique_transcripts_with_pairs = transcript_pairs["transcript_id"].nunique()
            
            print(f"  ├─ Transcripts with paired WT/truncated sequences: {unique_transcripts_with_pairs}")
            print(f"  └─ Total pairs: {len(transcript_pairs) // 2}")  # Divide by 2 since each pair has WT and truncated
    else:
        print("  └─ No valid sequences generated")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate amino acid sequences from truncated transcripts"
    )
    parser.add_argument("gene_list", help="Path to file containing gene names")
    parser.add_argument("output_dir", help="Directory to save output files")
    parser.add_argument(
        "--genome", default="../data/genome_data/hg38.fa", help="Path to genome FASTA"
    )
    parser.add_argument(
        "--annotation",
        default="../data/genome_data/hg38.ncbiRefSeq.gtf",
        help="Path to genome annotation",
    )
    parser.add_argument(
        "--bed",
        default="../data/ribosome_profiling/full_truncations_JL_cleaned.bed",
        help="Path to alternative isoform BED file",
    )
    parser.add_argument(
        "--include-canonical", action="store_true", help="Include canonical sequences"
    )
    parser.add_argument(
        "--pairs-only", 
        action="store_true", 
        help="Only include canonical/truncated pairs where truncation affects the transcript"
    )
    parser.add_argument(
        "--min-length", type=int, default=0, help="Minimum protein length"
    )
    parser.add_argument(
        "--max-length", type=int, default=10000, help="Maximum protein length"
    )
    parser.add_argument(
        "--format", default="fasta,csv", help="Output format: fasta, csv, or fasta,csv"
    )
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    asyncio.run(
        generate_protein_sequences(
            gene_list_path=args.gene_list,
            output_dir=args.output_dir,
            genome_path=args.genome,
            annotation_path=args.annotation,
            bed_path=args.bed,
            include_canonical=args.include_canonical,
            min_length=args.min_length,
            max_length=args.max_length,
            output_format=args.format,
            pairs_only=args.pairs_only,
        )
    )