#!/usr/bin/env python3
import asyncio
import pandas as pd
from pathlib import Path

from swissisoform.genome import GenomeHandler
from swissisoform.isoform import AlternativeIsoform
from swissisoform.translation import TruncatedProteinGenerator


async def generate_amino_acid_sequences(
    gene_list_path: str,
    output_dir: str,
    genome_path: str = '../data/genome_data/hg38.fa',
    annotation_path: str = '../data/genome_data/hg38.ncbiRefSeq.gtf',
    bed_path: str = '../data/ribosome_profiling/RiboTISHV6_MD2025_AnnoToTruncation_exonintersect.bed',
    include_canonical: bool = True,
    min_length: int = 50,
    max_length: int = 1000,
    output_format: str = 'fasta,csv'
):
    """
    Generate amino acid sequences from truncated transcripts for a list of genes.
    
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
    """
    print("Initializing components:")
    
    print("  ├─ Loading genome data...")
    genome = GenomeHandler(genome_path, annotation_path)
    
    print("  ├─ Loading alternative isoform data...")
    alt_isoforms = AlternativeIsoform()
    alt_isoforms.load_bed(bed_path)
    
    print("  └─ Initializing protein generator...")
    protein_generator = TruncatedProteinGenerator(
        genome_handler=genome,
        alt_isoform_handler=alt_isoforms,
        output_dir=output_dir
    )
    
    # Read gene list
    print(f"\nReading gene list from {gene_list_path}")
    with open(gene_list_path, 'r') as f:
        gene_names = [line.strip() for line in f if line.strip()]
    
    total_genes = len(gene_names)
    print(f"Found {total_genes} genes to process")
    
    # Generate dataset
    print("\nGenerating amino acid sequences for deep learning dataset")
    dataset = protein_generator.prepare_deep_learning_dataset(
        gene_list=gene_names,
        output_format=output_format,
        include_canonical=include_canonical,
        min_length=min_length,
        max_length=max_length
    )
    
    # Print summary
    print("\nDataset summary:")
    print(f"  ├─ Total sequences: {len(dataset)}")
    
    if not dataset.empty:
        truncated_count = dataset[dataset['is_truncated'] == 1].shape[0]
        canonical_count = dataset[dataset['is_truncated'] == 0].shape[0]
        
        print(f"  ├─ Truncated sequences: {truncated_count}")
        print(f"  ├─ Canonical sequences: {canonical_count}")
        print(f"  ├─ Average sequence length: {dataset['length'].mean():.1f}")
        print(f"  ├─ Minimum sequence length: {dataset['length'].min()}")
        print(f"  ├─ Maximum sequence length: {dataset['length'].max()}")
        
        genes_with_data = dataset['gene'].nunique()
        print(f"  └─ Genes with valid sequences: {genes_with_data}/{total_genes}")
    else:
        print("  └─ No valid sequences generated")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate amino acid sequences from truncated transcripts')
    parser.add_argument('gene_list', help='Path to file containing gene names')
    parser.add_argument('output_dir', help='Directory to save output files')
    parser.add_argument('--genome', default='../data/genome_data/hg38.fa', help='Path to genome FASTA')
    parser.add_argument('--annotation', default='../data/genome_data/hg38.ncbiRefSeq.gtf', help='Path to genome annotation')
    parser.add_argument('--bed', default='../data/ribosome_profiling/RiboTISHV6_MD2025_AnnoToTruncation_exonintersect.bed', 
                        help='Path to alternative isoform BED file')
    parser.add_argument('--include-canonical', action='store_true', help='Include canonical sequences')
    parser.add_argument('--min-length', type=int, default=50, help='Minimum protein length')
    parser.add_argument('--max-length', type=int, default=1000, help='Maximum protein length')
    parser.add_argument('--format', default='fasta,csv', help='Output format: fasta, csv, or fasta,csv')
    
    args = parser.parse_args()
    
    asyncio.run(generate_amino_acid_sequences(
        gene_list_path=args.gene_list,
        output_dir=args.output_dir,
        genome_path=args.genome,
        annotation_path=args.annotation,
        bed_path=args.bed,
        include_canonical=args.include_canonical,
        min_length=args.min_length,
        max_length=args.max_length,
        output_format=args.format
    ))