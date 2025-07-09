#!/usr/bin/env python3
"""Generate protein sequences from truncated transcripts.

This script generates amino acid sequences from truncated transcript variants
for a list of genes, creating datasets suitable for deep learning
applications.
"""

import asyncio
import pandas as pd
from datetime import datetime
from pathlib import Path
import argparse
import logging
import warnings
from typing import Optional, List, Dict

# Suppress pandas FutureWarning
warnings.filterwarnings("ignore", category=FutureWarning)

from swissisoform.genome import GenomeHandler
from swissisoform.alternative_isoforms import AlternativeIsoform
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


async def process_gene(
    gene_name: str,
    protein_generator: TruncatedProteinGenerator,
    preferred_transcripts: set = None,
    min_length: int = 10,
    max_length: int = 100000,
    pairs_only: bool = False,
    include_canonical: bool = True,
) -> dict:
    """Process a single gene to generate protein sequences.

    Args:
        gene_name: Name of the gene to process
        protein_generator: Initialized TruncatedProteinGenerator instance
        preferred_transcripts: Optional set of transcript IDs to prioritize
        min_length: Minimum protein length to include
        max_length: Maximum protein length to include
        pairs_only: Only include WT and truncated pairs where truncation affects the transcript
        include_canonical: Whether to include canonical transcripts

    Returns:
        Dictionary containing processing results
    """
    try:
        print(f"  ├─ Processing gene {gene_name}...")

        # Use the new extract_gene_proteins method which handles everything
        gene_result = protein_generator.extract_gene_proteins(
            gene_name, preferred_transcripts
        )

        if not gene_result:
            return {
                "gene_name": gene_name,
                "status": "no_data",
                "error": "Could not extract proteins for gene",
            }

        canonical_protein = gene_result['canonical']['protein']
        truncated_protein = gene_result['truncated']['protein']
        transcript_id = gene_result['transcript_id']

        # Check length constraints
        if not (min_length <= len(canonical_protein) <= max_length):
            return {
                "gene_name": gene_name,
                "status": "length_filter",
                "error": f"Canonical protein length {len(canonical_protein)} outside range {min_length}-{max_length}",
            }

        if not (min_length <= len(truncated_protein) <= max_length):
            return {
                "gene_name": gene_name,
                "status": "length_filter",
                "error": f"Truncated protein length {len(truncated_protein)} outside range {min_length}-{max_length}",
            }

        # Check if the sequences are actually different
        if truncated_protein == canonical_protein:
            return {
                "gene_name": gene_name,
                "status": "identical_sequences",
                "error": "Truncation results in identical protein",
            }

        # Create sequence entries
        all_sequences = []
        
        # Add canonical sequence
        if include_canonical or pairs_only:
            canonical_entry = {
                "gene": gene_name,
                "transcript_id": transcript_id,
                "variant_id": "canonical",
                "sequence": canonical_protein,
                "length": len(canonical_protein),
                "is_truncated": 0,
            }
            all_sequences.append(canonical_entry)

        # Add truncated sequence
        trunc_id = f"trunc_{gene_result['truncation']['start']}_{gene_result['truncation']['end']}"
        if "start_codon" in gene_result['truncation'] and not pd.isna(gene_result['truncation']["start_codon"]):
            trunc_id = f"trunc_{gene_result['truncation']['start_codon']}_{gene_result['truncation']['start']}_{gene_result['truncation']['end']}"

        truncated_entry = {
            "gene": gene_name,
            "transcript_id": transcript_id,
            "variant_id": trunc_id,
            "sequence": truncated_protein,
            "length": len(truncated_protein),
            "is_truncated": 1,
        }
        all_sequences.append(truncated_entry)

        print(f"  │  ├─ Generated pair: canonical ({len(canonical_protein)} aa) + {trunc_id} ({len(truncated_protein)} aa)")

        return {
            "gene_name": gene_name,
            "status": "success",
            "transcript_id": transcript_id,
            "canonical_length": len(canonical_protein),
            "truncated_length": len(truncated_protein),
            "sequences": all_sequences,
            "error": None,
        }

    except Exception as e:
        logger.error(f"Error processing gene {gene_name}: {str(e)}")
        print(f"  └─ Error: {str(e)}")
        return {"gene_name": gene_name, "status": "error", "error": str(e)}


async def main(
    gene_list_path: str,
    output_dir: str,
    genome_path: str,
    annotation_path: str,
    bed_path: str,
    preferred_transcripts_path: Optional[str] = None,
    include_canonical: bool = True,
    min_length: int = 10,
    max_length: int = 100000,
    output_format: str = "fasta,csv",
    pairs_only: bool = False,
):
    """Main function to process genes for protein sequence generation.

    Args:
        gene_list_path: Path to file containing gene names
        output_dir: Directory to save output files
        genome_path: Path to the genome FASTA file
        annotation_path: Path to the genome annotation GTF file
        bed_path: Path to the alternative isoform BED file
        preferred_transcripts_path: Optional path to file with preferred transcript IDs
        include_canonical: Whether to include canonical transcripts in the dataset
        min_length: Minimum protein length to include
        max_length: Maximum protein length to include
        output_format: Format to save sequences ('fasta', 'csv', or 'fasta,csv')
        pairs_only: Only include WT and truncated pairs where truncation affects the transcript
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
    print("  └─ Initializing protein generator...")
    protein_generator = TruncatedProteinGenerator(
        genome_handler=genome, alt_isoform_handler=alt_isoforms, output_dir=output_dir
    )

    # Load preferred transcripts if provided
    preferred_transcripts = None
    if preferred_transcripts_path:
        print(f"\nLoading preferred transcripts from {preferred_transcripts_path}")
        preferred_transcripts = load_preferred_transcripts(preferred_transcripts_path)
        print(f"Loaded {len(preferred_transcripts)} preferred transcript IDs")

    # Read gene list
    print(f"\nReading gene list from {gene_list_path}")
    gene_names = parse_gene_list(gene_list_path)
    total_genes = len(gene_names)
    print(f"\nStarting protein sequence generation for {total_genes} genes")
    print(f"Parameters: min_length={min_length}, max_length={max_length}, pairs_only={pairs_only}")
    if preferred_transcripts:
        print(f"Using {len(preferred_transcripts)} preferred transcripts when available")

    # Process all genes
    results = []
    all_sequences = []

    for idx, gene_name in enumerate(gene_names, 1):
        print(f"\nProcessing gene {idx}/{total_genes}: {gene_name}")
        result = await process_gene(
            gene_name=gene_name,
            protein_generator=protein_generator,
            preferred_transcripts=preferred_transcripts,
            min_length=min_length,
            max_length=max_length,
            pairs_only=pairs_only,
            include_canonical=include_canonical,
        )
        results.append(result)
        
        # Collect sequences from successful results
        if result["status"] == "success" and "sequences" in result:
            all_sequences.extend(result["sequences"])

    # Create final dataset
    dataset = pd.DataFrame(all_sequences)

    # Save results in the requested formats
    if not dataset.empty:
        if "fasta" in output_format.lower():
            # Save as FASTA
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            from Bio import SeqIO
            
            records = []
            for _, row in dataset.iterrows():
                record_id = f"{row['gene']}_{row['transcript_id']}_{row['variant_id']}"
                description = f"{'Truncated' if row['is_truncated'] else 'Canonical'} protein"
                
                records.append(
                    SeqRecord(Seq(row["sequence"]), id=record_id, description=description)
                )
            
            if pairs_only:
                output_file = output_dir_path / "protein_sequence_dataset_pairs.fasta"
            else:
                output_file = output_dir_path / "protein_sequence_dataset.fasta"
            
            SeqIO.write(records, output_file, "fasta")
            print(f"Saved dataset FASTA to {output_file}")

        if "csv" in output_format.lower():
            if pairs_only:
                output_file = output_dir_path / "protein_sequence_dataset_pairs.csv"
            else:
                output_file = output_dir_path / "protein_sequence_dataset.csv"
            
            dataset.to_csv(output_file, index=False)
            print(f"Saved dataset CSV to {output_file}")

    # Final summary
    end_time = datetime.now()
    duration = end_time - start_time

    # Create and print summary
    results_df = pd.DataFrame(results)
    print(f"\nProtein Sequence Generation Summary:")
    print(f"  ├─ Total genes processed: {len(results_df)}")
    print("\n  ├─ Status breakdown:")
    for status, count in results_df["status"].value_counts().items():
        print(f"  │  ├─ {status}: {count}")

    if not dataset.empty:
        truncated_count = dataset[dataset["is_truncated"] == 1].shape[0]
        canonical_count = dataset[dataset["is_truncated"] == 0].shape[0]

        print(f"\n  ├─ Dataset summary:")
        print(f"  │  ├─ Total sequences: {len(dataset)}")
        print(f"  │  ├─ Truncated sequences: {truncated_count}")
        print(f"  │  ├─ Canonical sequences: {canonical_count}")
        print(f"  │  ├─ Average sequence length: {dataset['length'].mean():.1f}")
        print(f"  │  ├─ Minimum sequence length: {dataset['length'].min()}")
        print(f"  │  ├─ Maximum sequence length: {dataset['length'].max()}")

        genes_with_data = dataset["gene"].nunique()
        print(f"  │  └─ Genes with valid sequences: {genes_with_data}/{total_genes}")

        if pairs_only:
            # Count the number of transcript-truncation pairs
            transcript_pairs = dataset.groupby("transcript_id").filter(
                lambda x: (x["is_truncated"] == 0).any()
                and (x["is_truncated"] == 1).any()
            )
            unique_transcripts_with_pairs = transcript_pairs["transcript_id"].nunique()

            print(f"  │  ├─ Transcripts with paired WT/truncated sequences: {unique_transcripts_with_pairs}")
            print(f"  │  └─ Total pairs: {len(transcript_pairs) // 2}")
    else:
        print("\n  └─ No valid sequences generated")

    print(f"\n  └─ Analysis completed in {duration} at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate amino acid sequences from truncated transcripts"
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
        "--include-canonical", action="store_true", help="Include canonical sequences"
    )
    parser.add_argument(
        "--pairs-only",
        action="store_true",
        help="Only include canonical/truncated pairs where truncation affects the transcript",
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
        "--preferred-transcripts",
        default="../data/genome_data/hela_top_transcript.txt",
        help="Path to file containing preferred transcript IDs",
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
            include_canonical=args.include_canonical,
            min_length=args.min_length,
            max_length=args.max_length,
            output_format=args.format,
            pairs_only=args.pairs_only,
        )
    )