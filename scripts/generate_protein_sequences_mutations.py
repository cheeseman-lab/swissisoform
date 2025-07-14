#!/usr/bin/env python3
"""Generate protein sequences with integrated mutations from truncated transcripts.

This script generates amino acid sequences from truncated transcript variants
while incorporating mutations found in the truncation regions.
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
from swissisoform.mutations import MutationHandler
from swissisoform.translation_mutations import (
    ValidatedMutationIntegratedProteinGenerator,
)
from swissisoform.utils import (
    parse_gene_list,
    load_preferred_transcripts,
)

# Configure logger
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


async def process_gene_with_mutations(
    gene_name: str,
    protein_generator: ValidatedMutationIntegratedProteinGenerator,
    preferred_transcripts: set = None,
    include_mutations: bool = True,
    impact_types: List[str] = None,
    min_length: int = 10,
    max_length: int = 100000,
    debug: bool = False,
) -> dict:
    """Process a single gene to generate protein sequences with mutations.

    Args:
        gene_name: Name of the gene
        protein_generator: Initialized MutationIntegratedProteinGenerator instance
        preferred_transcripts: Optional set of transcript IDs to prioritize
        include_mutations: Whether to include mutation variants
        impact_types: List of mutation impact types to include
        min_length: Minimum protein length to include
        max_length: Maximum protein length to include

    Returns:
        Dictionary containing analysis results
    """
    try:
        print(f"  ├─ Processing gene {gene_name}...")

        # Set debug mode for detailed output
        protein_generator.debug = debug

        # Extract enhanced protein pairs with validated mutations
        enhanced_pairs = (
            await protein_generator.extract_gene_proteins_with_validated_mutations(
                gene_name, preferred_transcripts, include_mutations, impact_types
            )
        )

        if not enhanced_pairs:
            return {
                "gene_name": gene_name,
                "status": "no_data",
                "error": "Could not extract proteins for gene",
            }

        # Process and validate sequences
        valid_pairs = 0
        skipped_pairs = 0
        total_mutations = 0
        all_sequences = []

        for pair in enhanced_pairs:
            canonical_protein = pair["canonical"]["protein"]
            truncated_protein = pair["truncated_base"]["protein"]
            transcript_id = pair["transcript_id"]
            truncation_id = pair["truncation_id"]

            # Check length constraints for base sequences
            if not (min_length <= len(canonical_protein) <= max_length):
                skipped_pairs += 1
                continue

            if not (min_length <= len(truncated_protein) <= max_length):
                skipped_pairs += 1
                continue

            # Check if sequences are different
            if truncated_protein == canonical_protein:
                skipped_pairs += 1
                continue

            # Add canonical sequence
            all_sequences.append(
                {
                    "variant_type": "canonical",
                    "sequence": canonical_protein,
                    "length": len(canonical_protein),
                }
            )

            # Add base truncated sequence
            all_sequences.append(
                {
                    "variant_type": "truncated",
                    "sequence": truncated_protein,
                    "length": len(truncated_protein),
                }
            )

            # Process mutation variants
            mutation_variants = pair.get("truncated_mutations", [])
            pair_mutations = 0

            print(
                f"  │  │  └─ DEBUG: Found {len(mutation_variants)} mutation variants from enhanced_pairs"
            )

            for mut_variant in mutation_variants:
                mutated_protein = mut_variant["protein"]

                # Check length constraints
                if min_length <= len(mutated_protein) <= max_length:
                    # Check if mutation actually changed the protein
                    if mutated_protein != truncated_protein:
                        all_sequences.append(
                            {
                                "variant_type": "truncated_mutated",
                                "sequence": mutated_protein,
                                "length": len(mutated_protein),
                                "mutation": mut_variant["mutation"],
                            }
                        )
                        pair_mutations += 1
                        print(
                            f"  │  │  └─ DEBUG: Added mutation variant {pair_mutations}"
                        )

            total_mutations += pair_mutations
            valid_pairs += 1

            print(
                f"  │  ├─ {transcript_id} × {truncation_id}: "
                f"{pair_mutations} mutation variants generated"
            )

        if not all_sequences:
            return {
                "gene_name": gene_name,
                "status": "length_filter",
                "error": f"No pairs passed filters (skipped {skipped_pairs} pairs)",
            }

        # Calculate summary statistics
        canonical_count = len(
            [s for s in all_sequences if s["variant_type"] == "canonical"]
        )
        truncated_count = len(
            [s for s in all_sequences if s["variant_type"] == "truncated"]
        )
        mutated_count = len(
            [s for s in all_sequences if s["variant_type"] == "truncated_mutated"]
        )

        avg_canonical_length = sum(
            s["length"] for s in all_sequences if s["variant_type"] == "canonical"
        ) / max(canonical_count, 1)
        avg_truncated_length = sum(
            s["length"] for s in all_sequences if s["variant_type"] == "truncated"
        ) / max(truncated_count, 1)
        avg_mutated_length = sum(
            s["length"]
            for s in all_sequences
            if s["variant_type"] == "truncated_mutated"
        ) / max(mutated_count, 1)

        print(f"  │  ├─ Generated {valid_pairs} valid transcript-truncation pairs")
        print(
            f"  │  ├─ Canonical: {canonical_count} sequences (avg length: {avg_canonical_length:.0f} aa)"
        )
        print(
            f"  │  ├─ Truncated: {truncated_count} sequences (avg length: {avg_truncated_length:.0f} aa)"
        )
        print(
            f"  │  ├─ Mutated: {mutated_count} sequences (avg length: {avg_mutated_length:.0f} aa)"
        )
        if skipped_pairs > 0:
            print(f"  │  └─ Skipped {skipped_pairs} pairs due to filters")

        return {
            "gene_name": gene_name,
            "status": "success",
            "valid_pairs": valid_pairs,
            "skipped_pairs": skipped_pairs,
            "canonical_sequences": canonical_count,
            "truncated_sequences": truncated_count,
            "mutated_sequences": mutated_count,
            "total_mutations": total_mutations,
            "total_sequences": len(all_sequences),
            "avg_canonical_length": avg_canonical_length,
            "avg_truncated_length": avg_truncated_length,
            "avg_mutated_length": avg_mutated_length,
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
    include_mutations: bool = True,
    impact_types: List[str] = None,
    min_length: int = 10,
    max_length: int = 100000,
    output_format: str = "fasta,csv",
):
    """Main function to process genes for mutation-integrated protein sequence generation.

    Args:
        gene_list_path: Path to file containing gene names
        output_dir: Directory to save output files
        genome_path: Path to the genome FASTA file
        annotation_path: Path to the genome annotation GTF file
        bed_path: Path to the alternative isoform BED file
        preferred_transcripts_path: Optional path to file with preferred transcript IDs
        include_mutations: Whether to include mutation variants
        impact_types: List of mutation impact types to include
        min_length: Minimum protein length to include
        max_length: Maximum protein length to include
        output_format: Format to save sequences ('fasta', 'csv', or 'fasta,csv')
    """
    start_time = datetime.now()
    print(
        f"Starting mutation-integrated protein sequence generation at {start_time.strftime('%Y-%m-%d %H:%M:%S')}"
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
    print("  └─ Initializing validated mutation-integrated protein generator...")
    protein_generator = ValidatedMutationIntegratedProteinGenerator(
        genome_handler=genome,
        alt_isoform_handler=alt_isoforms,
        output_dir=output_dir,
        mutation_handler=mutation_handler,
        debug=False,  # Can be enabled for detailed debugging
    )

    # Load preferred transcripts if provided
    preferred_transcripts = None
    if preferred_transcripts_path:
        print(f"\nLoading preferred transcripts from {preferred_transcripts_path}")
        preferred_transcripts = load_preferred_transcripts(preferred_transcripts_path)
        print(f"Loaded {len(preferred_transcripts)} preferred transcript IDs")

    # Set default impact types if not provided
    if impact_types is None:
        impact_types = ["missense variant", "nonsense variant", "frameshift variant"]

    # Debug: Print what impact types we actually received
    print(f"Received impact types: {impact_types}")

    # Fix common issues with impact type parsing
    if len(impact_types) > 3 and any(" " not in item for item in impact_types):
        # Likely the arguments got split incorrectly
        # Try to reconstruct proper impact types
        reconstructed = []
        i = 0
        while i < len(impact_types):
            if (
                impact_types[i] == "missense"
                and i + 1 < len(impact_types)
                and impact_types[i + 1] == "variant"
            ):
                reconstructed.append("missense variant")
                i += 2
            elif (
                impact_types[i] == "nonsense"
                and i + 1 < len(impact_types)
                and impact_types[i + 1] == "variant"
            ):
                reconstructed.append("nonsense variant")
                i += 2
            elif (
                impact_types[i] == "frameshift"
                and i + 1 < len(impact_types)
                and impact_types[i + 1] == "variant"
            ):
                reconstructed.append("frameshift variant")
                i += 2
            else:
                reconstructed.append(impact_types[i])
                i += 1
        impact_types = reconstructed
        print(f"Reconstructed impact types: {impact_types}")

    # Read gene list
    print(f"\nReading gene list from {gene_list_path}")
    gene_names = parse_gene_list(gene_list_path)
    total_genes = len(gene_names)

    print(
        f"\nStarting mutation-integrated protein sequence generation for {total_genes} genes"
    )
    print(f"Parameters: min_length={min_length}, max_length={max_length}")
    print(f"Include mutations: {include_mutations}")
    if include_mutations:
        print(f"Impact types: {impact_types}")
    if preferred_transcripts:
        print(
            f"Using {len(preferred_transcripts)} preferred transcripts when available"
        )

    # Process all genes
    results = []
    all_sequences = []

    for idx, gene_name in enumerate(gene_names, 1):
        print(f"\nProcessing gene {idx}/{total_genes}: {gene_name}")
        result = await process_gene_with_mutations(
            gene_name=gene_name,
            protein_generator=protein_generator,
            preferred_transcripts=preferred_transcripts,
            include_mutations=include_mutations,
            impact_types=impact_types,
            min_length=min_length,
            max_length=max_length,
            debug=True,  # Set to True for detailed debugging
        )
        results.append(result)

        # Collect sequences from successful results
        if result["status"] == "success" and "sequences" in result:
            # Add gene and metadata to each sequence
            for seq in result["sequences"]:
                seq.update(
                    {
                        "gene": gene_name,
                        "transcript_id": f"{gene_name}_transcript",  # Simplified
                        "truncation_id": f"{gene_name}_truncation",  # Simplified
                    }
                )
            all_sequences.extend(result["sequences"])

    # Create final dataset
    if all_sequences:
        dataset = pd.DataFrame(all_sequences)

        # Save results in the requested formats
        if "fasta" in output_format.lower():
            # Save as FASTA
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            from Bio import SeqIO

            records = []
            for idx, row in dataset.iterrows():
                record_id = f"{row['gene']}_{row['variant_type']}_{idx}"
                description = f"{row['variant_type']} protein"

                if row["variant_type"] == "truncated_mutated" and "mutation" in row:
                    mut_info = row["mutation"]
                    description += f" with mutation {mut_info['position']}_{mut_info['reference']}>{mut_info['alternate']}"

                records.append(
                    SeqRecord(
                        Seq(row["sequence"]), id=record_id, description=description
                    )
                )

            output_file = output_dir_path / "protein_sequences_with_mutations.fasta"
            SeqIO.write(records, output_file, "fasta")
            print(f"Saved dataset FASTA to {output_file}")

        if "csv" in output_format.lower():
            output_file = output_dir_path / "protein_sequences_with_mutations.csv"
            dataset.to_csv(output_file, index=False)
            print(f"Saved dataset CSV to {output_file}")

    # Final summary
    end_time = datetime.now()
    duration = end_time - start_time

    # Create and print summary
    results_df = pd.DataFrame(results)
    print(f"\nMutation-Integrated Protein Sequence Generation Summary:")
    print(f"  ├─ Total genes processed: {len(results_df)}")
    print("\n  ├─ Status breakdown:")
    for status, count in results_df["status"].value_counts().items():
        print(f"  │  ├─ {status}: {count}")

    if all_sequences:
        canonical_count = len(
            [s for s in all_sequences if s["variant_type"] == "canonical"]
        )
        truncated_count = len(
            [s for s in all_sequences if s["variant_type"] == "truncated"]
        )
        mutated_count = len(
            [s for s in all_sequences if s["variant_type"] == "truncated_mutated"]
        )

        print(f"\n  ├─ Dataset summary:")
        print(f"  │  ├─ Total sequences: {len(all_sequences)}")
        print(f"  │  ├─ Canonical sequences: {canonical_count}")
        print(f"  │  ├─ Truncated sequences: {truncated_count}")
        print(f"  │  ├─ Mutated truncated sequences: {mutated_count}")
        print(f"  │  ├─ Average sequence length: {dataset['length'].mean():.1f}")
        print(f"  │  ├─ Minimum sequence length: {dataset['length'].min()}")
        print(f"  │  └─ Maximum sequence length: {dataset['length'].max()}")

        genes_with_data = dataset["gene"].nunique()
        print(f"  │  └─ Genes with valid sequences: {genes_with_data}/{total_genes}")
    else:
        print("\n  └─ No valid sequences generated")

    print(
        f"\n  └─ Analysis completed in {duration} at {end_time.strftime('%Y-%m-%d %H:%M:%S')}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate amino acid sequences with mutations from truncated transcripts"
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
        "--include-mutations",
        action="store_true",
        default=True,
        help="Include mutation variants in output",
    )
    parser.add_argument(
        "--impact-types",
        nargs="+",
        default=None,  # Will be set in main function
        help="Mutation impact types to include (space-separated, use quotes for multi-word types)",
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
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable detailed debug output for mutation processing",
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
            include_mutations=args.include_mutations,
            impact_types=args.impact_types,
            min_length=args.min_length,
            max_length=args.max_length,
            output_format=args.format,
        )
    )
