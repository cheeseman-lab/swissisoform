#!/usr/bin/env python3
"""Analyze mutations in alternative isoform truncation regions.

This script batch processes a list of genes to identify mutations that occur
within alternative isoform truncation sites. It generates detailed statistical
analysis and visualizations of transcript-truncation pairs.
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
from swissisoform.utils import (
    parse_gene_list,
    save_gene_level_results,
    save_truncation_level_results,
    print_analysis_summary,
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
    visualize: bool = False,
    sources: List[str] = None,
    impact_types: List[str] = None,
):
    """Main function to process genes for mutation analysis.

    Args:
        gene_list_path: Path to file containing gene names
        output_dir: Directory to save output files
        genome_path: Path to the genome FASTA file
        annotation_path: Path to the genome annotation GTF file
        bed_path: Path to the alternative isoform BED file
        preferred_transcripts_path: Optional path to file with preferred transcript IDs
        visualize: Whether to generate visualizations
        include_unfiltered: Whether to include unfiltered mutation analysis
        sources: List of mutation sources to query
        impact_types: List of mutation impact types to include
    """
    start_time = datetime.now()
    print(f"Starting mutation analysis at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

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
    print("  └─ Initializing mutation handler...")
    mutation_handler = MutationHandler()

    # Load preferred transcripts if provided
    preferred_transcripts = None
    if preferred_transcripts_path:
        print(f"\nLoading preferred transcripts from {preferred_transcripts_path}")
        preferred_transcripts = load_preferred_transcripts(preferred_transcripts_path)
        print(f"Loaded {len(preferred_transcripts)} preferred transcript IDs")

    # Set default values and create impact_types dict for compatibility
    if sources is None:
        sources = ["clinvar"]
    if impact_types is None:
        impact_types = ["missense variant", "nonsense variant", "frameshift variant"]

    # Convert to the dict format expected by the mutation handler
    impact_types_dict = {sources[0]: impact_types}

    # Read gene list
    print(f"\nReading gene list from {gene_list_path}")
    gene_names = parse_gene_list(gene_list_path)

    total_genes = len(gene_names)
    print(f"\nStarting mutation analysis of {total_genes} genes")
    print(f"Configuration:")
    print(f"  ├─ Sources: {', '.join(sources)}")
    print(f"  ├─ Impact types: {', '.join(impact_types)}")
    print(f"  ├─ Visualizations: {visualize}")
    print(f"  ├─ Include unfiltered: {include_unfiltered}")
    if preferred_transcripts:
        print(f"  └─ Using {len(preferred_transcripts)} preferred transcripts")

    # Process all genes
    results = []

    for idx, gene_name in enumerate(gene_names, 1):
        print(f"\nProcessing gene {idx}/{total_genes}: {gene_name}")
        result = await mutation_handler.analyze_gene_mutations_comprehensive(
            gene_name=gene_name,
            genome_handler=genome,
            alt_isoform_handler=alt_isoforms,
            output_dir=output_dir,
            visualize=visualize,
            impact_types=impact_types_dict,
            preferred_transcripts=preferred_transcripts,
        )
        results.append(result)

        # Save both levels of results
        save_gene_level_results(results, output_dir)
        save_truncation_level_results(results, output_dir)

    # Final summary
    end_time = datetime.now()
    duration = end_time - start_time

    # Create and print summary
    results_df = pd.DataFrame(results)
    print_analysis_summary(results_df, output_dir)

    print(
        f"  └─ Analysis completed in {duration} at {end_time.strftime('%Y-%m-%d %H:%M:%S')}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze mutations in alternative isoform truncation regions"
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
        "--visualize", action="store_true", help="Generate visualizations for each gene"
    )
    parser.add_argument(
        "--sources",
        nargs="+",
        default=["clinvar"],
        help="Mutation sources to query (space-separated)",
    )
    parser.add_argument(
        "--impact-types",
        nargs="+",
        default=["missense variant", "nonsense variant", "frameshift variant"],
        help="Mutation impact types to include (space-separated)",
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
            visualize=args.visualize,
            sources=args.sources,
            impact_types=args.impact_types,
        )
    )
