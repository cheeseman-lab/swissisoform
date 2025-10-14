#!/usr/bin/env python3
"""Analyze mutations in alternative isoform truncation regions.

This script batch processes a list of genes to identify mutations that occur
within alternative isoform truncation sites. It generates detailed statistical
analysis and visualizations of transcript-truncation pairs.

Arguments:
    gene_list (str): Path to file containing gene names.
    output_dir (str): Directory to save output files.
    --genome (str): Path to genome FASTA file.
    --annotation (str): Path to genome annotation GTF file.
    --bed (str): Path to alternative isoform BED file.
    --visualize (bool): Generate visualizations for each gene.
    --sources (List[str]): Mutation sources to query.
    --impact-types (List[str]): Mutation impact types to include.
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
from tqdm import tqdm
import sys

# Suppress pandas FutureWarning
warnings.filterwarnings("ignore", category=FutureWarning)

from swissisoform.genome import GenomeHandler
from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.mutations import MutationHandler
from swissisoform.utils import (
    parse_gene_list,
    save_gene_level_results,
    save_isoform_level_results,
    print_mutation_summary,
)

# Configure logger
logger = logging.getLogger(__name__)


class TqdmLoggingHandler(logging.Handler):
    """Logging handler that uses tqdm.write() to avoid interfering with progress bars.

    Only uses tqdm.write() when output is to a TTY (interactive terminal).
    For file/pipe output, uses regular stderr to avoid repeated progress bars.
    """

    def __init__(self, level=logging.NOTSET):
        """Initialize the TqdmLoggingHandler.

        Args:
            level: Logging level (default: logging.NOTSET).
        """
        super().__init__(level)

    def emit(self, record):
        """Emit a log record using tqdm.write() or stderr.

        Args:
            record: LogRecord to emit.
        """
        try:
            msg = self.format(record)
            # Only use tqdm.write() if we're in an interactive terminal
            # Otherwise, write directly to stderr to avoid progress bar duplication in logs
            if sys.stderr.isatty():
                tqdm.write(msg, file=sys.stderr)
            else:
                sys.stderr.write(msg + "\n")
                sys.stderr.flush()
        except Exception:
            self.handleError(record)


async def main(
    gene_list_path: str,
    output_dir: str,
    genome_path: str,
    annotation_path: str,
    bed_path: str,
    visualize: bool = False,
    sources: List[str] = None,
    impact_types: List[str] = None,
    top_n_per_type: Optional[int] = None,
):
    """Main function to process genes for mutation analysis.

    Args:
        gene_list_path (str): Path to file containing gene names.
        output_dir (str): Directory to save output files.
        genome_path (str): Path to the genome FASTA file.
        annotation_path (str): Path to the genome annotation GTF file.
        bed_path (str): Path to the alternative isoform BED file.
        visualize (bool): Whether to generate visualizations.
        sources (Optional[List[str]]): List of mutation sources to query.
        impact_types (Optional[List[str]]): List of mutation impact types to include.
        top_n_per_type (Optional[int]): Number of top alternative start sites to keep per type per transcript (None = all).

    Returns:
        None

    Raises:
        Exception: If any gene fails processing, error is printed and gene is skipped.
    """
    start_time = datetime.now()
    logger.info(
        f"Starting mutation analysis at {start_time.strftime('%Y-%m-%d %H:%M:%S')}"
    )

    # Create output directory
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # Initialize handlers
    logger.info("Initializing components...")
    logger.info("  Loading genome data...")
    genome = GenomeHandler(genome_path, annotation_path)
    logger.info("  Loading alternative isoform data...")
    alt_isoforms = AlternativeIsoform()
    alt_isoforms.load_bed(bed_path)
    logger.info("  Initializing mutation handler...")
    mutation_handler = MutationHandler(genome_handler=genome)

    # Set default values and create impact_types dict for compatibility
    if sources is None:
        sources = ["clinvar"]
    if impact_types is None:
        impact_types = ["missense variant", "nonsense variant", "frameshift variant"]

    # Convert to the dict format expected by the mutation handler
    impact_types_dict = {sources[0]: impact_types}

    # Read gene list
    logger.info(f"Reading gene list from {gene_list_path}")
    gene_names = parse_gene_list(gene_list_path)

    total_genes = len(gene_names)
    logger.info(f"Starting mutation analysis of {total_genes} genes")
    logger.info(f"Configuration:")
    logger.info(f"  Sources: {', '.join(sources)}")
    logger.info(f"  Impact types: {', '.join(impact_types)}")
    logger.info(f"  Visualizations: {visualize}")

    # Process all genes
    results = []

    # Disable progress bar if output is not a TTY (e.g., piped to file)
    show_progress = sys.stderr.isatty()
    for idx, gene_name in enumerate(
        tqdm(
            gene_names,
            desc="Analyzing genes",
            unit="gene",
            disable=not show_progress,
            file=sys.stderr,
            leave=True,
        ),
        start=1,
    ):
        # Log progress for non-TTY output (e.g., when piped to file)
        if not show_progress:
            logger.info(f"Processing gene {idx}/{total_genes}: {gene_name}")

        result = await mutation_handler.analyze_gene_mutations_comprehensive(
            gene_name=gene_name,
            genome_handler=genome,
            alt_isoform_handler=alt_isoforms,
            output_dir=output_dir,
            visualize=visualize,
            impact_types=impact_types_dict,
            sources=sources,
            top_n_per_type_per_transcript=top_n_per_type,
        )
        results.append(result)

        # Save both levels of results
        save_gene_level_results(results, output_dir)
        save_isoform_level_results(results, output_dir)

    # Final summary
    end_time = datetime.now()
    duration = end_time - start_time

    # Create and print summary
    results_df = pd.DataFrame(results)
    print_mutation_summary(results_df, output_dir)

    logger.info(
        f"Analysis completed in {duration} at {end_time.strftime('%Y-%m-%d %H:%M:%S')}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze mutations in alternative isoform regions"
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
        default="../data/ribosome_profiling/isoforms_with_transcripts.bed",
        help="Path to alternative isoform BED file",
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
        default=[
            "missense variant",
            "nonsense variant",
            "frameshift variant",
        ],
        help="Mutation impact types to include (space-separated)",
    )
    parser.add_argument(
        "--top-n-per-type",
        type=int,
        default=None,
        help="Number of top alternative start sites to keep per type (Truncated/Extended) per transcript (default: None = all)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="count",
        default=0,
        help="Increase verbosity (use -v, -vv, or -vvv for more detail)",
    )

    args = parser.parse_args()

    # Configure logging based on verbosity
    if args.verbose == 0:
        log_level = logging.WARNING
    elif args.verbose == 1:
        log_level = logging.INFO
    else:
        log_level = logging.DEBUG

    # Configure root logger with TqdmLoggingHandler to avoid conflicts with progress bars
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)

    # Clear any existing handlers
    root_logger.handlers.clear()

    # Add TqdmLoggingHandler
    handler = TqdmLoggingHandler()
    handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s - %(levelname)s - %(message)s",
            datefmt="%H:%M:%S",
        )
    )
    root_logger.addHandler(handler)

    # Suppress verbose logging from external libraries
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    logging.getLogger("requests").setLevel(logging.WARNING)
    logging.getLogger("httpx").setLevel(logging.WARNING)
    logging.getLogger("httpcore").setLevel(logging.WARNING)

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
            visualize=args.visualize,
            sources=args.sources,
            impact_types=args.impact_types,
            top_n_per_type=args.top_n_per_type,
        )
    )
