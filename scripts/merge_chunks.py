#!/usr/bin/env python3
"""Merge chunked mutation analysis results.

This script handles the merge step after parallel chunk processing:
1. Merges gene-level results from all chunks
2. Merges isoform-level results from all chunks

Note: MANE annotation and filtering are now done in the summary step (summarize_results.py),
not during mutation analysis. This keeps raw data pristine and allows re-running summaries
with different filters/annotations without re-analyzing mutations.

Usage:
    python merge_chunks.py \\
        --output-dir ../results/hela/gnomad/mutations \\
        --num-chunks 8
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def merge_chunks(output_dir: Path, num_chunks: int, result_type: str) -> pd.DataFrame:
    """Merge results from all chunks.

    Args:
        output_dir: Base output directory containing chunk subdirectories
        num_chunks: Number of chunks to merge
        result_type: Either 'gene_level_results' or 'isoform_level_results'

    Returns:
        Merged DataFrame
    """
    logger.info(f"Merging {result_type} from {num_chunks} chunks...")

    merged_df = None

    for i in range(1, num_chunks + 1):
        chunk_file = output_dir / f"chunk_{i}" / f"{result_type}.csv"

        if not chunk_file.exists():
            logger.warning(f"Chunk {i} {result_type} missing: {chunk_file}")
            continue

        chunk_df = pd.read_csv(chunk_file)
        logger.info(f"  ├─ Chunk {i}: {len(chunk_df)} rows")

        if merged_df is None:
            merged_df = chunk_df
        else:
            merged_df = pd.concat([merged_df, chunk_df], ignore_index=True)

    if merged_df is None:
        raise ValueError(f"No {result_type} files found to merge")

    logger.info(f"  └─ Total merged: {len(merged_df)} rows")
    return merged_df


def main():
    """Main function for merging mutation analysis results."""
    parser = argparse.ArgumentParser(description="Merge chunked mutation results")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Base output directory (e.g., ../results/hela/gnomad/mutations)",
    )
    parser.add_argument(
        "--num-chunks",
        type=int,
        default=8,
        help="Number of chunks to merge (default: 8)",
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    if not output_dir.exists():
        logger.error(f"Output directory does not exist: {output_dir}")
        sys.exit(1)

    # Merge gene-level results
    gene_df = merge_chunks(output_dir, args.num_chunks, "gene_level_results")
    gene_output = output_dir / "gene_level_results.csv"
    gene_df.to_csv(gene_output, index=False)
    logger.info(f"Saved merged gene-level results: {gene_output}")

    # Merge isoform-level results
    isoform_df = merge_chunks(output_dir, args.num_chunks, "isoform_level_results")
    isoform_output = output_dir / "isoform_level_results.csv"
    isoform_df.to_csv(isoform_output, index=False)
    logger.info(f"Saved merged isoform-level results: {isoform_output}")

    logger.info("Merge complete!")


if __name__ == "__main__":
    main()
