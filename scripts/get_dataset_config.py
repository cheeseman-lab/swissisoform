#!/usr/bin/env python3
"""Get dataset configuration from YAML config file.

Reads dataset configuration and outputs the required paths as space-separated values
for use in bash scripts.

Usage:
    python get_dataset_config.py <config_file> <dataset_name>

Output (space-separated):
    <bed_file> <gtf_file> <gene_list>
"""

import sys
from pathlib import Path

import yaml


def get_dataset_config(config_file: str, dataset_name: str) -> tuple:
    """Read dataset configuration from YAML file.

    Args:
        config_file: Path to YAML configuration file
        dataset_name: Name of dataset to look up

    Returns:
        Tuple of (bed_file, gtf_file, gene_list)

    Raises:
        SystemExit: If dataset not found in config
    """
    with open(config_file) as f:
        config = yaml.safe_load(f)

    # Find the dataset
    dataset_config = None
    for ds in config["datasets"]:
        if ds["name"] == dataset_name:
            dataset_config = ds
            break

    if not dataset_config:
        print(f"ERROR: Dataset {dataset_name} not found in config", file=sys.stderr)
        sys.exit(1)

    # Get BED and gene list paths (use standard naming convention)
    bed_file = (
        f"../data/ribosome_profiling/{dataset_name}_isoforms_with_transcripts.bed"
    )
    gene_list = f"../data/ribosome_profiling/{dataset_name}_isoforms_gene_list.txt"

    # Get GTF path - prefer v47names version if it exists
    source_gtf_path = Path(dataset_config["source_gtf_path"])
    v47_gtf = source_gtf_path.parent / f"{source_gtf_path.stem}.v47names.gtf"

    if v47_gtf.exists():
        gtf_file = str(v47_gtf)
    else:
        gtf_file = dataset_config["source_gtf_path"]

    return bed_file, gtf_file, gene_list


def main():
    """Main function."""
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    config_file = sys.argv[1]
    dataset_name = sys.argv[2]

    bed_file, gtf_file, gene_list = get_dataset_config(config_file, dataset_name)

    # Output space-separated for bash
    print(f"{bed_file} {gtf_file} {gene_list}")


if __name__ == "__main__":
    main()
