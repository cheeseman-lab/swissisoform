"""Utility functions for the SwissIsoform package.

This module provides various helper functions for mutation analysis,
file handling, and data processing used throughout the package.
"""

import pandas as pd
from typing import List, Dict, Optional, Union
import asyncio
import logging
from pathlib import Path
from biomart import BiomartServer


# Configure logger
logger = logging.getLogger(__name__)


def save_analysis_results(
    results: List[Dict],
    output_dir: Union[str, Path],
    filename: str = "analysis_results.csv",
) -> None:
    """Save analysis results to CSV file.

    Args:
        results: List of result dictionaries to save
        output_dir: Directory to save the results file
        filename: Name of the results file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results_df = pd.DataFrame(results)
    output_path = output_dir / filename

    results_df.to_csv(output_path, index=False)
    print(f"Results saved to {output_path}")


def validate_gene_list(
    gene_list_path: Union[str, Path], genome_handler: "GenomeHandler"
) -> List[str]:
    """Validate a list of genes against available genomic data.

    Args:
        gene_list_path: Path to file containing gene names
        genome_handler: Initialized GenomeHandler instance

    Returns:
        List of valid gene names
    """
    gene_list_path = Path(gene_list_path)

    if not gene_list_path.exists():
        raise FileNotFoundError(f"Gene list file not found: {gene_list_path}")

    # Read gene list
    with open(gene_list_path, "r") as f:
        gene_names = [line.strip() for line in f if line.strip()]

    # Validate genes
    valid_genes = []
    invalid_genes = []

    for gene_name in gene_names:
        features = genome_handler.find_gene_features(gene_name)
        if features.empty:
            invalid_genes.append(gene_name)
        else:
            valid_genes.append(gene_name)

    if invalid_genes:
        print(f"Warning: {len(invalid_genes)} genes not found in genome annotations:")
        print(
            ", ".join(invalid_genes[:10]) + (", ..." if len(invalid_genes) > 10 else "")
        )

    print(f"Found {len(valid_genes)} valid genes")
    return valid_genes


def parse_gene_list(gene_list_path: Union[str, Path]) -> List[str]:
    """Parse a file containing a list of gene names.

    Args:
        gene_list_path: Path to file containing gene names

    Returns:
        List of gene names
    """
    gene_list_path = Path(gene_list_path)

    if not gene_list_path.exists():
        raise FileNotFoundError(f"Gene list file not found: {gene_list_path}")

    with open(gene_list_path, "r") as f:
        gene_names = [line.strip() for line in f if line.strip()]

    print(f"Read {len(gene_names)} genes from {gene_list_path}")
    return gene_names


def get_ensembl_reference():
    """Fetch current ENSEMBL ID to gene name mappings using biomart.

    Returns:
        dict: Mapping of ENSEMBL IDs to gene names
    """
    server = BiomartServer("http://www.ensembl.org/biomart")
    database = server.databases["ENSEMBL_MART_ENSEMBL"]
    dataset = database.datasets["hsapiens_gene_ensembl"]

    attributes = ["ensembl_gene_id", "external_gene_name"]
    response = dataset.search({"attributes": attributes})

    reference_data = {}
    for line in response.iter_lines():
        ensembl_id, gene_name = line.decode("utf-8").split("\t")
        if ensembl_id and gene_name:
            reference_data[ensembl_id] = gene_name

    return reference_data


def parse_bed_line(line):
    """Parse a single BED file line to extract ENSEMBL ID and gene name."""
    fields = line.strip().split("\t")
    if len(fields) < 4:
        return None

    name_parts = fields[3].split("_")
    if len(name_parts) < 3:
        return None

    # Get base ENSEMBL ID without version number
    full_ensembl_id = name_parts[0]
    ensembl_id = full_ensembl_id.split(".")[0]
    gene_name = name_parts[1]

    return {
        "ensembl_id": ensembl_id,
        "full_ensembl_id": full_ensembl_id,
        "gene_name": gene_name,
        "fields": fields,
        "line": line.strip(),
    }


def process_bed_file(bed_file_path, reference_data, callback):
    """Process BED file with a custom callback function for each line."""
    stats = {
        "total": 0,
        "processed": 0,
        "mismatches": [],
        "invalid_format": [],
        "not_found": [],
    }

    with open(bed_file_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            if not line.strip():
                continue

            stats["total"] += 1
            parsed = parse_bed_line(line)

            if not parsed:
                stats["invalid_format"].append((line_num, line.strip()))
                continue

            if not parsed["ensembl_id"].startswith("ENSG"):
                stats["invalid_format"].append((line_num, line.strip()))
                continue

            if parsed["ensembl_id"] not in reference_data:
                stats["not_found"].append(
                    (line_num, parsed["ensembl_id"], parsed["gene_name"])
                )
                continue

            ref_gene_name = reference_data[parsed["ensembl_id"]]
            if ref_gene_name.upper() != parsed["gene_name"].upper():
                stats["mismatches"].append(
                    (line_num, parsed["ensembl_id"], parsed["gene_name"], ref_gene_name)
                )

            callback(parsed, ref_gene_name)
            stats["processed"] += 1

    return stats


def validate_ensembl_mappings(bed_file_path):
    """Validate ENSEMBL ID to gene name mappings in a BED file."""
    print("Fetching current Ensembl reference data...")
    reference_data = get_ensembl_reference()
    print(f"Retrieved {len(reference_data)} reference mappings")

    valid_entries = []

    def validation_callback(parsed, ref_gene_name):
        valid_entries.append(
            (parsed["ensembl_id"], parsed["gene_name"], parsed["line"])
        )

    stats = process_bed_file(bed_file_path, reference_data, validation_callback)
    return {**stats, "valid_entries": valid_entries}


def update_bed_with_reference_names(bed_file_path, output_path):
    """Create a new BED file with gene names updated to match Ensembl reference."""
    print("Fetching current Ensembl reference data...")
    reference_data = get_ensembl_reference()
    print(f"Retrieved {len(reference_data)} reference mappings")

    with open(output_path, "w") as f_out:

        def update_callback(parsed, ref_gene_name):
            fields = parsed["fields"].copy()
            name_parts = fields[3].split("_")
            name_parts[1] = ref_gene_name
            fields[3] = "_".join(name_parts)
            f_out.write("\t".join(fields) + "\n")

        stats = process_bed_file(bed_file_path, reference_data, update_callback)

    return stats


def print_validation_results(results):
    """Print validation results in a readable format."""
    print(f"\nValidation Results:")
    print(f"Total entries processed: {results['total']}")
    print(f"Successfully validated: {results['processed']}")

    if results["invalid_format"]:
        print(f"\nInvalid format entries ({len(results['invalid_format'])} found):")
        for line_num, line in results["invalid_format"]:
            print(f"Line {line_num}: {line}")

    if results["not_found"]:
        print(
            f"\nENSEMBL IDs not found in reference ({len(results['not_found'])} found):"
        )
        for line_num, ensembl_id, gene_name in results["not_found"]:
            print(f"Line {line_num}: {ensembl_id} ({gene_name})")

    if results["mismatches"]:
        print(f"\nReference mismatches ({len(results['mismatches'])} found):")
        for line_num, ensembl_id, bed_gene, ref_gene in results["mismatches"]:
            print(f"Line {line_num}: {ensembl_id}")
            print(f"  BED gene name: {bed_gene}")
            print(f"  Reference gene name: {ref_gene}\n")


def print_update_results(stats):
    """Print update operation results."""
    print(f"\nUpdate Results:")
    print(f"Total entries: {stats['total']}")
    print(f"Successfully updated: {stats['processed']}")
    print(f"Invalid format entries: {len(stats['invalid_format'])}")
    print(f"ENSEMBL IDs not found: {len(stats['not_found'])}")
    print(f"Gene names that were different: {len(stats['mismatches'])}")
