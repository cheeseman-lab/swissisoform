#!/usr/bin/env python3
"""Cleanup ribosome profiling data - simplified for transcript-based data."""

import yaml
import logging
from pathlib import Path
from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.utils import (
    simple_transcript_based_cleanup,
    update_gencode_gene_names,
)

# Configure logging to show INFO messages
logging.basicConfig(level=logging.INFO, format="%(message)s")

print("=" * 60)
print("SwissIsoform Data Cleanup Pipeline")
print("=" * 60)

# Load dataset configuration
config_path = "../data/ribosome_profiling/dataset_config.yaml"
with open(config_path) as f:
    config = yaml.safe_load(f)

print(f"\nStep 0: Updating GTFs with GENCODE v47 gene names")
print("=" * 60)

genome_data_dir = Path("../data/genome_data")
v47_gtf = genome_data_dir / "gencode.v47.annotation.gtf"

if not v47_gtf.exists():
    print(f"ERROR: GENCODE v47 GTF not found: {v47_gtf}")
    print("Cannot proceed without v47 reference for gene name updates.")
    exit(1)

# Update v24 and v25 GTFs with v47 gene names
gtf_versions = {
    "v24": {
        "input": genome_data_dir / "gencode.v24.annotation.gtf",
        "output": genome_data_dir / "gencode.v24.annotation.v47names.gtf",
    },
    "v25": {
        "input": genome_data_dir / "gencode.v25.annotation.gtf",
        "output": genome_data_dir / "gencode.v25.annotation.v47names.gtf",
    },
}

for version, paths in gtf_versions.items():
    if not paths["input"].exists():
        print(f"WARNING: {version} GTF not found: {paths['input']}")
        continue

    # Check if updated GTF already exists
    if paths["output"].exists():
        print(f"✓ {version} GTF with v47 names already exists: {paths['output'].name}")
    else:
        print(f"Updating {version} GTF with v47 gene names...")
        stats = update_gencode_gene_names(
            input_gtf_path=paths["input"],
            output_gtf_path=paths["output"],
            reference_gtf_path=v47_gtf,
            verbose=True,
        )
        print(f"✓ {version} update complete: {stats['genes_updated']} genes updated\n")

print(f"\n{'=' * 60}")
print("Processing Datasets")
print(f"{'=' * 60}\n")

# Process each dataset
all_results = []

for dataset in config["datasets"]:
    dataset_name = dataset["name"]
    bed_file = dataset["bed_file"]
    source_gtf = dataset["source_gtf_path"]

    # Use v47-updated GTF if it exists, otherwise fall back to original
    source_gtf_path = Path(source_gtf)
    v47_updated_gtf = source_gtf_path.parent / f"{source_gtf_path.stem}.v47names.gtf"

    if v47_updated_gtf.exists():
        actual_gtf = str(v47_updated_gtf)
    else:
        actual_gtf = source_gtf

    # Run simplified cleanup (verbose output comes from utils.py)
    input_bed_path = f"../data/ribosome_profiling/{bed_file}"
    output_bed_path = (
        f"../data/ribosome_profiling/{dataset_name}_isoforms_with_transcripts.bed"
    )
    refseq_gtf = "../data/genome_data/GRCh38_latest_genomic.gtf"

    stats = simple_transcript_based_cleanup(
        input_bed_path=input_bed_path,
        gtf_path=actual_gtf,  # Use v47-updated GTF
        output_bed_path=output_bed_path,
        refseq_gtf_path=refseq_gtf,
        verbose=True,
    )

    # Generate gene list
    alt_isoforms = AlternativeIsoform()
    alt_isoforms.load_bed(output_bed_path)
    gene_list = alt_isoforms.get_gene_list()

    # Write gene list
    gene_list_path = f"../data/ribosome_profiling/{dataset_name}_isoforms_gene_list.txt"
    with open(gene_list_path, "w") as f:
        for gene in gene_list:
            f.write(gene + "\n")

    all_results.append({"dataset": dataset_name, "genes": len(gene_list), **stats})

# Final summary
print(f"\n{'=' * 60}")
print("Cleanup Completed Successfully")
print(f"{'=' * 60}")

print(f"\nProcessed {len(all_results)} datasets:")
for result in all_results:
    print(
        f"  {result['dataset']}: {result['genes']:,} genes, {result['final_entries']:,} entries"
    )

print(f"\nNext step: sbatch 2_analyze_mutations.sh")
