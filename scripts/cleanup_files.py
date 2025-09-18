#!/usr/bin/env python3
"""Cleanup files for the ribosome profiling data.

This script cleans up the GTF file and BED files for the ribosome profiling data.
It updates the gene names in the GTF file to match the latest version of the GTF file.
It performs comprehensive BED cleanup with transcript mapping.
"""

from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.utils import (
    update_gencode_gene_names,
    comprehensive_cleanup_bed_with_transcripts,
    subset_gene_list,
)

# GTF: update gene names
input_gtf = "../data/genome_data/gencode.v25.annotation.gtf"
output_gtf = "../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf"
reference_gtf = "../data/genome_data/gencode.v47.annotation.gtf"
preferred_transcripts = "../data/genome_data/hela_top_transcript.txt"

print("=== SwissIsoform Data Cleanup Pipeline ===")
print()

# Step 1: Update gene names in GTF
print("Step 1: Updating GTF gene names...")
update_gencode_gene_names(
    input_gtf_path=input_gtf,
    output_gtf_path=output_gtf,
    reference_gtf_path=reference_gtf,
    verbose=True,
)

print("\n" + "=" * 60)

# Step 2: Comprehensive BED cleanup with transcript mapping
input_bed = "../data/ribosome_profiling/Ly_2024b_TableS2_formatted.bed"
enhanced_bed = "../data/ribosome_profiling/isoforms_with_transcripts.bed"

print("Step 2: Comprehensive BED cleanup with transcript mapping...")
cleanup_stats = comprehensive_cleanup_bed_with_transcripts(
    input_bed_path=input_bed,
    gtf_path=output_gtf,
    preferred_transcripts_path=preferred_transcripts,
    output_bed_path=enhanced_bed,
    verbose=True,
)

print("\n" + "=" * 60)

# Step 3: Generate gene lists
print("Step 3: Generating gene lists...")

# Load the enhanced BED file to get gene list
alt_isoforms = AlternativeIsoform()
alt_isoforms.load_bed(enhanced_bed)
gene_list = alt_isoforms.get_gene_list()

# Write full gene list
gene_list_path = "../data/ribosome_profiling/isoforms_gene_list.txt"
with open(gene_list_path, "w") as f:
    for gene in gene_list:
        f.write(gene + "\n")

print(f"  ├─ Full gene list: {len(gene_list)} genes → {gene_list_path}")

# Generate reduced gene list
reduced_gene_list = subset_gene_list(gene_list_path=gene_list_path)

reduced_gene_list_path = "../data/ribosome_profiling/isoforms_gene_list_reduced.txt"
with open(reduced_gene_list_path, "w") as f:
    for gene in reduced_gene_list:
        f.write(gene + "\n")

print(
    f"  └─ Reduced gene list: {len(reduced_gene_list)} genes → {reduced_gene_list_path}"
)

print("\n" + "=" * 60)

# Step 4: Final summary
print("Final Summary:")
print(f"  ├─ Input entries: {cleanup_stats['total_entries']}")
print(f"  ├─ Valid entries: {cleanup_stats['valid_entries']}")
print(f"  ├─ Enhanced entries: {cleanup_stats['enhanced_entries']}")
print(f"  ├─ Final genes: {cleanup_stats['final_genes']}")
print(f"  ├─ Added annotated starts: {cleanup_stats['added_annotated']}")
print(f"  ├─ Transcript assignments: {cleanup_stats['transcript_assignments']}")
print(f"  └─ Reduced analysis set: {len(reduced_gene_list)} genes")

print(f"\nOutput files:")
print(f"  ├─ Enhanced GTF: {output_gtf}")
print(f"  ├─ Enhanced BED with transcripts: {enhanced_bed}")
print(f"  ├─ Full gene list: {gene_list_path}")
print(f"  └─ Reduced gene list: {reduced_gene_list_path}")

print(f"\nReady for analysis!")
print(f"  ├─ Each BED entry now has transcript mapping")
print(f"  ├─ Biologically-relevant canonical start selection")
print(f"  ├─ Complete genes only (annotated + alternatives)")
print(f"  └─ No separate transcript files needed")

# Optional: Run quality check
print(f"\nOptional quality check:")
print(f"  Run: alt_isoforms.check_data_quality()")
print(f"  This will verify the enhanced dataset")

print("\n🎉 Comprehensive cleanup completed successfully!")
