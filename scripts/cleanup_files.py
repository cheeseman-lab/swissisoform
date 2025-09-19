#!/usr/bin/env python3
"""Cleanup files for the ribosome profiling data with RefSeq transcript mapping."""

from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.utils import (
    update_gencode_gene_names,
    comprehensive_cleanup_bed_with_refseq_gtf_direct,
    subset_gene_list,
)

# File paths
input_gtf = "../data/genome_data/gencode.v25.annotation.gtf"
output_gtf = "../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf"
reference_gtf = "../data/genome_data/gencode.v47.annotation.gtf"
preferred_transcripts = "../data/genome_data/hela_top_transcript.txt"
input_bed = "../data/ribosome_profiling/Ly_2024b_TableS2_formatted.bed"
enhanced_bed = "../data/ribosome_profiling/isoforms_with_transcripts.bed"

print("=== SwissIsoform Data Cleanup Pipeline ===")

# Step 1: Update gene names
print("Step 1: Updating GTF gene names...")
update_gencode_gene_names(
    input_gtf_path=input_gtf,
    output_gtf_path=output_gtf,
    reference_gtf_path=reference_gtf,
    verbose=False,  # Reduce verbosity
)
print("  ✓ Gene names updated")

# Step 2: Enhanced BED cleanup with RefSeq mapping
print("Step 2: Processing BED file with RefSeq mapping...")
cleanup_stats = comprehensive_cleanup_bed_with_refseq_gtf_direct(
    input_bed_path=input_bed,
    gtf_path=output_gtf,
    preferred_transcripts_path=preferred_transcripts,
    output_bed_path=enhanced_bed,
    verbose=False,  # Reduce verbosity
)
print(f"  ✓ Enhanced {cleanup_stats['enhanced_entries']} entries")
print(f"  ✓ RefSeq mapping: {cleanup_stats['refseq_mapping_rate']:.1f}%")

# Step 3: Generate gene lists
print("Step 3: Generating gene lists...")
alt_isoforms = AlternativeIsoform()
alt_isoforms.load_bed(enhanced_bed)
gene_list = alt_isoforms.get_gene_list()

# Write gene lists
gene_list_path = "../data/ribosome_profiling/isoforms_gene_list.txt"
with open(gene_list_path, "w") as f:
    for gene in gene_list:
        f.write(gene + "\n")

reduced_gene_list = subset_gene_list(gene_list_path=gene_list_path)
reduced_gene_list_path = "../data/ribosome_profiling/isoforms_gene_list_reduced.txt"
with open(reduced_gene_list_path, "w") as f:
    for gene in reduced_gene_list:
        f.write(gene + "\n")

print(f"  ✓ Full gene list: {len(gene_list)} genes")
print(f"  ✓ Reduced gene list: {len(reduced_gene_list)} genes")

# Final summary
print("\n=== Summary ===")
print(f"Total entries processed: {cleanup_stats['enhanced_entries']}")
print(f"RefSeq mapping success: {cleanup_stats['refseq_mapping_rate']:.1f}%")
print(f"Final gene count: {cleanup_stats['final_genes']}")

# Assessment
refseq_rate = cleanup_stats.get("refseq_mapping_rate", 0)
if refseq_rate > 80:
    print("✓ Excellent RefSeq coverage")
elif refseq_rate > 60:
    print("✓ Good RefSeq coverage")
elif refseq_rate > 40:
    print("⚠ Moderate RefSeq coverage")
else:
    print("❌ Low RefSeq coverage")

print(f"\nOutput files:")
print(f"  ├─ Enhanced BED: {enhanced_bed}")
print(f"  ├─ Gene lists: {len(gene_list)} total, {len(reduced_gene_list)} reduced")
print(f"  └─ Ready for analysis")

print("\n✓ Cleanup completed successfully!")
