#!/usr/bin/env python3

"""Cleanup files for the ribosome profiling data.

This script cleans up the GTF file and BED files for the ribosome profiling data.
It updates the gene names in the GTF file to match the latest version of the GTF file.
It also cleans up the BED files to update gene names and remove duplicates.
"""

from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.utils import cleanup_bed, update_gencode_gene_names

# GTF: update gene names 
input_gtf = "../data/genome_data/gencode.v25.annotation.gtf"
output_gtf = "../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf"
reference_gtf = "../data/genome_data/gencode.v47.annotation.gtf"

update_gencode_gene_names(
    input_gtf_path=input_gtf,
    output_gtf_path=output_gtf,
    reference_gtf_path=reference_gtf,
    verbose=True,
)

# All truncations: bed cleanup
input_bed = "../data/ribosome_profiling/full_truncations_JL.bed"
output_bed = "../data/ribosome_profiling/full_truncations_JL_cleaned.bed"

cleanup_bed(input_bed, output_bed, gtf_path=output_gtf, verbose=True)

alt_isoforms = AlternativeIsoform()
alt_isoforms.load_bed("../data/ribosome_profiling/full_truncations_JL_cleaned.bed")
gene_list = alt_isoforms.get_gene_list()

with open("../data/ribosome_profiling/gene_list.txt", "w") as f:
    for gene in gene_list:
        f.write(gene + "\n")

# Selected truncations: bed cleanup
input_bed = "../data/ribosome_profiling/selected_truncations_JL.bed"
output_bed = "../data/ribosome_profiling/selected_truncations_JL_cleaned.bed"

cleanup_bed(input_bed, output_bed, gtf_path=output_gtf, verbose=True)

alt_isoforms = AlternativeIsoform()
alt_isoforms.load_bed("../data/ribosome_profiling/selected_truncations_JL_cleaned.bed")
gene_list = alt_isoforms.get_gene_list()

with open("../data/ribosome_profiling/gene_list_reduced.txt", "w") as f:
    for gene in gene_list:
        f.write(gene + "\n")