#!/bin/bash

# 1_cleanup_files.sh
# Cleans up ribosome profiling data and generates gene lists

set -e  # Exit on any error

echo "=============================================="
echo "SwissIsoform Pipeline Step 1: Cleanup Files"
echo "=============================================="

# Check if required input files exist
echo "Checking for required input files..."

GENOME_DIR="../data/genome_data"
RIBOPROF_DIR="../data/ribosome_profiling"

required_genome_files=(
    "$GENOME_DIR/gencode.v25.annotation.gtf"
    "$GENOME_DIR/gencode.v47.annotation.gtf"
)

required_riboprof_files=(
    "$RIBOPROF_DIR/Ly_2024b_TableS2_formatted.bed"
)

echo ""
echo "Checking genome files..."
for file in "${required_genome_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ $(basename $file)"
    else
        echo "✗ $file missing"
        echo "❌ Run 0_download_genome.sh first"
        exit 1
    fi
done

echo ""
echo "Checking ribosome profiling files..."
mkdir -p "$RIBOPROF_DIR"  # Create directory if it doesn't exist

missing_riboprof=false
for file in "${required_riboprof_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ $(basename $file)"
    else
        echo "✗ $(basename $file) missing"
        missing_riboprof=true
    fi
done

if [ "$missing_riboprof" = true ]; then
    echo ""
    echo "❌ Missing ribosome profiling BED files!"
    echo "Please place your experimental data files in $RIBOPROF_DIR:"
    echo "  - isoforms.bed (all detected isoforms sites)"
    exit 1
fi

# Activate conda environment
echo ""
echo "Activating conda environment..."
source ~/.bashrc
conda activate swissisoform || {
    echo "❌ Failed to activate swissisoform conda environment"
    echo "Please run: conda env create --file=../environment.yml"
    exit 1
}

# Run the cleanup Python script
echo ""
echo "Running cleanup script..."
python3 cleanup_files.py

# Verify outputs
echo "Verifying cleanup outputs..."

# Check GTF file
if [ -f "../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf" ]; then
    echo "✓ gencode.v25.annotation.ensembl_cleaned.gtf"
else
    echo "✗ gencode.v25.annotation.ensembl_cleaned.gtf"
    exit 1
fi

# Check filtered BED file
if [ -f "../data/ribosome_profiling/isoforms_with_transcripts.bed" ]; then
    echo "✓ isoforms_with_transcripts.bed"
else
    echo "✗ isoforms_with_transcripts.bed"
    exit 1
fi

# Check gene list file and validate gene count
if [ -f "../data/ribosome_profiling/isoforms_gene_list.txt" ]; then
    GENE_COUNT_TXT=$(wc -l < "../data/ribosome_profiling/isoforms_gene_list.txt")
    GENE_COUNT_BED=$(cut -f4 "../data/ribosome_profiling/isoforms_with_transcripts.bed" | cut -d'_' -f1 | sort | uniq | wc -l)
    
    if [ "$GENE_COUNT_TXT" -eq "$GENE_COUNT_BED" ]; then
        echo "✓ isoforms_gene_list.txt ($GENE_COUNT_TXT genes, matches BED file)"
    else
        echo "✗ isoforms_gene_list.txt ($GENE_COUNT_TXT genes, but BED has $GENE_COUNT_BED genes)"
        echo "  Gene count mismatch detected - check AlternativeIsoform.get_gene_list() method"
        exit 1
    fi
else
    echo "✗ isoforms_gene_list.txt"
    exit 1
fi

# Check reduced gene list
if [ -f "../data/ribosome_profiling/isoforms_gene_list_reduced.txt" ]; then
    REDUCED_COUNT=$(wc -l < "../data/ribosome_profiling/isoforms_gene_list_reduced.txt")
    echo "✓ isoforms_gene_list_reduced.txt ($REDUCED_COUNT genes)"
    
    if [ "$REDUCED_COUNT" -eq 0 ]; then
        echo "  WARNING: No genes found in reduced list - predefined genes may not exist in dataset"
        echo "  Consider updating the subset_gene_list() function with genes from your actual dataset"
    fi
else
    echo "✗ isoforms_gene_list_reduced.txt"
    exit 1
fi

echo "🎉 File cleanup completed successfully!"
echo "Generated files:"
echo "  ├─ Cleaned GTF annotation"
echo "  ├─ Cleaned and filtered BED files"
echo "  └─ Gene lists for analysis"

# Add diagnostic info
echo ""
echo "Summary:"
echo "  ├─ Total genes with alternatives: $GENE_COUNT_BED"
echo "  ├─ Genes in reduced list: $REDUCED_COUNT"
if [ "$REDUCED_COUNT" -gt 0 ]; then
    echo "  └─ Ready for analysis"
else
    echo "  └─ WARNING: No genes available for reduced analysis"
    echo ""
    echo "Suggested next steps:"
    echo "  1. Check which genes are available: head -20 ../data/ribosome_profiling/isoforms_gene_list.txt"
    echo "  2. Update subset_gene_list() with genes from your dataset"
fi

echo ""
echo "Next step:"
echo "  Run: sbatch 2_analyze_mutations.sh"
echo "  Run: sbatch 3_generate_proteins.sh"