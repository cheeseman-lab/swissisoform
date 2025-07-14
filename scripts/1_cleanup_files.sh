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
    "$RIBOPROF_DIR/full_truncations_JL.bed"
    "$RIBOPROF_DIR/selected_truncations_JL.bed"
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
    echo "  - full_truncations_JL.bed (all detected truncation sites)"
    echo "  - selected_truncations_JL.bed (curated truncation sites)"
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
echo ""
echo "Verifying cleanup outputs..."

expected_outputs=(
    "$GENOME_DIR/gencode.v25.annotation.ensembl_cleaned.gtf"
    "$RIBOPROF_DIR/full_truncations_JL_cleaned.bed"
    "$RIBOPROF_DIR/selected_truncations_JL_cleaned.bed"
    "$RIBOPROF_DIR/gene_list.txt"
    "$RIBOPROF_DIR/gene_list_reduced.txt"
)

all_outputs_present=true
for file in "${expected_outputs[@]}"; do
    if [ -f "$file" ]; then
        if [[ "$file" == *.txt ]]; then
            count=$(wc -l < "$file")
            echo "✓ $(basename $file) ($count genes)"
        else
            echo "✓ $(basename $file)"
        fi
    else
        echo "✗ $(basename $file) missing"
        all_outputs_present=false
    fi
done

if [ "$all_outputs_present" = true ]; then
    echo ""
    echo "🎉 File cleanup completed successfully!"
    echo ""
    echo "Generated files:"
    echo "  ├─ Cleaned GTF annotation"
    echo "  ├─ Cleaned BED files (full and selected truncations)"
    echo "  └─ Gene lists for analysis"
    echo ""
    echo "Next step:"
    echo "  Run: sbatch 3_generate_proteins.sh"
else
    echo ""
    echo "❌ Cleanup failed. Some output files are missing."
    echo "Check the cleanup_files.py script for errors."
    exit 1
fi