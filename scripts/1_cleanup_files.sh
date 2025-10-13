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
CONFIG_FILE="$RIBOPROF_DIR/dataset_config.yaml"

# Check for config file
if [ ! -f "$CONFIG_FILE" ]; then
    echo "✗ $CONFIG_FILE missing"
    echo "❌ Configuration file not found!"
    echo "Please create dataset_config.yaml or run 0_download_genome.sh"
    exit 1
fi

echo "✓ dataset_config.yaml found"

# Check for GTF files referenced in config
echo ""
echo "Checking genome files..."

if [ -f "$GENOME_DIR/gencode.v25.annotation.gtf" ]; then
    echo "✓ gencode.v25.annotation.gtf"
else
    echo "✗ gencode.v25.annotation.gtf missing"
fi

if [ -f "$GENOME_DIR/gencode.v24.annotation.gtf" ]; then
    echo "✓ gencode.v24.annotation.gtf"
else
    echo "✗ gencode.v24.annotation.gtf missing"
fi

# Check for ribosome profiling BED files
echo ""
echo "Checking ribosome profiling files..."
mkdir -p "$RIBOPROF_DIR"  # Create directory if it doesn't exist

missing_riboprof=false

# Check for the new data files
for file in "HELA_Ly2024.bed" "HFF_Chen2020.bed" "IPSC_Chen2020.bed"; do
    if [ -f "$RIBOPROF_DIR/$file" ]; then
        echo "✓ $file"
    else
        echo "✗ $file missing"
        missing_riboprof=true
    fi
done

if [ "$missing_riboprof" = true ]; then
    echo ""
    echo "❌ Missing ribosome profiling BED files!"
    echo "Please place your experimental data files in $RIBOPROF_DIR"
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

# Verify outputs for each dataset
echo ""
echo "Verifying cleanup outputs..."

# Parse datasets from config (simplified check)
datasets=("hela" "hff" "ipsc")

all_outputs_exist=true

for dataset in "${datasets[@]}"; do
    output_bed="$RIBOPROF_DIR/${dataset}_isoforms_with_transcripts.bed"
    gene_list="$RIBOPROF_DIR/${dataset}_isoforms_gene_list.txt"

    if [ -f "$output_bed" ]; then
        echo "✓ ${dataset}_isoforms_with_transcripts.bed"
    else
        echo "✗ ${dataset}_isoforms_with_transcripts.bed"
        all_outputs_exist=false
    fi

    if [ -f "$gene_list" ]; then
        GENE_COUNT=$(wc -l < "$gene_list")
        echo "✓ ${dataset}_isoforms_gene_list.txt ($GENE_COUNT genes)"
    else
        echo "✗ ${dataset}_isoforms_gene_list.txt"
        all_outputs_exist=false
    fi
done

if [ "$all_outputs_exist" = false ]; then
    echo ""
    echo "❌ Some output files are missing!"
    exit 1
fi

echo ""
echo "🎉 File cleanup completed successfully!"
echo "Generated files per dataset:"
echo "  ├─ {dataset}_isoforms_with_transcripts.bed"
echo "  └─ {dataset}_isoforms_gene_list.txt"

echo ""
echo "Next step:"
echo "  Run: sbatch 2_analyze_mutations.sh"
