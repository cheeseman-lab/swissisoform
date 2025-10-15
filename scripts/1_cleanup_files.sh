#!/bin/bash
#
# SwissIsoform Pipeline Step 1: Cleanup Files
#
# This script cleans up ribosome profiling data and generates gene lists.
# It processes BED files and validates them against reference GTF files.
#
# Usage:
#   bash 1_cleanup_files.sh
#
# Prerequisites:
#   - dataset_config.yaml must exist (created by 0_download_genome.sh)
#   - GTF files referenced in config must be present
#   - BED files must be in ../data/ribosome_profiling/
#

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Start timing
START_TIME=$(date +%s)

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   SwissIsoform Pipeline Step 1: Cleanup Files                ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "This script processes ribosome profiling data and generates gene lists."
echo ""

# ============================================================================
# Validation
# ============================================================================

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Validation                                                  ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

GENOME_DIR="../data/genome_data"
RIBOPROF_DIR="../data/ribosome_profiling"
CONFIG_FILE="$RIBOPROF_DIR/dataset_config.yaml"

# Check for config file
echo -e "${YELLOW}→${NC} Checking for configuration file..."
if [ ! -f "$CONFIG_FILE" ]; then
    echo -e "${RED}✗${NC} $CONFIG_FILE missing"
    echo ""
    echo -e "${RED}Configuration file not found!${NC}"
    echo "Please create dataset_config.yaml or run 0_download_genome.sh"
    exit 1
fi
echo -e "${GREEN}✓${NC} dataset_config.yaml found"

# Check for GTF files referenced in config
echo ""
echo -e "${YELLOW}→${NC} Checking genome files..."
if [ -f "$GENOME_DIR/gencode.v25.annotation.gtf" ]; then
    SIZE=$(du -h "$GENOME_DIR/gencode.v25.annotation.gtf" | cut -f1)
    echo -e "${GREEN}✓${NC} gencode.v25.annotation.gtf (${SIZE})"
else
    echo -e "${YELLOW}⚠${NC} gencode.v25.annotation.gtf missing"
fi

if [ -f "$GENOME_DIR/gencode.v24.annotation.gtf" ]; then
    SIZE=$(du -h "$GENOME_DIR/gencode.v24.annotation.gtf" | cut -f1)
    echo -e "${GREEN}✓${NC} gencode.v24.annotation.gtf (${SIZE})"
else
    echo -e "${YELLOW}⚠${NC} gencode.v24.annotation.gtf missing"
fi

# Check for ribosome profiling BED files
echo ""
echo -e "${YELLOW}→${NC} Checking ribosome profiling files..."
mkdir -p "$RIBOPROF_DIR"  # Create directory if it doesn't exist

missing_riboprof=false

# Check for the new data files
for file in "HELA_Ly2024.bed" "HFF_Chen2020.bed" "IPSC_Chen2020.bed"; do
    if [ -f "$RIBOPROF_DIR/$file" ]; then
        SIZE=$(du -h "$RIBOPROF_DIR/$file" | cut -f1)
        echo -e "${GREEN}✓${NC} $file (${SIZE})"
    else
        echo -e "${RED}✗${NC} $file missing"
        missing_riboprof=true
    fi
done

if [ "$missing_riboprof" = true ]; then
    echo ""
    echo -e "${RED}Missing ribosome profiling BED files!${NC}"
    echo "Please place your experimental data files in $RIBOPROF_DIR"
    exit 1
fi

# ============================================================================
# Environment Setup
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Environment Setup                                           ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${YELLOW}→${NC} Activating conda environment..."
source ~/.bashrc
conda activate swissisoform || {
    echo -e "${RED}✗${NC} Failed to activate swissisoform conda environment"
    echo ""
    echo "Please run: conda env create --file=../environment.yml"
    exit 1
}
echo -e "${GREEN}✓${NC} Environment activated"

# ============================================================================
# Cleanup Processing
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Processing Datasets                                         ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${YELLOW}→${NC} Processing ribosome profiling datasets..."
echo ""
python3 cleanup_files.py

# ============================================================================
# Verification
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Verification                                                ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Parse datasets from config (simplified check)
datasets=("hela" "hff" "ipsc")

all_outputs_exist=true

echo -e "${YELLOW}→${NC} Checking output files..."
echo ""

for dataset in "${datasets[@]}"; do
    output_bed="$RIBOPROF_DIR/${dataset}_isoforms_with_transcripts.bed"
    gene_list="$RIBOPROF_DIR/${dataset}_isoforms_gene_list.txt"

    if [ -f "$output_bed" ]; then
        SIZE=$(du -h "$output_bed" | cut -f1)
        echo -e "${GREEN}✓${NC} ${dataset}_isoforms_with_transcripts.bed (${SIZE})"
    else
        echo -e "${RED}✗${NC} ${dataset}_isoforms_with_transcripts.bed missing"
        all_outputs_exist=false
    fi

    if [ -f "$gene_list" ]; then
        GENE_COUNT=$(wc -l < "$gene_list")
        echo -e "${GREEN}✓${NC} ${dataset}_isoforms_gene_list.txt ($GENE_COUNT genes)"
    else
        echo -e "${RED}✗${NC} ${dataset}_isoforms_gene_list.txt missing"
        all_outputs_exist=false
    fi
done

if [ "$all_outputs_exist" = false ]; then
    echo ""
    echo -e "${RED}Some output files are missing!${NC}"
    exit 1
fi

# Calculate duration
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
if [ $DURATION -lt 60 ]; then
    DURATION_STR="${DURATION}s"
elif [ $DURATION -lt 3600 ]; then
    MINUTES=$((DURATION / 60))
    SECS=$((DURATION % 60))
    DURATION_STR="${MINUTES}m ${SECS}s"
else
    HOURS=$((DURATION / 3600))
    MINUTES=$(((DURATION % 3600) / 60))
    SECS=$((DURATION % 60))
    DURATION_STR="${HOURS}h ${MINUTES}m ${SECS}s"
fi

# ============================================================================
# Summary
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Cleanup Complete                                            ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${GREEN}✓ File cleanup completed successfully in $DURATION_STR${NC}"
echo ""
echo "Generated files per dataset:"
echo "  ├─ {dataset}_isoforms_with_transcripts.bed"
echo "  └─ {dataset}_isoforms_gene_list.txt"
echo ""
echo -e "${BLUE}Next step:${NC}"
echo "  Run: sbatch 2_analyze_mutations.sh"
echo ""
