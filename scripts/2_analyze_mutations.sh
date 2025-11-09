#!/bin/bash
#
# SwissIsoform Pipeline Step 2: Analyze Mutations (Parallel)
#
# This script analyzes mutations in alternative isoform truncation regions with
# gene-level parallelization. It reads dataset configuration from YAML file and
# processes chunks in parallel using SLURM array jobs.
#
# Usage:
#   sbatch 2_analyze_mutations.sh
#   sbatch --export=DATASET=hela 2_analyze_mutations.sh
#   sbatch --export=DATASET=hela,SOURCES="clinvar|cosmic" 2_analyze_mutations.sh
#   sbatch --export=DATASET=hela,SOURCES="clinvar",IMPACT_TYPES="missense variant|nonsense variant" 2_analyze_mutations.sh
#   sbatch --export=DATASET=hela,CUSTOM_PARQUET="/path/to/custom_mutations.parquet" 2_analyze_mutations.sh
#   sbatch --export=DATASET=hela,OUTPUT_NAME="hela_bch",CUSTOM_PARQUET="/path/to/bch.parquet" 2_analyze_mutations.sh
#
# Environment Variables:
#   DATASET - Dataset to process (default: hela)
#   OUTPUT_NAME - Custom output folder name (default: same as DATASET)
#   SOURCES - Mutation databases to query (default: clinvar, pipe-separated: clinvar|gnomad|cosmic)
#   IMPACT_TYPES - Mutation types to analyze (default: all types, pipe-separated)
#   CUSTOM_PARQUET - Path to custom parquet file with mutation data (optional)
#
# Prerequisites:
#   - 1_cleanup_files.sh must have been run
#   - Gene lists must exist for the dataset
#   - Reference files must be present
#

#SBATCH --job-name=mutations               # Job name
#SBATCH --partition=20                     # Partition name
#SBATCH --array=1-8                        # 8 chunks for parallel processing
#SBATCH --cpus-per-task=4                  # CPUs per task
#SBATCH --mem=64G                          # Memory per task
#SBATCH --time=24:00:00                    # Time limit (hrs:min:sec)
#SBATCH --output=out/mutations-%A_%a.out   # %A = job ID, %a = array task ID

set -e  # Exit on error

# Source shared utilities (colors, helper functions)
# Use relative path from scripts directory
UTILS_PATH="$(dirname "$0")/utils.sh"
if [ -f "$UTILS_PATH" ]; then
    source "$UTILS_PATH"
elif [ -f "utils.sh" ]; then
    source "utils.sh"
elif [ -f "scripts/utils.sh" ]; then
    source "scripts/utils.sh"
else
    # Fallback: define colors inline
    if [ -t 1 ]; then
        RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
        BLUE='\033[0;34m'; CYAN='\033[0;36m'; NC='\033[0m'
    else
        RED=''; GREEN=''; YELLOW=''; BLUE=''; CYAN=''; NC=''
    fi
fi

# Dataset selection (default: hela)
DATASET="${DATASET:-hela}"

# Output name (default: same as DATASET, but can be overridden for custom output folders)
OUTPUT_NAME="${OUTPUT_NAME:-$DATASET}"

# Sources selection (default: clinvar)
# Convert pipe-separated string to space-separated for command line
if [ -z "$SOURCES" ]; then
    SOURCES_ARGS=("clinvar")
else
    IFS='|' read -ra SOURCES_ARGS <<< "$SOURCES"
fi

# Impact types selection (default: all types)
# Convert pipe-separated string to space-separated for command line
if [ -z "$IMPACT_TYPES" ]; then
    IMPACT_TYPES_ARGS=("missense variant" "nonsense variant" "frameshift variant" "synonymous variant" "inframe deletion" "inframe insertion")
else
    IFS='|' read -ra IMPACT_TYPES_ARGS <<< "$IMPACT_TYPES"
fi

# Custom parquet file (optional)
CUSTOM_PARQUET="${CUSTOM_PARQUET:-}"

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   SwissIsoform Pipeline Step 2: Analyze Mutations            ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
echo "Dataset: $DATASET"
echo "Output name: $OUTPUT_NAME"
echo "Sources: ${SOURCES_ARGS[@]}"
echo "Impact types: ${IMPACT_TYPES_ARGS[@]}"
if [ -n "$CUSTOM_PARQUET" ]; then
    echo "Custom parquet: $CUSTOM_PARQUET"
fi
echo ""

# ============================================================================
# Helper Functions
# ============================================================================

# Function to split gene list into chunks
split_gene_list() {
    local input_file=$1
    local total_chunks=$2
    local chunk_id=$3
    local output_file=$4

    total_genes=$(wc -l < "$input_file")
    genes_per_chunk=$(( (total_genes + total_chunks - 1) / total_chunks ))
    start_line=$(( (chunk_id - 1) * genes_per_chunk + 1 ))
    end_line=$(( chunk_id * genes_per_chunk ))

    sed -n "${start_line},${end_line}p" "$input_file" > "$output_file"
}

# ============================================================================
# Configuration and Validation (Task 1 Only)
# ============================================================================

if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Configuration                                               ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    echo -e "${YELLOW}→${NC} Reading dataset configuration..."
    CONFIG_FILE="../data/ribosome_profiling/dataset_config.yaml"

    if [ ! -f "$CONFIG_FILE" ]; then
        echo -e "${RED}✗${NC} Dataset configuration file not found: $CONFIG_FILE"
        exit 1
    fi

    # Use Python to parse YAML and get dataset-specific paths
    read -r DATASET_BED DATASET_GTF GENE_LIST <<< $(python3 -c "
import yaml
import sys
from pathlib import Path

config_file = '$CONFIG_FILE'
dataset = '$DATASET'

with open(config_file) as f:
    config = yaml.safe_load(f)

# Find the dataset
dataset_config = None
for ds in config['datasets']:
    if ds['name'] == dataset:
        dataset_config = ds
        break

if not dataset_config:
    print(f'ERROR: Dataset {dataset} not found in config', file=sys.stderr)
    sys.exit(1)

# Get paths
bed_file = f\"../data/ribosome_profiling/{dataset}_isoforms_with_transcripts.bed\"
gene_list = f\"../data/ribosome_profiling/{dataset}_isoforms_gene_list.txt\"

# Get GTF path - prefer v47names version if it exists
source_gtf_path = Path(dataset_config['source_gtf_path'])
v47_gtf = source_gtf_path.parent / f\"{source_gtf_path.stem}.v47names.gtf\"

if v47_gtf.exists():
    gtf_file = str(v47_gtf)
else:
    gtf_file = dataset_config['source_gtf_path']

print(bed_file, gtf_file, gene_list)
")

    if [ $? -ne 0 ]; then
        echo -e "${RED}✗${NC} Failed to read dataset configuration"
        exit 1
    fi

    echo -e "${GREEN}✓${NC} Configuration loaded"
    echo ""
    echo "  BED file: $(basename $DATASET_BED)"
    echo "  GTF file: $(basename $DATASET_GTF)"
    echo "  Gene list: $(basename $GENE_LIST)"

    # Check if required input files exist
    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Validation                                                  ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    GENOME_DIR="../data/genome_data"
    GENOME_PATH="$GENOME_DIR/GRCh38.p7.genome.fa"

    required_files=(
        "$GENOME_PATH"
        "$DATASET_GTF"
        "$DATASET_BED"
        "$GENE_LIST"
    )

    echo -e "${YELLOW}→${NC} Checking required files..."
    echo ""
    missing_files=false
    for file in "${required_files[@]}"; do
        if [ -f "$file" ]; then
            SIZE=$(du -h "$file" | cut -f1)
            echo -e "${GREEN}✓${NC} $(basename $file) (${SIZE})"
        else
            echo -e "${RED}✗${NC} $(basename $file) missing"
            missing_files=true
        fi
    done

    if [ "$missing_files" = true ]; then
        echo ""
        echo -e "${RED}Missing required files! Run 1_cleanup_files.sh first${NC}"
        exit 1
    fi

    # Create results directory structure
    echo ""
    echo -e "${YELLOW}→${NC} Creating output directories..."
    mkdir -p "../results/${OUTPUT_NAME}/mutations"
    mkdir -p ../results/temp/chunks
    echo -e "${GREEN}✓${NC} Directories created"
else
    # Other tasks need to read the config too
    CONFIG_FILE="../data/ribosome_profiling/dataset_config.yaml"
    read -r DATASET_BED DATASET_GTF GENE_LIST <<< $(python3 -c "
import yaml
from pathlib import Path

with open('$CONFIG_FILE') as f:
    config = yaml.safe_load(f)

dataset_config = next(ds for ds in config['datasets'] if ds['name'] == '$DATASET')
bed_file = f'../data/ribosome_profiling/${DATASET}_isoforms_with_transcripts.bed'
gene_list = f'../data/ribosome_profiling/${DATASET}_isoforms_gene_list.txt'

source_gtf_path = Path(dataset_config['source_gtf_path'])
v47_gtf = source_gtf_path.parent / f\"{source_gtf_path.stem}.v47names.gtf\"
gtf_file = str(v47_gtf) if v47_gtf.exists() else dataset_config['source_gtf_path']

print(bed_file, gtf_file, gene_list)
")
fi

# Wait for setup to complete and add staggered delays to avoid API rate limits
sleep 2
STAGGER_DELAY=$(( (SLURM_ARRAY_TASK_ID - 1) * 3 ))
echo -e "${YELLOW}→${NC} Adding ${STAGGER_DELAY}s staggered delay for task ${SLURM_ARRAY_TASK_ID}"
sleep $STAGGER_DELAY

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
# Mutation Analysis
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Processing Chunk ${SLURM_ARRAY_TASK_ID}                                          ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Define common paths
GENOME_PATH="../data/genome_data/GRCh38.p7.genome.fa"
ANNOTATION_PATH="$DATASET_GTF"
TRUNCATIONS_PATH="$DATASET_BED"
OUTPUT_DIR="../results/${OUTPUT_NAME}/mutations/chunk_${SLURM_ARRAY_TASK_ID}"
CHUNK_ID=$SLURM_ARRAY_TASK_ID
TOTAL_CHUNKS=8

# Create task-specific output directory
mkdir -p "$OUTPUT_DIR"

# Create gene chunk for this task
CHUNK_FILE="../results/temp/chunks/${DATASET}_chunk_${CHUNK_ID}.txt"
split_gene_list "$GENE_LIST" "$TOTAL_CHUNKS" "$CHUNK_ID" "$CHUNK_FILE"

# Skip if no genes in chunk
if [ ! -s "$CHUNK_FILE" ]; then
    echo -e "${YELLOW}⚠${NC} No genes in chunk ${CHUNK_ID}, exiting"
    exit 0
fi

GENE_COUNT=$(wc -l < "$CHUNK_FILE")
echo -e "${YELLOW}→${NC} Processing ${DATASET} dataset chunk ${CHUNK_ID}"
echo "  Gene list: $GENE_COUNT genes"
echo ""

# Run analysis on the gene chunk
if [ -n "$CUSTOM_PARQUET" ]; then
    python3 analyze_mutations.py "$CHUNK_FILE" "$OUTPUT_DIR" \
      --genome "$GENOME_PATH" \
      --annotation "$ANNOTATION_PATH" \
      --bed "$TRUNCATIONS_PATH" \
      --sources "${SOURCES_ARGS[@]}" \
      --custom-parquet "$CUSTOM_PARQUET" \
      --impact-types "${IMPACT_TYPES_ARGS[@]}" \
      --visualize \
      -v
else
    python3 analyze_mutations.py "$CHUNK_FILE" "$OUTPUT_DIR" \
      --genome "$GENOME_PATH" \
      --annotation "$ANNOTATION_PATH" \
      --bed "$TRUNCATIONS_PATH" \
      --sources "${SOURCES_ARGS[@]}" \
      --impact-types "${IMPACT_TYPES_ARGS[@]}" \
      --visualize \
      -v
fi

echo ""
echo -e "${GREEN}✓${NC} Completed ${DATASET} dataset chunk ${CHUNK_ID}"

# Clean up chunk file
rm -f "$CHUNK_FILE"

# ============================================================================
# Result Merging and Verification (Task 8 Only)
# ============================================================================

if [ "$SLURM_ARRAY_TASK_ID" -eq 8 ]; then
    # Wait for other tasks to finish
    echo ""
    echo -e "${YELLOW}→${NC} Waiting for all tasks to complete..."
    sleep 30

    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Merging Results                                             ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    # Create final output directory
    FINAL_OUTPUT_DIR="../results/${OUTPUT_NAME}/mutations"
    mkdir -p "$FINAL_OUTPUT_DIR"

    # Merge gene level results
    echo -e "${YELLOW}→${NC} Merging gene level results..."
    first_file=true
    for i in {1..8}; do
        chunk_file="../results/${OUTPUT_NAME}/mutations/chunk_${i}/gene_level_results.csv"
        if [ -f "$chunk_file" ]; then
            if [ "$first_file" = true ]; then
                cp "$chunk_file" "$FINAL_OUTPUT_DIR/gene_level_results.csv"
                first_file=false
                rows=$(wc -l < "$chunk_file")
                echo "  ├─ Added chunk $i with header ($rows rows)"
            else
                tail -n +2 "$chunk_file" >> "$FINAL_OUTPUT_DIR/gene_level_results.csv"
                rows=$(($(wc -l < "$chunk_file") - 1))
                echo "  ├─ Added chunk $i ($rows rows)"
            fi
        else
            echo -e "  ├─ ${YELLOW}⚠${NC} Warning: chunk $i gene results missing"
        fi
    done

    # Merge isoform level results
    echo ""
    echo -e "${YELLOW}→${NC} Merging isoform level results..."
    first_file=true
    for i in {1..8}; do
        chunk_file="../results/${OUTPUT_NAME}/mutations/chunk_${i}/isoform_level_results.csv"
        if [ -f "$chunk_file" ]; then
            if [ "$first_file" = true ]; then
                cp "$chunk_file" "$FINAL_OUTPUT_DIR/isoform_level_results.csv"
                first_file=false
                rows=$(wc -l < "$chunk_file")
                echo "  ├─ Added chunk $i with header ($rows rows)"
            else
                tail -n +2 "$chunk_file" >> "$FINAL_OUTPUT_DIR/isoform_level_results.csv"
                rows=$(($(wc -l < "$chunk_file") - 1))
                echo "  ├─ Added chunk $i ($rows rows)"
            fi
        else
            echo -e "  ├─ ${YELLOW}⚠${NC} Warning: chunk $i isoform results missing"
        fi
    done

    # Merge visualization files
    echo ""
    echo -e "${YELLOW}→${NC} Merging visualization files..."
    for i in {1..8}; do
        chunk_viz_dir="../results/${OUTPUT_NAME}/mutations/chunk_${i}"
        if [ -d "$chunk_viz_dir" ]; then
            find "$chunk_viz_dir" -maxdepth 1 -type d -not -name "chunk_*" -exec mv {} "$FINAL_OUTPUT_DIR/" \; 2>/dev/null || true
        fi
    done

    # Verification
    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Verification                                                ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    expected_files=(
        "../results/${OUTPUT_NAME}/mutations/gene_level_results.csv"
        "../results/${OUTPUT_NAME}/mutations/isoform_level_results.csv"
    )

    all_files_present=true
    for file in "${expected_files[@]}"; do
        if [ -f "$file" ]; then
            count=$(($(wc -l < "$file") - 1))
            echo -e "${GREEN}✓${NC} $(basename $file) ($count rows)"
        else
            echo -e "${RED}✗${NC} $(basename $file) missing"
            all_files_present=false
        fi
    done

    # Check for visualization outputs
    echo ""
    vis_count=$(find "../results/${OUTPUT_NAME}/mutations" -name "*.pdf" -not -path "*/chunk_*/*" 2>/dev/null | wc -l)
    echo -e "${GREEN}✓${NC} ${OUTPUT_NAME} dataset visualizations: $vis_count PDFs"

    # Summary
    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Analysis Complete                                           ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    if [ "$all_files_present" = true ]; then
        echo -e "${GREEN}✓ ${OUTPUT_NAME} mutation analysis completed successfully!${NC}"
        echo ""
        echo "Generated analysis results:"
        echo "  └─ ${OUTPUT_NAME}/mutations/"
        echo "     ├─ gene_level_results.csv"
        echo "     ├─ isoform_level_results.csv"
        echo "     └─ [gene_name]/ (visualizations)"
        echo ""
        echo -e "${BLUE}Next step:${NC}"
        echo "  Run: sbatch --export=DATASET=${DATASET},OUTPUT_NAME=${OUTPUT_NAME} scripts/3_generate_proteins.sh"
        echo ""

        # Clean up chunk directories
        echo -e "${YELLOW}→${NC} Cleaning up chunk directories..."
        rm -rf "../results/${OUTPUT_NAME}/mutations/chunk_"*
        echo -e "${GREEN}✓${NC} Cleanup complete"
    else
        echo -e "${RED}✗ ${OUTPUT_NAME} mutation analysis failed${NC}"
        echo "Some output files are missing. Check the logs for errors."
        exit 1
    fi
fi

# Final cleanup
if [ "$SLURM_ARRAY_TASK_ID" -eq 8 ]; then
    sleep 30
    rm -rf ../results/temp/chunks
fi
