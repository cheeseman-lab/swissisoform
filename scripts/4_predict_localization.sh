#!/bin/bash
#
# SwissIsoform Pipeline Step 4: Predict Localization (Parallel)
#
# This script generates subcellular localization predictions using DeepLoc.
# It processes protein sequences in parallel with both Fast and Accurate modes.
#
# Usage:
#   sbatch 4_predict_localization.sh
#   sbatch --export=DATASET=hela 4_predict_localization.sh
#   sbatch --export=DATASET=hela_bch 4_predict_localization.sh
#
# Environment Variables:
#   DATASET - Dataset to process (default: hela)
#
# Prerequisites:
#   - 3_generate_proteins.sh must have been run for the dataset
#   - Protein sequence files must exist
#   - DeepLoc conda environment must be installed
#

#SBATCH --job-name=deeploc                 # Job name
#SBATCH --partition=nvidia-A4000-20        # GPU partition
#SBATCH --array=1-4                        # 4 tasks for parallel processing
#SBATCH --cpus-per-task=4                  # CPUs per task
#SBATCH --mem=36G                          # Memory per task
#SBATCH --gres=gpu:1                       # 1 GPU per task
#SBATCH --time=36:00:00                    # Time limit
#SBATCH --output=out/deeploc-%A_%a.out     # %A = job ID, %a = array task ID

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

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   SwissIsoform Pipeline Step 4: Predict Localization         ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
echo "Dataset: $DATASET"
echo ""

# ============================================================================
# GPU Check
# ============================================================================

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  GPU Availability                                            ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${YELLOW}→${NC} Checking GPU availability..."
nvidia-smi

# Set environment variables
export HF_HOME="../.cache/huggingface"
mkdir -p "$HF_HOME"

# ============================================================================
# Validation (Task 1 Only)
# ============================================================================

if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Validation                                                  ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    echo -e "${YELLOW}→${NC} Checking for input files..."
    echo ""

    missing_files=()

    pairs_file="../results/$DATASET/proteins/protein_sequences_pairs.fasta"
    mutations_file="../results/$DATASET/proteins/protein_sequences_with_mutations.fasta"

    if [ -f "$pairs_file" ]; then
        count=$(grep -c '^>' "$pairs_file")
        echo -e "${GREEN}✓${NC} $DATASET/protein_sequences_pairs.fasta ($count sequences)"
    else
        echo -e "${RED}✗${NC} $DATASET/protein_sequences_pairs.fasta missing"
        missing_files+=("$DATASET:pairs")
    fi

    if [ -f "$mutations_file" ]; then
        count=$(grep -c '^>' "$mutations_file")
        echo -e "${GREEN}✓${NC} $DATASET/protein_sequences_with_mutations.fasta ($count sequences)"
    else
        echo -e "${RED}✗${NC} $DATASET/protein_sequences_with_mutations.fasta missing"
        missing_files+=("$DATASET:mutations")
    fi

    if [ ${#missing_files[@]} -gt 0 ]; then
        echo ""
        echo -e "${RED}Missing required files!${NC}"
        echo "Run 3_generate_proteins.sh first for dataset: $DATASET"
        exit 1
    fi

    # Create localization output directories
    echo ""
    echo -e "${YELLOW}→${NC} Creating output directories..."
    mkdir -p ../results/$DATASET/localization
    echo -e "${GREEN}✓${NC} Directories created"
fi

# ============================================================================
# Environment Setup
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Environment Setup                                           ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${YELLOW}→${NC} Activating DeepLoc environment..."
source ~/.bashrc
conda activate deeploc || {
    echo -e "${RED}✗${NC} Failed to activate deeploc environment"
    echo ""
    echo "To install DeepLoc:"
    echo "  1. Download DeepLoc 2.1 from: https://services.healthtech.dtu.dk/services/DeepLoc-2.0/"
    echo "  2. Place DeepLoc-2.1.0.tar.gz in the main swissisoform/ directory"
    echo "  3. conda create -n deeploc python=3.8"
    echo "  4. conda activate deeploc"
    echo "  5. pip install DeepLoc-2.1.0.tar.gz"
    exit 1
}
echo -e "${GREEN}✓${NC} Environment activated"

# ============================================================================
# Task Assignment
# ============================================================================

# Define which file and mode each task processes (4 tasks total)
case $SLURM_ARRAY_TASK_ID in
    1)
        file_type="pairs"; mode="Fast"
        input_file="../results/$DATASET/proteins/protein_sequences_pairs.fasta"
        ;;
    2)
        file_type="pairs"; mode="Accurate"
        input_file="../results/$DATASET/proteins/protein_sequences_pairs.fasta"
        ;;
    3)
        file_type="mutations"; mode="Fast"
        input_file="../results/$DATASET/proteins/protein_sequences_with_mutations.fasta"
        ;;
    4)
        file_type="mutations"; mode="Accurate"
        input_file="../results/$DATASET/proteins/protein_sequences_with_mutations.fasta"
        ;;
    *)
        echo -e "${RED}✗${NC} Unknown array task ID: $SLURM_ARRAY_TASK_ID"
        exit 1
        ;;
esac

# ============================================================================
# DeepLoc Prediction
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Processing Task ${SLURM_ARRAY_TASK_ID}                                            ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${YELLOW}→${NC} Processing $DATASET dataset ($file_type sequences, $mode mode)"
echo "  Input file: $(basename $input_file)"
echo ""

if [ -f "$input_file" ]; then
    # Remove outputs directory if it exists
    if [ -d "outputs/" ]; then
        echo -e "${YELLOW}→${NC} Removing existing outputs/ directory..."
        rm -rf outputs/
    fi

    # Create descriptive temporary subfolder for this specific run
    temp_subdir="../results/$DATASET/localization/${DATASET}_${file_type}_${mode}_temp_$$"
    mkdir -p "$temp_subdir"

    echo -e "${YELLOW}→${NC} Starting DeepLoc $mode mode for $DATASET ($file_type) at $(date)"
    echo ""

    # Set GPU memory growth to avoid OOM errors
    export TF_FORCE_GPU_ALLOW_GROWTH=true
    export CUDA_VISIBLE_DEVICES=0

    deeploc2 -f "$input_file" -m "$mode" -o "$temp_subdir/" -d cuda

    # Find and move the results file
    result_file=$(find "$temp_subdir" -name "results_*.csv" | head -n 1)
    if [ -n "$result_file" ] && [ -f "$result_file" ]; then
        output_file="../results/$DATASET/localization/protein_sequences_${file_type}_${mode}_results.csv"
        mv "$result_file" "$output_file"
        echo ""
        echo -e "${GREEN}✓${NC} Moved $mode results to protein_sequences_${file_type}_${mode}_results.csv"
    else
        echo ""
        echo -e "${RED}✗${NC} $mode results not found"
        exit 1
    fi

    # Clean up temp subfolder
    rm -rf "$temp_subdir"

    # Remove outputs directory if it exists
    if [ -d "outputs/" ]; then
        echo -e "${YELLOW}→${NC} Removing outputs/ directory..."
        rm -rf outputs/
    fi

    echo ""
    echo -e "${GREEN}✓${NC} Completed $DATASET ($file_type, $mode) at $(date)"
else
    echo -e "${YELLOW}⚠${NC} Skipping $DATASET ($file_type, $mode) - input file not found"
    exit 1
fi

# ============================================================================
# Verification (Task 8 Only)
# ============================================================================

if [ "$SLURM_ARRAY_TASK_ID" -eq 4 ]; then
    # Wait for file system sync
    sleep 10

    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Verification                                                ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    echo -e "${YELLOW}→${NC} Verifying DeepLoc outputs for $DATASET..."
    echo ""

    # Check for outputs from all tasks
    file_types=("pairs" "mutations")
    modes=("Fast" "Accurate")
    found_outputs=()

    for file_type in "${file_types[@]}"; do
        for mode in "${modes[@]}"; do
            output_file="../results/$DATASET/localization/protein_sequences_${file_type}_${mode}_results.csv"

            if [ -f "$output_file" ]; then
                count=$(($(wc -l < "$output_file") - 1))
                echo -e "${GREEN}✓${NC} $DATASET/$(basename $output_file) ($count predictions)"
                found_outputs+=("$output_file")
            else
                echo -e "${YELLOW}⚠${NC} $DATASET/$(basename $output_file) missing"
            fi
        done
    done

    # Summary
    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Localization Prediction Complete                           ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    if [ ${#found_outputs[@]} -eq 0 ]; then
        echo -e "${RED}✗ No DeepLoc outputs were generated!${NC}"
        echo "Check the SLURM logs for DeepLoc error messages."

        # Debug: List what files exist in localization directories
        echo ""
        echo "Debug: Files in localization directory:"
        localization_dir="../results/$DATASET/localization"
        if [ -d "$localization_dir" ]; then
            echo "  $DATASET/localization/:"
            ls -la "$localization_dir" | sed 's/^/    /'
        fi

        exit 1
    else
        echo -e "${GREEN}✓ DeepLoc predictions completed successfully for $DATASET!${NC}"
        echo ""
        echo "Generated predictions:"
        echo "  └─ $DATASET/localization/"
        for output in "${found_outputs[@]}"; do
            echo "     ├─ $(basename $output)"
        done
        echo ""
        echo "Summary: Generated ${#found_outputs[@]} prediction files"
        echo ""
        echo -e "${BLUE}Next step:${NC}"
        echo "  Run: sbatch --export=DATASET=${DATASET} scripts/5_summarize_results.sh"
        echo ""
    fi
fi
