#!/bin/bash
#
# SwissIsoform Pipeline Step 4: Predict Localization (Parallel)
#
# This script generates subcellular localization predictions using DeepLoc.
# It processes protein sequences in parallel across multiple datasets and modes.
#
# Usage:
#   sbatch scripts/4_predict_localization.sh
#
# Prerequisites:
#   - 3_generate_proteins.sh must have been run
#   - Protein sequence files must exist
#   - DeepLoc conda environment must be installed
#

#SBATCH --job-name=deeploc                 # Job name
#SBATCH --partition=nvidia-A4000-20        # GPU partition
#SBATCH --array=1-8                        # 8 tasks for parallel processing
#SBATCH --cpus-per-task=4                  # CPUs per task
#SBATCH --mem=36G                          # Memory per task
#SBATCH --gres=gpu:1                       # 1 GPU per task
#SBATCH --time=36:00:00                    # Time limit
#SBATCH --output=out/deeploc-%A_%a.out     # %A = job ID, %a = array task ID

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   SwissIsoform Pipeline Step 4: Predict Localization         ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
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

    datasets=("reduced" "full")
    available_files=()
    missing_files=()

    for dataset in "${datasets[@]}"; do
        pairs_file="../results/$dataset/proteins/protein_sequences_pairs.fasta"
        mutations_file="../results/$dataset/proteins/protein_sequences_with_mutations.fasta"

        if [ -f "$pairs_file" ]; then
            count=$(grep -c '^>' "$pairs_file")
            echo -e "${GREEN}✓${NC} $dataset/protein_sequences_pairs.fasta ($count sequences)"
        else
            echo -e "${YELLOW}⚠${NC} $dataset/protein_sequences_pairs.fasta missing"
            missing_files+=("$dataset:pairs")
        fi

        if [ -f "$mutations_file" ]; then
            count=$(grep -c '^>' "$mutations_file")
            echo -e "${GREEN}✓${NC} $dataset/protein_sequences_with_mutations.fasta ($count sequences)"
        else
            echo -e "${YELLOW}⚠${NC} $dataset/protein_sequences_with_mutations.fasta missing"
            missing_files+=("$dataset:mutations")
        fi
    done

    # Create localization output directories
    echo ""
    echo -e "${YELLOW}→${NC} Creating output directories..."
    for dataset in "${datasets[@]}"; do
        mkdir -p ../results/$dataset/localization
    done
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

# Define which file and mode each task processes (8 tasks total)
case $SLURM_ARRAY_TASK_ID in
    1)
        dataset="reduced"; file_type="pairs"; mode="Fast"
        input_file="../results/$dataset/proteins/protein_sequences_pairs.fasta"
        ;;
    2)
        dataset="reduced"; file_type="pairs"; mode="Accurate"
        input_file="../results/$dataset/proteins/protein_sequences_pairs.fasta"
        ;;
    3)
        dataset="reduced"; file_type="mutations"; mode="Fast"
        input_file="../results/$dataset/proteins/protein_sequences_with_mutations.fasta"
        ;;
    4)
        dataset="reduced"; file_type="mutations"; mode="Accurate"
        input_file="../results/$dataset/proteins/protein_sequences_with_mutations.fasta"
        ;;
    5)
        dataset="full"; file_type="pairs"; mode="Fast"
        input_file="../results/$dataset/proteins/protein_sequences_pairs.fasta"
        ;;
    6)
        dataset="full"; file_type="pairs"; mode="Accurate"
        input_file="../results/$dataset/proteins/protein_sequences_pairs.fasta"
        ;;
    7)
        dataset="full"; file_type="mutations"; mode="Fast"
        input_file="../results/$dataset/proteins/protein_sequences_with_mutations.fasta"
        ;;
    8)
        dataset="full"; file_type="mutations"; mode="Accurate"
        input_file="../results/$dataset/proteins/protein_sequences_with_mutations.fasta"
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

echo -e "${YELLOW}→${NC} Processing $dataset dataset ($file_type sequences, $mode mode)"
echo "  Input file: $(basename $input_file)"
echo ""

if [ -f "$input_file" ]; then
    # Remove outputs directory if it exists
    if [ -d "outputs/" ]; then
        echo -e "${YELLOW}→${NC} Removing existing outputs/ directory..."
        rm -rf outputs/
    fi

    # Create descriptive temporary subfolder for this specific run
    temp_subdir="../results/$dataset/localization/${dataset}_${file_type}_${mode}_temp_$$"
    mkdir -p "$temp_subdir"

    echo -e "${YELLOW}→${NC} Starting DeepLoc $mode mode for $dataset ($file_type) at $(date)"
    echo ""

    # Set GPU memory growth to avoid OOM errors
    export TF_FORCE_GPU_ALLOW_GROWTH=true
    export CUDA_VISIBLE_DEVICES=0

    deeploc2 -f "$input_file" -m "$mode" -o "$temp_subdir/" -d cuda

    # Find and move the results file
    result_file=$(find "$temp_subdir" -name "results_*.csv" | head -n 1)
    if [ -n "$result_file" ] && [ -f "$result_file" ]; then
        output_file="../results/$dataset/localization/protein_sequences_${file_type}_${mode}_results.csv"
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
    echo -e "${GREEN}✓${NC} Completed $dataset ($file_type, $mode) at $(date)"
else
    echo -e "${YELLOW}⚠${NC} Skipping $dataset ($file_type, $mode) - input file not found"
    exit 1
fi

# ============================================================================
# Verification (Task 8 Only)
# ============================================================================

if [ "$SLURM_ARRAY_TASK_ID" -eq 8 ]; then
    # Wait for file system sync
    sleep 10

    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Verification                                                ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    echo -e "${YELLOW}→${NC} Verifying DeepLoc outputs..."
    echo ""

    # Check for outputs from all tasks
    datasets=("reduced" "full")
    file_types=("pairs" "mutations")
    modes=("Fast" "Accurate")
    found_outputs=()

    for dataset in "${datasets[@]}"; do
        for file_type in "${file_types[@]}"; do
            for mode in "${modes[@]}"; do
                output_file="../results/$dataset/localization/protein_sequences_${file_type}_${mode}_results.csv"

                if [ -f "$output_file" ]; then
                    count=$(($(wc -l < "$output_file") - 1))
                    echo -e "${GREEN}✓${NC} $dataset/$(basename $output_file) ($count predictions)"
                    found_outputs+=("$output_file")
                else
                    echo -e "${YELLOW}⚠${NC} $dataset/$(basename $output_file) missing"
                fi
            done
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
        echo "Debug: Files in localization directories:"
        for dataset in "${datasets[@]}"; do
            localization_dir="../results/$dataset/localization"
            if [ -d "$localization_dir" ]; then
                echo "  $dataset/localization/:"
                ls -la "$localization_dir" | sed 's/^/    /'
            fi
        done

        exit 1
    else
        echo -e "${GREEN}✓ DeepLoc predictions completed successfully!${NC}"
        echo ""
        echo "Generated predictions:"
        echo ""

        # Group outputs by dataset
        for dataset in "${datasets[@]}"; do
            dataset_outputs=($(printf '%s\n' "${found_outputs[@]}" | grep "$dataset/localization/"))

            if [ ${#dataset_outputs[@]} -gt 0 ]; then
                echo "  └─ $dataset/localization/"
                for output in "${dataset_outputs[@]}"; do
                    echo "     ├─ $(basename $output)"
                done
                echo ""
            fi
        done

        echo "Summary: Generated ${#found_outputs[@]} prediction files"
        echo ""
        echo -e "${BLUE}Next step:${NC}"
        echo "  Analyze localization predictions for your research questions"
        echo ""
    fi
fi
