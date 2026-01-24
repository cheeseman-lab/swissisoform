#!/bin/bash
#
# SwissIsoform Pipeline Step 4: Predict Localization (Parallel)
#
# This script generates subcellular localization predictions using DeepLoc.
# It processes protein sequences in parallel with both Fast and Accurate modes.
#
# Usage:
#   # Single source (backward compatible) - 4 tasks
#   sbatch --export=DATASET=hela,SOURCE=gnomad 4_predict_localization.sh
#
#   # Multiple sources in parallel (RECOMMENDED) - AUTOMATICALLY includes base pairs!
#   # For 5 sources: (2 base + 5 sources × 2) = 12 tasks
#   sbatch --export=DATASET=hela,SOURCES="gnomad|clinvar|cosmic|custom_bch|custom_msk" 4_predict_localization.sh
#
#   # This will generate:
#   #   Task 1-2:  results/hela/default/localization/ (base pairs Fast + Accurate)
#   #   Task 3-4:  results/hela/gnomad/localization/ (mutations Fast + Accurate)
#   #   Task 5-6:  results/hela/clinvar/localization/ (mutations Fast + Accurate)
#   #   Task 7-8:  results/hela/cosmic/localization/ (mutations Fast + Accurate)
#   #   Task 9-10: results/hela/custom_bch/localization/ (mutations Fast + Accurate)
#   #   Task 11-12: results/hela/custom_msk/localization/ (mutations Fast + Accurate)
#
# Environment Variables:
#   DATASET - Dataset to process (default: hela)
#   SOURCE - Single mutation source (default: gnomad)
#   SOURCES - Multiple sources for parallel processing (pipe-separated: "gnomad|clinvar|cosmic")
#             ALWAYS predicts base pairs first (tasks 1-2), then source mutations
#             Array size: 2 + (num_sources × 2) tasks (e.g., 5 sources = 12 tasks)
#
# Prerequisites:
#   - 3_generate_proteins.sh must have been run for the dataset
#   - Protein sequence files must exist
#   - DeepLoc conda environment must be installed
#

#SBATCH --job-name=deeploc                 # Job name
#SBATCH --partition=nvidia-A6000-20        # GPU partition
#SBATCH --array=1-16                       # Max array size (2 base + up to 7 sources × 2)
#SBATCH --cpus-per-task=4                  # CPUs per task
#SBATCH --mem=36G                          # Memory per task
#SBATCH --gres=gpu:1                       # 1 GPU per task
#SBATCH --time=36:00:00                    # Time limit
#SBATCH --output=out/deeploc-%A_%a.out     # %A = job ID, %a = array task ID

set -e  # Exit on error

# Constants
TASKS_FOR_BASE=2  # Fast + Accurate for base pairs
TASKS_PER_SOURCE=2  # Fast + Accurate for mutations only

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

# ============================================================================
# Multi-Source Parallelization Setup
# ============================================================================
# ALWAYS predicts base pairs first (tasks 1-2: Fast + Accurate)
# Then predicts mutations for each source (2 tasks per source)
#
# Task mapping for 5 sources:
#   Tasks 1-2: base pairs (default/localization/)
#   Tasks 3-4: source 1 mutations only (Fast + Accurate)
#   Tasks 5-6: source 2 mutations only (Fast + Accurate)
#   ... etc

if [ -n "$SOURCES" ]; then
    # Multi-source mode: base (tasks 1-2) + sources (tasks 3+)
    IFS='|' read -ra SOURCES_ARRAY <<< "$SOURCES"
    NUM_SOURCES=${#SOURCES_ARRAY[@]}
    TOTAL_TASKS=$(( TASKS_FOR_BASE + NUM_SOURCES * TASKS_PER_SOURCE ))

    # Exit early if this task is beyond needed range
    if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_TASKS" ]; then
        echo "Task $SLURM_ARRAY_TASK_ID not needed (only $TOTAL_TASKS tasks required: 2 base + $NUM_SOURCES sources × 2)"
        exit 0
    fi

    TASK_ID=$SLURM_ARRAY_TASK_ID

    if [ "$TASK_ID" -le "$TASKS_FOR_BASE" ]; then
        # Tasks 1-2: base pairs only
        MODE="base"
        SOURCE="default"
        PROTEINS_DIR="../results/${DATASET}/default/proteins"
        LOCALIZATION_DIR="../results/${DATASET}/default/localization"

        # Task 1 = Fast, Task 2 = Accurate
        if [ "$TASK_ID" -eq 1 ]; then
            DEEPLOC_MODE="Fast"
        else
            DEEPLOC_MODE="Accurate"
        fi

        FILE_TYPE="pairs"
        INPUT_FILE="${PROTEINS_DIR}/protein_sequences_pairs.fasta"
    else
        # Tasks 3+: source mutations only (no pairs prediction)
        MODE="source"
        ADJUSTED_TASK_ID=$(( TASK_ID - TASKS_FOR_BASE ))
        SOURCE_INDEX=$(( (ADJUSTED_TASK_ID - 1) / TASKS_PER_SOURCE ))
        TASK_WITHIN_SOURCE=$(( ((ADJUSTED_TASK_ID - 1) % TASKS_PER_SOURCE) + 1 ))

        SOURCE="${SOURCES_ARRAY[$SOURCE_INDEX]}"
        PROTEINS_DIR="../results/${DATASET}/${SOURCE}/proteins"
        LOCALIZATION_DIR="../results/${DATASET}/${SOURCE}/localization"

        # Task 1 = Fast, Task 2 = Accurate (within source)
        if [ "$TASK_WITHIN_SOURCE" -eq 1 ]; then
            DEEPLOC_MODE="Fast"
        else
            DEEPLOC_MODE="Accurate"
        fi

        FILE_TYPE="mutations"
        INPUT_FILE="${PROTEINS_DIR}/protein_sequences_with_mutations.fasta"
    fi
else
    # Single-source mode (backward compatible)
    SOURCE="${SOURCE:-gnomad}"
    NUM_SOURCES=1
    MODE="source"

    PROTEINS_DIR="../results/${DATASET}/${SOURCE}/proteins"
    LOCALIZATION_DIR="../results/${DATASET}/${SOURCE}/localization"

    # Old behavior: 4 tasks (pairs Fast/Accurate + mutations Fast/Accurate)
    case $SLURM_ARRAY_TASK_ID in
        1)
            FILE_TYPE="pairs"; DEEPLOC_MODE="Fast"
            INPUT_FILE="${PROTEINS_DIR}/protein_sequences_pairs.fasta"
            ;;
        2)
            FILE_TYPE="pairs"; DEEPLOC_MODE="Accurate"
            INPUT_FILE="${PROTEINS_DIR}/protein_sequences_pairs.fasta"
            ;;
        3)
            FILE_TYPE="mutations"; DEEPLOC_MODE="Fast"
            INPUT_FILE="${PROTEINS_DIR}/protein_sequences_with_mutations.fasta"
            ;;
        4)
            FILE_TYPE="mutations"; DEEPLOC_MODE="Accurate"
            INPUT_FILE="${PROTEINS_DIR}/protein_sequences_with_mutations.fasta"
            ;;
        *)
            echo "Single-source mode only uses tasks 1-4"
            exit 0
            ;;
    esac
fi

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   SwissIsoform Pipeline Step 4: Predict Localization         ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
if [ -n "$SOURCES" ]; then
    echo "Multi-source mode: base + ${NUM_SOURCES} sources (${TOTAL_TASKS} tasks total)"
    echo "  Task 1-2: base pairs (default/localization/) [Fast, Accurate]"
    task_num=3
    for src in "${SOURCES_ARRAY[@]}"; do
        echo "  Task ${task_num}-$(( task_num + 1 )): ${src} mutations [Fast, Accurate]"
        task_num=$(( task_num + 2 ))
    done
    echo ""
    echo "  This task: MODE=${MODE}, SOURCE=${SOURCE}, FILE=${FILE_TYPE}, DEEPLOC=${DEEPLOC_MODE}"
else
    echo "Single-source mode: SOURCE=${SOURCE}"
    echo "  This task: FILE=${FILE_TYPE}, DEEPLOC=${DEEPLOC_MODE}"
fi
echo "Dataset: $DATASET"
echo "Output directory: ${LOCALIZATION_DIR}/"
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

    pairs_file="${PROTEINS_DIR}/protein_sequences_pairs.fasta"

    if [ -f "$pairs_file" ]; then
        count=$(grep -c '^>' "$pairs_file")
        echo -e "${GREEN}✓${NC} protein_sequences_pairs.fasta ($count sequences)"
    else
        echo -e "${RED}✗${NC} protein_sequences_pairs.fasta missing"
        missing_files+=("pairs")
    fi

    # Only check for mutations file when MODE != base
    if [ "$MODE" != "base" ]; then
        mutations_file="${PROTEINS_DIR}/protein_sequences_with_mutations.fasta"
        if [ -f "$mutations_file" ]; then
            count=$(grep -c '^>' "$mutations_file")
            echo -e "${GREEN}✓${NC} protein_sequences_with_mutations.fasta ($count sequences)"
        else
            echo -e "${RED}✗${NC} protein_sequences_with_mutations.fasta missing"
            missing_files+=("mutations")
        fi
    fi

    if [ ${#missing_files[@]} -gt 0 ]; then
        echo ""
        echo -e "${RED}Missing required files!${NC}"
        if [ "$MODE" = "base" ]; then
            echo "Run: sbatch --export=DATASET=${DATASET},MODE=base scripts/3_generate_proteins.sh"
        else
            echo "Run: sbatch --export=DATASET=${DATASET},SOURCE=${SOURCE} scripts/3_generate_proteins.sh"
        fi
        exit 1
    fi

    # Create localization output directories
    echo ""
    echo -e "${YELLOW}→${NC} Creating output directories..."
    mkdir -p "${LOCALIZATION_DIR}"
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
# DeepLoc Prediction
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Processing Task ${SLURM_ARRAY_TASK_ID}                                            ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${YELLOW}→${NC} Processing ${FILE_TYPE} sequences in ${DEEPLOC_MODE} mode"
echo "  Input file: $(basename $INPUT_FILE)"
echo ""

if [ -f "$INPUT_FILE" ]; then
    # Remove outputs directory if it exists
    if [ -d "outputs/" ]; then
        echo -e "${YELLOW}→${NC} Removing existing outputs/ directory..."
        rm -rf outputs/
    fi

    # Create descriptive temporary subfolder for this specific run
    temp_subdir="${LOCALIZATION_DIR}/${FILE_TYPE}_${DEEPLOC_MODE}_temp_$$"
    mkdir -p "$temp_subdir"

    echo -e "${YELLOW}→${NC} Starting DeepLoc ${DEEPLOC_MODE} mode for ${FILE_TYPE} at $(date)"
    echo ""

    # Set GPU memory growth to avoid OOM errors
    export TF_FORCE_GPU_ALLOW_GROWTH=true
    export CUDA_VISIBLE_DEVICES=0

    deeploc2 -f "$INPUT_FILE" -m "$DEEPLOC_MODE" -o "$temp_subdir/" -d cuda

    # Find and move the results file
    result_file=$(find "$temp_subdir" -name "results_*.csv" | head -n 1)
    if [ -n "$result_file" ] && [ -f "$result_file" ]; then
        output_file="${LOCALIZATION_DIR}/protein_sequences_${FILE_TYPE}_${DEEPLOC_MODE}_results.csv"
        mv "$result_file" "$output_file"
        echo ""
        echo -e "${GREEN}✓${NC} Moved ${DEEPLOC_MODE} results to protein_sequences_${FILE_TYPE}_${DEEPLOC_MODE}_results.csv"
    else
        echo ""
        echo -e "${RED}✗${NC} ${DEEPLOC_MODE} results not found"
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
    echo -e "${GREEN}✓${NC} Completed ${FILE_TYPE} (${DEEPLOC_MODE}) at $(date)"
else
    echo -e "${YELLOW}⚠${NC} Skipping ${FILE_TYPE} (${DEEPLOC_MODE}) - input file not found"
    exit 1
fi

# ============================================================================
# Verification (Last Task Only)
# ============================================================================

# Only verify on the very last task in multi-source mode, or task 4 in single-source mode
IS_LAST_TASK=false
if [ -n "$SOURCES" ] && [ "$TASK_ID" -eq "$TOTAL_TASKS" ]; then
    IS_LAST_TASK=true
elif [ -z "$SOURCES" ] && [ "$SLURM_ARRAY_TASK_ID" -eq 4 ]; then
    IS_LAST_TASK=true
fi

if [ "$IS_LAST_TASK" = true ]; then
    # Wait for file system sync
    sleep 10

    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Verification Summary                                        ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    if [ -n "$SOURCES" ]; then
        echo -e "${GREEN}✓ Multi-source localization prediction completed!${NC}"
        echo ""
        echo "Generated predictions:"
        echo "  └─ ${DATASET}/"
        echo "     ├─ default/localization/"
        echo "     │  ├─ protein_sequences_pairs_Fast_results.csv"
        echo "     │  └─ protein_sequences_pairs_Accurate_results.csv"
        for src in "${SOURCES_ARRAY[@]}"; do
            echo "     ├─ ${src}/localization/"
            echo "     │  ├─ protein_sequences_mutations_Fast_results.csv"
            echo "     │  └─ protein_sequences_mutations_Accurate_results.csv"
        done
        echo ""
        echo -e "${BLUE}Next step:${NC}"
        echo "  Run: DATASET=${DATASET} SOURCES=\"${SOURCES}\" bash scripts/5_summarize_results.sh"
    else
        echo -e "${GREEN}✓ Single-source localization prediction completed for ${DATASET}/${SOURCE}!${NC}"
        echo ""
        echo "Generated predictions:"
        echo "  └─ ${DATASET}/${SOURCE}/localization/"
        echo "     ├─ protein_sequences_pairs_Fast_results.csv"
        echo "     ├─ protein_sequences_pairs_Accurate_results.csv"
        echo "     ├─ protein_sequences_mutations_Fast_results.csv"
        echo "     └─ protein_sequences_mutations_Accurate_results.csv"
        echo ""
        echo -e "${BLUE}Next step:${NC}"
        echo "  Run: DATASET=${DATASET} SOURCE=${SOURCE} bash scripts/5_summarize_results.sh"
    fi
    echo ""
fi
