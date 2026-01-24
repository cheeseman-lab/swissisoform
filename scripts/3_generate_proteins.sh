#!/bin/bash
#
# SwissIsoform Pipeline Step 3: Generate Proteins (Parallel)
#
# This script generates protein sequences using chunked parallel processing with
# pre-validated missense mutations from step 2. It processes both canonical/alternative
# pairs and mutated variants.
#
# Usage:
#   # Single source (backward compatible) - 8 tasks
#   sbatch --export=DATASET=hela,SOURCE=gnomad 3_generate_proteins.sh
#
#   # Multiple sources in parallel (RECOMMENDED) - AUTOMATICALLY includes base proteins!
#   # For 5 sources: (1 base + 5 sources) × 8 = 48 tasks
#   sbatch --export=DATASET=hela,SOURCES="gnomad|clinvar|cosmic|custom_bch|custom_msk",CUSTOM_PARQUETS="custom_bch:/path/bch.parquet|custom_msk:/path/msk.parquet" 3_generate_proteins.sh
#
#   # This will generate:
#   #   Tasks 1-8:   results/hela/default/proteins/ (base - canonical + alternative)
#   #   Tasks 9-16:  results/hela/gnomad/proteins/ (with gnomAD mutations)
#   #   Tasks 17-24: results/hela/clinvar/proteins/ (with ClinVar mutations)
#   #   Tasks 25-32: results/hela/cosmic/proteins/ (with COSMIC mutations)
#   #   Tasks 33-40: results/hela/custom_bch/proteins/ (with BCH mutations)
#   #   Tasks 41-48: results/hela/custom_msk/proteins/ (with MSK mutations)
#
# Environment Variables:
#   DATASET - Dataset to process (default: hela)
#   SOURCE - Single mutation source for output directory (default: gnomad)
#   SOURCES - Multiple sources for parallel processing (pipe-separated: "gnomad|clinvar|cosmic")
#             ALWAYS generates base proteins first (tasks 1-8), then sources
#             Array size: (1 + num_sources) × 8 chunks (e.g., 5 sources = 48 tasks)
#   CUSTOM_PARQUET - Path to custom parquet file for single-source mode (optional)
#   CUSTOM_PARQUETS - Mapping of sources to parquet files for multi-source mode (optional)
#                     Format: "source1:/path/file1.parquet|source2:/path/file2.parquet"
#
# Note:
#   Only missense variants are processed. The mutation CSV from step 2 must contain
#   a bed_name column and ids_<source>_missense_variant columns.
#
# Prerequisites:
#   - 1_cleanup_files.sh must have been run
#   - 2_analyze_mutations.sh must have been run
#   - Mutation results must exist for the dataset
#

#SBATCH --job-name=proteins                # Job name
#SBATCH --partition=24                     # Partition name
#SBATCH --array=1-64                       # Max array size (supports up to 7 sources + 1 base × 8 chunks)
#SBATCH --cpus-per-task=4                  # CPUs per task
#SBATCH --mem=45G                          # Memory per task (48 tasks × 45G = 2160G total)
#SBATCH --time=24:00:00                    # Time limit (hrs:min:sec)
#SBATCH --output=out/proteins-%A_%a.out    # %A = job ID, %a = array task ID

set -e  # Exit on error

# Constants
CHUNKS_PER_SOURCE=8

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

# Mode selection: "base" or "source" (default: source)
# base = generate only canonical+alternative (no mutations) -> results/DATASET/default/
# source = generate with mutations -> results/DATASET/SOURCE/
MODE="${MODE:-source}"

# ============================================================================
# Multi-Source Parallelization Setup
# ============================================================================
# ALWAYS generates base proteins (default/) as tasks 1-8
# Then generates source-specific proteins with mutations
#
# Task mapping:
#   Tasks 1-8: base proteins (default/)
#   Tasks 9-16: source 1
#   Tasks 17-24: source 2
#   ... etc

if [ -n "$SOURCES" ]; then
    # Multi-source mode: base (tasks 1-8) + sources (tasks 9+)
    IFS='|' read -ra SOURCES_ARRAY <<< "$SOURCES"
    NUM_SOURCES=${#SOURCES_ARRAY[@]}
    TOTAL_TASKS=$(( (NUM_SOURCES + 1) * CHUNKS_PER_SOURCE ))  # +1 for base

    # Exit early if this task is beyond needed range
    if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_TASKS" ]; then
        echo "Task $SLURM_ARRAY_TASK_ID not needed (only $TOTAL_TASKS tasks required: 8 base + $NUM_SOURCES sources × 8)"
        exit 0
    fi

    # Map task ID: 1-8 = base, 9+ = sources
    TASK_ID=$SLURM_ARRAY_TASK_ID

    if [ "$TASK_ID" -le "$CHUNKS_PER_SOURCE" ]; then
        # Tasks 1-8: base mode
        MODE="base"
        SOURCE="default"
        CHUNK_ID=$TASK_ID
        OUTPUT_DIR_BASE="../results/${DATASET}/default"

        IS_MERGE_TASK=false
        if [ "$CHUNK_ID" -eq 8 ]; then
            IS_MERGE_TASK=true
        fi
    else
        # Tasks 9+: source mode
        MODE="source"
        ADJUSTED_TASK_ID=$(( TASK_ID - CHUNKS_PER_SOURCE ))  # Subtract base tasks
        SOURCE_INDEX=$(( (ADJUSTED_TASK_ID - 1) / CHUNKS_PER_SOURCE ))
        CHUNK_ID=$(( ((ADJUSTED_TASK_ID - 1) % CHUNKS_PER_SOURCE) + 1 ))
        SOURCE="${SOURCES_ARRAY[$SOURCE_INDEX]}"
        OUTPUT_DIR_BASE="../results/${DATASET}/${SOURCE}"

        # Determine if this is the last task for this source (responsible for merging)
        LAST_TASK_FOR_SOURCE=$(( CHUNKS_PER_SOURCE + (SOURCE_INDEX + 1) * CHUNKS_PER_SOURCE ))
        IS_MERGE_TASK=false
        if [ "$TASK_ID" -eq "$LAST_TASK_FOR_SOURCE" ]; then
            IS_MERGE_TASK=true
        fi
    fi
else
    # Single-source mode: use SOURCE variable (backward compatibility)
    SOURCE="${SOURCE:-gnomad}"
    NUM_SOURCES=1
    TOTAL_TASKS=$CHUNKS_PER_SOURCE
    CHUNK_ID=$SLURM_ARRAY_TASK_ID
    MODE="source"

    IS_MERGE_TASK=false
    if [ "$CHUNK_ID" -eq 8 ]; then
        IS_MERGE_TASK=true
    fi

    OUTPUT_DIR_BASE="../results/${DATASET}/${SOURCE}"
fi

# Sources selection (default: clinvar)
# Convert pipe-separated string to space-separated for command line
if [ -z "$SOURCES" ]; then
    SOURCES_ARGS=("clinvar")
else
    IFS='|' read -ra SOURCES_ARGS <<< "$SOURCES"
fi

# Custom parquet files (optional)
CUSTOM_PARQUET="${CUSTOM_PARQUET:-}"
CUSTOM_PARQUETS="${CUSTOM_PARQUETS:-}"

# If using multi-source mode with custom parquets, parse the mapping
if [ -n "$CUSTOM_PARQUETS" ] && [ "$MODE" != "base" ]; then
    # Parse CUSTOM_PARQUETS into associative array
    # Format: "source1:/path1|source2:/path2"
    declare -A PARQUET_MAP
    IFS='|' read -ra PARQUET_ENTRIES <<< "$CUSTOM_PARQUETS"
    for entry in "${PARQUET_ENTRIES[@]}"; do
        source_name="${entry%%:*}"
        parquet_path="${entry#*:}"
        PARQUET_MAP["$source_name"]="$parquet_path"
    done

    # Set CUSTOM_PARQUET for this specific source if it exists in the map
    if [ -n "${PARQUET_MAP[$SOURCE]}" ]; then
        CUSTOM_PARQUET="${PARQUET_MAP[$SOURCE]}"
    else
        CUSTOM_PARQUET=""
    fi
fi

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   SwissIsoform Pipeline Step 3: Generate Proteins            ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
if [ -n "$SOURCES" ]; then
    echo "Multi-source mode: base + ${NUM_SOURCES} sources ($(( (NUM_SOURCES + 1) * 8 )) tasks total)"
    echo "  Task 1-8: base (default/)"
    task_num=9
    for src in "${SOURCES_ARRAY[@]}"; do
        echo "  Task ${task_num}-$(( task_num + 7 )): ${src}"
        task_num=$(( task_num + 8 ))
    done
    echo ""
    echo "  This task: MODE=${MODE}, SOURCE=${SOURCE}, CHUNK=${CHUNK_ID}/8"
else
    echo "Single-source mode: SOURCE=${SOURCE}"
fi
echo "Dataset: $DATASET"
echo "Output directory: ${OUTPUT_DIR_BASE}/proteins/"
if [ "$MODE" != "base" ]; then
    echo "Mutation databases: ${SOURCES_ARGS[@]}"
    if [ -n "$CUSTOM_PARQUET" ]; then
        echo "Custom parquet for ${SOURCE}: $CUSTOM_PARQUET"
    fi
    if [ -n "$CUSTOM_PARQUETS" ] && [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
        echo ""
        echo "Custom parquet mappings:"
        for entry in "${PARQUET_ENTRIES[@]}"; do
            echo "  ${entry%%:*} → ${entry#*:}"
        done
    fi
    echo "Impact types: missense variant (only)"
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
    read -r DATASET_BED DATASET_GTF GENE_LIST <<< $(python3 get_dataset_config.py "$CONFIG_FILE" "$DATASET")

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

    # Only require mutation file when MODE != base
    if [ "$MODE" != "base" ]; then
        MUTATIONS_FILE="${OUTPUT_DIR_BASE}/mutations/isoform_level_results.csv"
        required_files+=("$MUTATIONS_FILE")
    fi

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
        echo -e "${RED}Missing required files!${NC}"
        if [ "$MODE" != "base" ]; then
            echo "Run 1_cleanup_files.sh and 2_analyze_mutations.sh first"
        else
            echo "Run 1_cleanup_files.sh first"
        fi
        exit 1
    fi

    # Check that we have mutation results (only when MODE != base)
    if [ "$MODE" != "base" ]; then
        mutation_count=$(tail -n +2 "$MUTATIONS_FILE" | wc -l)
        if [ "$mutation_count" -eq 0 ]; then
            echo ""
            echo -e "${RED}✗${NC} No mutation results found!"
            echo "Run 2_analyze_mutations.sh with SOURCE=${SOURCE} first"
            exit 1
        fi

        echo ""
        echo -e "${GREEN}✓${NC} Found $mutation_count validated mutation results"
    else
        echo ""
        echo -e "${GREEN}✓${NC} Base mode: skipping mutation validation"
    fi

    # Create results directory structure
    echo ""
    echo -e "${YELLOW}→${NC} Creating output directories..."
    mkdir -p "${OUTPUT_DIR_BASE}/proteins"
    if [ "$MODE" = "base" ]; then
        mkdir -p "../results/temp/protein_chunks_${DATASET}_base"
        mkdir -p "../results/temp/protein_markers_${DATASET}_base"
    else
        mkdir -p "../results/temp/protein_chunks_${DATASET}_${SOURCE}"
        mkdir -p "../results/temp/protein_markers_${DATASET}_${SOURCE}"
    fi
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

# Wait for setup to complete and add staggered delays
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
# Protein Generation
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Processing Chunk ${CHUNK_ID}                                          ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Define common paths
GENOME_PATH="../data/genome_data/GRCh38.p7.genome.fa"
ANNOTATION_PATH="$DATASET_GTF"
TRUNCATIONS_PATH="$DATASET_BED"
OUTPUT_DIR="${OUTPUT_DIR_BASE}/proteins/chunk_${CHUNK_ID}"
TOTAL_CHUNKS=8
MIN_LENGTH=10
MAX_LENGTH=100000
FORMAT="fasta,csv"

# Set mutations file only when MODE != base
if [ "$MODE" != "base" ]; then
    MUTATIONS_FILE="${OUTPUT_DIR_BASE}/mutations/isoform_level_results.csv"
fi

# Create task-specific output directory
mkdir -p "$OUTPUT_DIR"

# Create gene chunk for this task (use different naming for base vs source mode)
if [ "$MODE" = "base" ]; then
    CHUNK_FILE="../results/temp/protein_chunks_${DATASET}_base/chunk_${CHUNK_ID}.txt"
    MARKER_DIR="../results/temp/protein_markers_${DATASET}_base"
else
    CHUNK_FILE="../results/temp/protein_chunks_${DATASET}_${SOURCE}/chunk_${CHUNK_ID}.txt"
    MARKER_DIR="../results/temp/protein_markers_${DATASET}_${SOURCE}"
fi

# Ensure temp directories exist for this source
mkdir -p "$(dirname "$CHUNK_FILE")"
mkdir -p "$MARKER_DIR"

split_gene_list "$GENE_LIST" "$TOTAL_CHUNKS" "$CHUNK_ID" "$CHUNK_FILE"

# Skip if no genes in chunk
if [ ! -s "$CHUNK_FILE" ]; then
    echo -e "${YELLOW}⚠${NC} No genes in chunk ${CHUNK_ID}, exiting"
    exit 0
fi

GENE_COUNT=$(wc -l < "$CHUNK_FILE")
echo -e "${YELLOW}→${NC} Processing ${DATASET} dataset chunk ${CHUNK_ID}"
echo "  Gene list: $GENE_COUNT genes"
echo "  Mode: $MODE"
echo ""
echo "Configuration:"
if [ "$MODE" != "base" ]; then
    echo "  ├─ Pre-validated mutations: $(basename $MUTATIONS_FILE)"
    echo "  ├─ Sources: ${SOURCES_ARGS[@]}"
    if [ -n "$CUSTOM_PARQUET" ]; then
        echo "  ├─ Custom parquet: $(basename $CUSTOM_PARQUET)"
    fi
    echo "  ├─ Impact types: missense variant (only)"
else
    echo "  ├─ Generating base proteins only (no mutations)"
fi
echo "  ├─ Length range: $MIN_LENGTH-$MAX_LENGTH amino acids"
echo "  └─ Output format: $FORMAT"
echo ""

# Generate protein sequences
echo -e "${YELLOW}→${NC} Generating protein sequences..."
echo ""

# Build command based on mode
if [ "$MODE" = "base" ]; then
    # Base mode: only canonical + alternative, no mutations
    PYTHON_CMD=(python3 generate_proteins.py "$CHUNK_FILE" "$OUTPUT_DIR"
      --genome "$GENOME_PATH"
      --annotation "$ANNOTATION_PATH"
      --bed "$TRUNCATIONS_PATH"
      --min-length "$MIN_LENGTH"
      --max-length "$MAX_LENGTH"
      --format "$FORMAT"
      --base-only
      -v)
else
    # Source mode: with mutations
    PYTHON_CMD=(python3 generate_proteins.py "$CHUNK_FILE" "$OUTPUT_DIR"
      --mutations-file "$MUTATIONS_FILE"
      --genome "$GENOME_PATH"
      --annotation "$ANNOTATION_PATH"
      --bed "$TRUNCATIONS_PATH"
      --sources "${SOURCES_ARGS[@]}"
      --min-length "$MIN_LENGTH"
      --max-length "$MAX_LENGTH"
      --format "$FORMAT"
      -v)

    # Add custom parquet path if specified
    if [ -n "$CUSTOM_PARQUET" ]; then
        PYTHON_CMD+=(--custom-parquet "$CUSTOM_PARQUET")
    fi
fi

"${PYTHON_CMD[@]}"

exit_code=$?

if [ $exit_code -ne 0 ]; then
    echo ""
    echo -e "${RED}✗${NC} Protein generation failed for chunk ${CHUNK_ID}"
    exit 1
fi

echo ""
echo -e "${GREEN}✓${NC} Completed chunk ${CHUNK_ID}"

# Clean up chunk file
rm -f "$CHUNK_FILE"

# Create completion marker for this task (MARKER_DIR already set and created earlier)
touch "$MARKER_DIR/chunk_${CHUNK_ID}.done"
echo -e "${GREEN}✓${NC} Created completion marker for chunk ${CHUNK_ID}"

# ============================================================================
# Result Merging and Verification (Last Task Per Source)
# ============================================================================

if [ "$IS_MERGE_TASK" = true ]; then
    # Wait for other chunks of this source to finish
    echo ""
    if [ "$MODE" = "base" ]; then
        echo -e "${YELLOW}→${NC} Waiting for all base protein chunks to complete..."
    else
        echo -e "${YELLOW}→${NC} Waiting for all chunks of source ${SOURCE} to complete..."
    fi

    # Wait for all completion markers for this source (chunks 1-8)
    MAX_WAIT=259200  # 72 hours maximum wait (matches job time limit)
    WAIT_INTERVAL=1800  # Check every 30 minutes
    elapsed=0

    while [ $elapsed -lt $MAX_WAIT ]; do
        all_done=true
        missing_chunks=""

        for i in $(seq 1 $CHUNKS_PER_SOURCE); do
            if [ ! -f "$MARKER_DIR/chunk_${i}.done" ]; then
                all_done=false
                missing_chunks="$missing_chunks $i"
            fi
        done

        if [ "$all_done" = true ]; then
            if [ "$MODE" = "base" ]; then
                echo -e "${GREEN}✓${NC} All base chunks completed!"
            else
                echo -e "${GREEN}✓${NC} All chunks for ${SOURCE} completed!"
            fi
            break
        fi

        if [ "$MODE" = "base" ]; then
            echo "  Waiting for base chunks:$missing_chunks (${elapsed}s elapsed)"
        else
            echo "  Waiting for ${SOURCE} chunks:$missing_chunks (${elapsed}s elapsed)"
        fi
        sleep $WAIT_INTERVAL
        elapsed=$((elapsed + WAIT_INTERVAL))
    done

    if [ "$all_done" != true ]; then
        if [ "$MODE" = "base" ]; then
            echo -e "${RED}✗${NC} Timeout waiting for base chunks to complete"
        else
            echo -e "${RED}✗${NC} Timeout waiting for ${SOURCE} chunks to complete"
        fi
        echo "Missing chunks:$missing_chunks"
        exit 1
    fi

    sleep 5  # Brief additional wait for filesystem sync

    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Merging Results                                             ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    # Create final output directory
    FINAL_OUTPUT_DIR="${OUTPUT_DIR_BASE}/proteins"
    mkdir -p "$FINAL_OUTPUT_DIR"

    # Merge pairs datasets
    echo -e "${YELLOW}→${NC} Merging pairs datasets..."
    echo ""

    # Merge FASTA files (pairs)
    echo "  ├─ Merging protein_sequences_pairs.fasta..."
    > "$FINAL_OUTPUT_DIR/protein_sequences_pairs.fasta"
    for i in $(seq 1 $CHUNKS_PER_SOURCE); do
        chunk_fasta="${OUTPUT_DIR_BASE}/proteins/chunk_${i}/protein_sequences_pairs.fasta"
        if [ -f "$chunk_fasta" ]; then
            cat "$chunk_fasta" >> "$FINAL_OUTPUT_DIR/protein_sequences_pairs.fasta"
            count=$(grep -c '^>' "$chunk_fasta" 2>/dev/null || echo 0)
            echo "  │  ├─ Added chunk $i ($count sequences)"
        else
            echo -e "  │  ├─ ${YELLOW}⚠${NC} Warning: chunk $i pairs FASTA missing"
        fi
    done

    # Merge CSV files (pairs)
    echo "  ├─ Merging protein_sequences_pairs.csv..."
    first_file=true
    for i in $(seq 1 $CHUNKS_PER_SOURCE); do
        chunk_csv="${OUTPUT_DIR_BASE}/proteins/chunk_${i}/protein_sequences_pairs.csv"
        if [ -f "$chunk_csv" ]; then
            if [ "$first_file" = true ]; then
                cp "$chunk_csv" "$FINAL_OUTPUT_DIR/protein_sequences_pairs.csv"
                first_file=false
                rows=$(wc -l < "$chunk_csv")
                echo "  │  ├─ Added chunk $i with header ($rows rows)"
            else
                tail -n +2 "$chunk_csv" >> "$FINAL_OUTPUT_DIR/protein_sequences_pairs.csv"
                rows=$(($(wc -l < "$chunk_csv") - 1))
                echo "  │  ├─ Added chunk $i ($rows rows)"
            fi
        else
            echo -e "  │  ├─ ${YELLOW}⚠${NC} Warning: chunk $i pairs CSV missing"
        fi
    done

    # Merge mutations datasets (only when MODE != base)
    if [ "$MODE" != "base" ]; then
        echo ""
        echo -e "${YELLOW}→${NC} Merging mutations datasets..."
        echo ""

        # Merge FASTA files (mutations)
        echo "  ├─ Merging protein_sequences_with_mutations.fasta..."
        > "$FINAL_OUTPUT_DIR/protein_sequences_with_mutations.fasta"
        for i in $(seq 1 $CHUNKS_PER_SOURCE); do
            chunk_fasta="${OUTPUT_DIR_BASE}/proteins/chunk_${i}/protein_sequences_with_mutations.fasta"
        if [ -f "$chunk_fasta" ]; then
            cat "$chunk_fasta" >> "$FINAL_OUTPUT_DIR/protein_sequences_with_mutations.fasta"
            count=$(grep -c '^>' "$chunk_fasta" 2>/dev/null || echo 0)
            echo "  │  ├─ Added chunk $i ($count sequences)"
        else
            echo -e "  │  ├─ ${YELLOW}⚠${NC} Warning: chunk $i mutations FASTA missing"
        fi
    done

        # Merge CSV files (mutations)
        echo "  ├─ Merging protein_sequences_with_mutations.csv..."
        first_file=true
        for i in $(seq 1 $CHUNKS_PER_SOURCE); do
            chunk_csv="${OUTPUT_DIR_BASE}/proteins/chunk_${i}/protein_sequences_with_mutations.csv"
        if [ -f "$chunk_csv" ]; then
            if [ "$first_file" = true ]; then
                cp "$chunk_csv" "$FINAL_OUTPUT_DIR/protein_sequences_with_mutations.csv"
                first_file=false
                rows=$(wc -l < "$chunk_csv")
                echo "  │  ├─ Added chunk $i with header ($rows rows)"
            else
                tail -n +2 "$chunk_csv" >> "$FINAL_OUTPUT_DIR/protein_sequences_with_mutations.csv"
                rows=$(($(wc -l < "$chunk_csv") - 1))
                echo "  │  ├─ Added chunk $i ($rows rows)"
            fi
        else
            echo -e "  │  ├─ ${YELLOW}⚠${NC} Warning: chunk $i mutations CSV missing"
            fi
        done
    fi  # End MODE != base check for mutations datasets

    # Verification
    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Verification                                                ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    expected_files=(
        "${OUTPUT_DIR_BASE}/proteins/protein_sequences_pairs.fasta"
        "${OUTPUT_DIR_BASE}/proteins/protein_sequences_pairs.csv"
    )

    # Only check for mutations files when MODE != base
    if [ "$MODE" != "base" ]; then
        expected_files+=(
            "${OUTPUT_DIR_BASE}/proteins/protein_sequences_with_mutations.fasta"
            "${OUTPUT_DIR_BASE}/proteins/protein_sequences_with_mutations.csv"
        )
    fi

    all_files_present=true
    for file in "${expected_files[@]}"; do
        if [ -f "$file" ]; then
            if [[ "$file" == *.fasta ]]; then
                count=$(grep -c '^>' "$file" 2>/dev/null || echo 0)
                echo -e "${GREEN}✓${NC} $(basename $file) ($count sequences)"
            elif [[ "$file" == *.csv ]]; then
                count=$(($(wc -l < "$file") - 1))
                echo -e "${GREEN}✓${NC} $(basename $file) ($count rows)"
            fi
        else
            echo -e "${RED}✗${NC} $(basename $file) missing"
            all_files_present=false
        fi
    done

    # Summary
    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Protein Generation Complete                                 ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    if [ "$all_files_present" = true ]; then
        if [ "$MODE" = "base" ]; then
            echo -e "${GREEN}✓ ${DATASET}/default protein sequence generation completed!${NC}"
            echo ""
            echo "Generated datasets:"
            echo "  └─ ${DATASET}/default/proteins/"
            echo "     ├─ protein_sequences_pairs.fasta"
            echo "     └─ protein_sequences_pairs.csv"
            echo ""
            echo "Performance benefits:"
            echo "  ├─ Parallel processing with 8 chunks"
            echo "  └─ Base proteins can be reused across all sources"
            echo ""
            echo -e "${BLUE}Next step:${NC}"
            echo "  Run: sbatch --export=DATASET=${DATASET},MODE=base scripts/4_predict_localization.sh"
        else
            echo -e "${GREEN}✓ ${DATASET}/${SOURCE} protein sequence generation completed!${NC}"
            echo ""
            echo "Generated datasets:"
            echo "  └─ ${DATASET}/${SOURCE}/proteins/"
            echo "     ├─ protein_sequences_pairs.* (canonical + alternative pairs)"
            echo "     └─ protein_sequences_with_mutations.* (with mutations applied)"
            echo ""
            echo "Performance benefits:"
            echo "  ├─ Parallel processing with 8 chunks"
            echo "  ├─ No mutation re-fetching (used cached results)"
            echo "  ├─ No mutation re-validation (used pre-validated impacts)"
            echo "  └─ Direct mutation application from step 2 results"
            echo ""
            echo -e "${BLUE}Next step:${NC}"
            if [ "$NUM_SOURCES" -gt 1 ]; then
                echo "  Multi-source job will automatically process all sources in parallel"
                echo "  Then run: sbatch --export=DATASET=${DATASET},SOURCES=\"${SOURCES}\" scripts/4_predict_localization.sh"
            else
                echo "  Run: sbatch --export=DATASET=${DATASET},SOURCE=${SOURCE} scripts/4_predict_localization.sh"
            fi
        fi
        echo ""

        # Clean up chunk directories
        echo -e "${YELLOW}→${NC} Cleaning up chunk directories..."
        rm -rf "${OUTPUT_DIR_BASE}/proteins/chunk_"*
        echo -e "${GREEN}✓${NC} Cleanup complete"
    else
        if [ "$MODE" = "base" ]; then
            echo -e "${RED}✗ ${DATASET}/default protein generation failed${NC}"
        else
            echo -e "${RED}✗ ${DATASET}/${SOURCE} protein generation failed${NC}"
        fi
        echo "Some output files are missing. Check the logs for errors."
        exit 1
    fi
fi

# Final cleanup (only remove this source's files, not other concurrent sources)
if [ "$IS_MERGE_TASK" = true ]; then
    sleep 5
    if [ "$MODE" = "base" ]; then
        rm -rf "../results/temp/protein_chunks_${DATASET}_base"
    else
        rm -rf "../results/temp/protein_chunks_${DATASET}_${SOURCE}"
    fi
    rm -rf "$MARKER_DIR"
fi
