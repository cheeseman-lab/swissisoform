#!/bin/bash
#
# SwissIsoform Pipeline Step 2: Analyze Mutations (Parallel)
#
# This script analyzes mutations in alternative isoform truncation regions with
# gene-level parallelization. It reads dataset configuration from YAML file and
# processes chunks in parallel using SLURM array jobs.
#
# Usage:
#   # Single source (backward compatible)
#   sbatch --export=DATASET=hela,SOURCE=gnomad 2_analyze_mutations.sh
#
#   # Multiple sources in parallel (RECOMMENDED)
#   sbatch --export=DATASET=hela,SOURCES="gnomad|clinvar|cosmic" 2_analyze_mutations.sh
#
#   # Single custom mutation source
#   sbatch --export=DATASET=hela,SOURCE=custom_bch,CUSTOM_PARQUET="/path/to/bch.parquet" 2_analyze_mutations.sh
#
#   # Multiple sources including custom parquets (RECOMMENDED for multiple custom sources)
#   sbatch --export=DATASET=hela,SOURCES="gnomad|custom_bch|custom_msk",CUSTOM_PARQUETS="custom_bch:/path/to/bch.parquet|custom_msk:/path/to/msk.parquet" 2_analyze_mutations.sh
#
# Environment Variables:
#   DATASET - Dataset to process (default: hela)
#   SOURCE - Single mutation source for output directory (default: gnomad)
#   SOURCES - Multiple sources for parallel processing (pipe-separated: "gnomad|clinvar|cosmic")
#             When SOURCES is set, processes all sources in parallel (SOURCE is ignored)
#             Array size: num_sources × 8 chunks (e.g., 3 sources = 24 tasks)
#   SOURCES (mutation databases) - Which databases to query for each source (default: clinvar, pipe-separated)
#   IMPACT_TYPES - Mutation types to analyze (default: all types, pipe-separated)
#   CUSTOM_PARQUET - Path to custom parquet file for single-source mode (optional)
#   CUSTOM_PARQUETS - Mapping of sources to parquet files for multi-source mode (optional)
#                     Format: "source1:/path/file1.parquet|source2:/path/file2.parquet"
#                     Example: "custom_bch:/lab/data/bch.parquet|custom_msk:/lab/data/msk.parquet"
#
# Notes:
#   - MANE Select annotations are added automatically (no toggle needed)
#   - Output: results/DATASET/SOURCE/mutations/ with gene_level_results.csv and isoform_level_results.csv
#   - Each custom source creates its own output directory (e.g., results/hela/custom_bch/, results/hela/custom_msk/)
#
# Prerequisites:
#   - 1_cleanup_files.sh must have been run
#   - Gene lists must exist for the dataset
#   - Reference files must be present
#

#SBATCH --job-name=mutations               # Job name
#SBATCH --partition=24                     # Partition name
#SBATCH --array=1-64                       # Max array size (supports up to 8 sources × 8 chunks)
#SBATCH --cpus-per-task=4                  # CPUs per task
#SBATCH --mem=55G                          # Memory per task (40 tasks × 55G = 2200G total)
#SBATCH --time=96:00:00                    # Time limit (hrs:min:sec)
#SBATCH --output=out/mutations-%A_%a.out   # %A = job ID, %a = array task ID

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

# ============================================================================
# Multi-Source Parallelization Setup
# ============================================================================
# SOURCES (plural, pipe-separated): Determines parallelization across sources
# If SOURCES is provided, runs all sources in parallel (SOURCE is ignored)
# If SOURCES not provided, falls back to single SOURCE

if [ -n "$SOURCES" ]; then
    # Multi-source mode: parse sources and calculate task mapping
    IFS='|' read -ra SOURCES_ARRAY <<< "$SOURCES"
    NUM_SOURCES=${#SOURCES_ARRAY[@]}
    TOTAL_TASKS=$((NUM_SOURCES * CHUNKS_PER_SOURCE))

    # Exit early if this task is beyond needed range
    if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_TASKS" ]; then
        echo "Task $SLURM_ARRAY_TASK_ID not needed (only $TOTAL_TASKS tasks required for $NUM_SOURCES sources)"
        exit 0
    fi

    # Map task ID to (source_index, chunk_id)
    TASK_ID=$SLURM_ARRAY_TASK_ID
    SOURCE_INDEX=$(( (TASK_ID - 1) / CHUNKS_PER_SOURCE ))
    CHUNK_ID=$(( ((TASK_ID - 1) % CHUNKS_PER_SOURCE) + 1 ))
    SOURCE="${SOURCES_ARRAY[$SOURCE_INDEX]}"

    # Determine if this is the last task for this source (responsible for merging)
    LAST_TASK_FOR_SOURCE=$(( (SOURCE_INDEX + 1) * CHUNKS_PER_SOURCE ))
    IS_MERGE_TASK=false
    if [ "$TASK_ID" -eq "$LAST_TASK_FOR_SOURCE" ]; then
        IS_MERGE_TASK=true
    fi
else
    # Single-source mode: use SOURCE variable (backward compatibility)
    SOURCE="${SOURCE:-gnomad}"
    NUM_SOURCES=1
    TOTAL_TASKS=$CHUNKS_PER_SOURCE
    CHUNK_ID=$SLURM_ARRAY_TASK_ID

    # In single-source mode, task 8 is the merge task
    IS_MERGE_TASK=false
    if [ "$CHUNK_ID" -eq 8 ]; then
        IS_MERGE_TASK=true
    fi
fi

# Output directory: results/DATASET/SOURCE/
OUTPUT_DIR_BASE="../results/${DATASET}/${SOURCE}"

# Determine which mutation database(s) to query for this specific source
# In multi-source mode: each source only queries itself (clean output)
# In single-source mode: can query multiple databases (legacy support)
if [ -n "$SOURCES" ]; then
    # Multi-source mode: only query the current source
    # For custom sources, use "custom" as the database name
    if [[ "$SOURCE" == custom_* ]]; then
        SOURCES_ARGS=("custom")
    else
        SOURCES_ARGS=("$SOURCE")
    fi
else
    # Single-source mode: query specified databases (backward compatible)
    # This is the legacy SOURCES variable (confusingly named)
    if [ -z "$MUTATION_SOURCES" ]; then
        SOURCES_ARGS=("clinvar")
    else
        IFS='|' read -ra SOURCES_ARGS <<< "$MUTATION_SOURCES"
    fi
fi

# Impact types selection (default: all types)
# Convert pipe-separated string to space-separated for command line
if [ -z "$IMPACT_TYPES" ]; then
    IMPACT_TYPES_ARGS=("missense variant" "nonsense variant" "frameshift variant" "synonymous variant" "inframe deletion" "inframe insertion")
else
    IFS='|' read -ra IMPACT_TYPES_ARGS <<< "$IMPACT_TYPES"
fi

# Custom parquet files (optional)
# Single parquet: CUSTOM_PARQUET="/path/to/file.parquet"
# Multiple parquets: CUSTOM_PARQUETS="custom_bch:/path/to/bch.parquet|custom_msk:/path/to/msk.parquet"
CUSTOM_PARQUET="${CUSTOM_PARQUET:-}"
CUSTOM_PARQUETS="${CUSTOM_PARQUETS:-}"

# If using multi-source mode with custom parquets, parse the mapping
if [ -n "$CUSTOM_PARQUETS" ]; then
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
echo -e "${BLUE}║   SwissIsoform Pipeline Step 2: Analyze Mutations            ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
if [ "$NUM_SOURCES" -gt 1 ]; then
    echo "Multi-source mode: processing ${NUM_SOURCES} sources in parallel"
    echo "  Sources: ${SOURCES_ARRAY[@]}"
    echo "  This task: SOURCE=${SOURCE}, CHUNK=${CHUNK_ID}/8"
else
    echo "Single-source mode: SOURCE=${SOURCE}"
fi
echo "Dataset: $DATASET"
echo "Output directory: results/${DATASET}/${SOURCE}/"
echo "Mutation databases: ${SOURCES_ARGS[@]}"
echo "Impact types: ${IMPACT_TYPES_ARGS[@]}"
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
    mkdir -p "${OUTPUT_DIR_BASE}/mutations"
    mkdir -p ../results/temp/chunks
    echo -e "${GREEN}✓${NC} Directories created"
else
    # Other tasks need to read the config too
    CONFIG_FILE="../data/ribosome_profiling/dataset_config.yaml"
    read -r DATASET_BED DATASET_GTF GENE_LIST <<< $(python3 get_dataset_config.py "$CONFIG_FILE" "$DATASET")
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
echo -e "${BLUE}║  Processing Chunk ${CHUNK_ID}                                          ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Define common paths
GENOME_PATH="../data/genome_data/GRCh38.p7.genome.fa"
ANNOTATION_PATH="$DATASET_GTF"
TRUNCATIONS_PATH="$DATASET_BED"
OUTPUT_DIR="${OUTPUT_DIR_BASE}/mutations/chunk_${CHUNK_ID}"
TOTAL_CHUNKS=$CHUNKS_PER_SOURCE

# Create task-specific output directory
mkdir -p "$OUTPUT_DIR"

# Create gene chunk for this task (use DATASET_SOURCE for isolation between concurrent jobs)
CHUNK_FILE="../results/temp/chunks/${DATASET}_${SOURCE}_chunk_${CHUNK_ID}.txt"
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

# Create completion marker for this task (use DATASET_SOURCE for isolation between concurrent jobs)
MARKER_DIR="../results/temp/markers/${DATASET}_${SOURCE}"
mkdir -p "$MARKER_DIR"
touch "$MARKER_DIR/chunk_${CHUNK_ID}.done"
echo -e "${GREEN}✓${NC} Created completion marker for chunk ${CHUNK_ID}"

# ============================================================================
# Result Merging and Verification (Last Task Per Source)
# ============================================================================

if [ "$IS_MERGE_TASK" = true ]; then
    # Wait for other chunks of this source to finish
    echo ""
    echo -e "${YELLOW}→${NC} Waiting for all chunks of source ${SOURCE} to complete..."

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
            echo -e "${GREEN}✓${NC} All chunks for ${SOURCE} completed!"
            break
        fi

        echo "  Waiting for ${SOURCE} chunks:$missing_chunks (${elapsed}s elapsed)"
        sleep $WAIT_INTERVAL
        elapsed=$((elapsed + WAIT_INTERVAL))
    done

    if [ "$all_done" != true ]; then
        echo -e "${RED}✗${NC} Timeout waiting for ${SOURCE} chunks to complete"
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
    FINAL_OUTPUT_DIR="${OUTPUT_DIR_BASE}/mutations"
    mkdir -p "$FINAL_OUTPUT_DIR"

    # Move gene visualization folders from chunks to main directory
    echo ""
    echo -e "${YELLOW}→${NC} Moving gene visualization folders..."
    moved_count=0
    for i in $(seq 1 $CHUNKS_PER_SOURCE); do
        chunk_dir="${OUTPUT_DIR_BASE}/mutations/chunk_${i}"
        if [ -d "$chunk_dir" ]; then
            for gene_dir in "$chunk_dir"/*/; do
                if [ -d "$gene_dir" ]; then
                    mv "$gene_dir" "$FINAL_OUTPUT_DIR/" 2>/dev/null && ((moved_count++)) || true
                fi
            done
        fi
    done
    echo -e "${GREEN}✓${NC} Moved $moved_count gene folders"

    # Merge chunks
    echo ""
    echo -e "${YELLOW}→${NC} Merging chunks..."
    python3 merge_chunks.py \
        --output-dir "${OUTPUT_DIR_BASE}/mutations" \
        --num-chunks ${CHUNKS_PER_SOURCE}

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓${NC} Chunk merge completed"
    else
        echo -e "${RED}✗${NC} Merge failed"
        exit 1
    fi

    # Verification
    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Verification                                                ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    expected_files=(
        "${OUTPUT_DIR_BASE}/mutations/gene_level_results.csv"
        "${OUTPUT_DIR_BASE}/mutations/isoform_level_results.csv"
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
    vis_count=$(find "${OUTPUT_DIR_BASE}/mutations" -name "*.pdf" -not -path "*/chunk_*/*" 2>/dev/null | wc -l)
    echo -e "${GREEN}✓${NC} ${DATASET}/${SOURCE} visualizations: $vis_count PDFs"

    # Summary
    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Analysis Complete                                           ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    if [ "$all_files_present" = true ]; then
        echo -e "${GREEN}✓ ${DATASET}/${SOURCE} mutation analysis completed successfully!${NC}"
        echo ""
        echo "Generated analysis results:"
        echo "  └─ ${DATASET}/${SOURCE}/mutations/"
        echo "     ├─ gene_level_results.csv"
        echo "     ├─ isoform_level_results.csv (with MANE annotations)"
        echo "     ├─ isoform_level_results.xlsx"
        echo "     └─ [gene_name]/ (visualizations)"
        echo ""
        echo -e "${BLUE}Next step:${NC}"
        if [ "$NUM_SOURCES" -gt 1 ]; then
            echo "  Multi-source job will automatically process all sources in parallel"
            echo "  Then run: sbatch --export=DATASET=${DATASET},SOURCES=\"${SOURCES}\" scripts/3_generate_proteins.sh"
        else
            echo "  Run: sbatch --export=DATASET=${DATASET},SOURCE=${SOURCE} scripts/3_generate_proteins.sh"
        fi
        echo ""

        # Clean up chunk directories
        echo -e "${YELLOW}→${NC} Cleaning up chunk directories..."
        rm -rf "${OUTPUT_DIR_BASE}/mutations/chunk_"*
        echo -e "${GREEN}✓${NC} Cleanup complete"
    else
        echo -e "${RED}✗ ${DATASET}/${SOURCE} mutation analysis failed${NC}"
        echo "Some output files are missing. Check the logs for errors."
        exit 1
    fi
fi

# Final cleanup (only remove this source's files, not other concurrent sources)
if [ "$IS_MERGE_TASK" = true ]; then
    sleep 5
    rm -f ../results/temp/chunks/${DATASET}_${SOURCE}_chunk_*.txt
    rm -rf "$MARKER_DIR"
fi
