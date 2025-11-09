#!/bin/bash
#
# SwissIsoform Pipeline Step 3: Generate Proteins (Parallel)
#
# This script generates protein sequences using chunked parallel processing with
# pre-validated missense mutations from step 2. It processes both canonical/alternative
# pairs and mutated variants.
#
# Usage:
#   sbatch 3_generate_proteins.sh
#   sbatch --export=DATASET=hela 3_generate_proteins.sh
#   sbatch --export=DATASET=hela,SOURCES="gnomad|cosmic" 3_generate_proteins.sh
#
# Environment Variables:
#   DATASET - Dataset to process (default: hela)
#   SOURCES - Mutation databases to use (default: clinvar, pipe-separated: clinvar|gnomad|cosmic)
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
#SBATCH --partition=20                     # Partition name
#SBATCH --array=1-8                        # 8 chunks for parallel processing
#SBATCH --cpus-per-task=4                  # CPUs per task
#SBATCH --mem=64G                          # Memory per task
#SBATCH --time=24:00:00                    # Time limit (hrs:min:sec)
#SBATCH --output=out/proteins-%A_%a.out    # %A = job ID, %a = array task ID

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

# Sources selection (default: clinvar)
# Convert pipe-separated string to space-separated for command line
if [ -z "$SOURCES" ]; then
    SOURCES_ARGS=("clinvar")
else
    IFS='|' read -ra SOURCES_ARGS <<< "$SOURCES"
fi

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   SwissIsoform Pipeline Step 3: Generate Proteins            ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
echo "Dataset: $DATASET"
echo "Sources: ${SOURCES_ARGS[@]}"
echo "Impact types: missense variant (only)"
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
    MUTATIONS_FILE="../results/${DATASET}/mutations/isoform_level_results.csv"

    required_files=(
        "$GENOME_PATH"
        "$DATASET_GTF"
        "$DATASET_BED"
        "$GENE_LIST"
        "$MUTATIONS_FILE"
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
        echo -e "${RED}Missing required files!${NC}"
        echo "Run 1_cleanup_files.sh and 2_analyze_mutations.sh first"
        exit 1
    fi

    # Check that we have mutation results
    mutation_count=$(tail -n +2 "$MUTATIONS_FILE" | wc -l)
    if [ "$mutation_count" -eq 0 ]; then
        echo ""
        echo -e "${RED}✗${NC} No mutation results found!"
        echo "Run 2_analyze_mutations.sh first"
        exit 1
    fi

    echo ""
    echo -e "${GREEN}✓${NC} Found $mutation_count validated mutation results"

    # Create results directory structure
    echo ""
    echo -e "${YELLOW}→${NC} Creating output directories..."
    mkdir -p "../results/${DATASET}/proteins"
    mkdir -p ../results/temp/protein_chunks
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
echo -e "${BLUE}║  Processing Chunk ${SLURM_ARRAY_TASK_ID}                                          ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Define common paths
GENOME_PATH="../data/genome_data/GRCh38.p7.genome.fa"
ANNOTATION_PATH="$DATASET_GTF"
TRUNCATIONS_PATH="$DATASET_BED"
MUTATIONS_FILE="../results/${DATASET}/mutations/isoform_level_results.csv"
OUTPUT_DIR="../results/${DATASET}/proteins/chunk_${SLURM_ARRAY_TASK_ID}"
CHUNK_ID=$SLURM_ARRAY_TASK_ID
TOTAL_CHUNKS=8
MIN_LENGTH=10
MAX_LENGTH=100000
FORMAT="fasta,csv"

# Create task-specific output directory
mkdir -p "$OUTPUT_DIR"

# Create gene chunk for this task
CHUNK_FILE="../results/temp/protein_chunks/chunk_${CHUNK_ID}.txt"
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
echo "Configuration:"
echo "  ├─ Pre-validated mutations: $(basename $MUTATIONS_FILE)"
echo "  ├─ Sources: ${SOURCES_ARGS[@]}"
echo "  ├─ Impact types: missense variant (only)"
echo "  ├─ Length range: $MIN_LENGTH-$MAX_LENGTH amino acids"
echo "  └─ Output format: $FORMAT"
echo ""

# Generate both pairs and mutations datasets for this chunk
echo -e "${YELLOW}→${NC} Generating protein sequences..."
echo ""
python3 generate_proteins.py "$CHUNK_FILE" "$OUTPUT_DIR" \
  --mutations-file "$MUTATIONS_FILE" \
  --genome "$GENOME_PATH" \
  --annotation "$ANNOTATION_PATH" \
  --bed "$TRUNCATIONS_PATH" \
  --sources "${SOURCES_ARGS[@]}" \
  --min-length "$MIN_LENGTH" \
  --max-length "$MAX_LENGTH" \
  --format "$FORMAT" \
  -v

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
    FINAL_OUTPUT_DIR="../results/${DATASET}/proteins"
    mkdir -p "$FINAL_OUTPUT_DIR"

    # Merge pairs datasets
    echo -e "${YELLOW}→${NC} Merging pairs datasets..."
    echo ""

    # Merge FASTA files (pairs)
    echo "  ├─ Merging protein_sequences_pairs.fasta..."
    > "$FINAL_OUTPUT_DIR/protein_sequences_pairs.fasta"
    for i in {1..8}; do
        chunk_fasta="../results/${DATASET}/proteins/chunk_${i}/protein_sequences_pairs.fasta"
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
    for i in {1..8}; do
        chunk_csv="../results/${DATASET}/proteins/chunk_${i}/protein_sequences_pairs.csv"
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

    # Merge mutations datasets
    echo ""
    echo -e "${YELLOW}→${NC} Merging mutations datasets..."
    echo ""

    # Merge FASTA files (mutations)
    echo "  ├─ Merging protein_sequences_with_mutations.fasta..."
    > "$FINAL_OUTPUT_DIR/protein_sequences_with_mutations.fasta"
    for i in {1..8}; do
        chunk_fasta="../results/${DATASET}/proteins/chunk_${i}/protein_sequences_with_mutations.fasta"
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
    for i in {1..8}; do
        chunk_csv="../results/${DATASET}/proteins/chunk_${i}/protein_sequences_with_mutations.csv"
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

    # Verification
    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Verification                                                ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    expected_files=(
        "../results/${DATASET}/proteins/protein_sequences_pairs.fasta"
        "../results/${DATASET}/proteins/protein_sequences_pairs.csv"
        "../results/${DATASET}/proteins/protein_sequences_with_mutations.fasta"
        "../results/${DATASET}/proteins/protein_sequences_with_mutations.csv"
    )

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
        echo -e "${GREEN}✓ ${DATASET} dataset protein sequence generation completed!${NC}"
        echo ""
        echo "Generated datasets:"
        echo "  └─ ${DATASET}/proteins/"
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
        echo "  Run: sbatch --export=DATASET=${DATASET} scripts/4_predict_localization.sh"
        echo ""

        # Clean up chunk directories
        echo -e "${YELLOW}→${NC} Cleaning up chunk directories..."
        rm -rf "../results/${DATASET}/proteins/chunk_"*
        echo -e "${GREEN}✓${NC} Cleanup complete"
    else
        echo -e "${RED}✗ ${DATASET} dataset protein generation failed${NC}"
        echo "Some output files are missing. Check the logs for errors."
        exit 1
    fi
fi

# Final cleanup
if [ "$SLURM_ARRAY_TASK_ID" -eq 8 ]; then
    sleep 30
    rm -rf ../results/temp/protein_chunks
fi
