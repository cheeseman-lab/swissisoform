#!/bin/bash

#SBATCH --job-name=mutations               # Job name
#SBATCH --partition=20                     # Partition name
#SBATCH --array=1-8                        # 8 chunks for parallel processing
#SBATCH --cpus-per-task=4                  # Reduced CPUs per task
#SBATCH --mem=16G                          # Reduced memory per task
#SBATCH --time=24:00:00                    # Time limit (hrs:min:sec)
#SBATCH --output=out/mutations-%A_%a.out   # %A = job ID, %a = array task ID

# 2_analyze_mutations.sh
# Analyzes mutations in alternative isoform truncation regions with gene-level parallelization
# Reads dataset configuration from YAML file

echo "========================================================"
echo "SwissIsoform Pipeline Step 2: Analyze Mutations (Parallel)"
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
echo "========================================================"

# Dataset selection (default: hela)
DATASET="${DATASET:-hela}"
echo "Dataset: $DATASET"

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

# Only run setup checks on the first task to avoid race conditions
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    # Read dataset configuration from YAML
    echo "Reading dataset configuration..."
    CONFIG_FILE="../data/ribosome_profiling/dataset_config.yaml"

    if [ ! -f "$CONFIG_FILE" ]; then
        echo "❌ Dataset configuration file not found: $CONFIG_FILE"
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
        echo "❌ Failed to read dataset configuration"
        exit 1
    fi

    echo "Dataset configuration:"
    echo "  BED file: $DATASET_BED"
    echo "  GTF file: $DATASET_GTF"
    echo "  Gene list: $GENE_LIST"

    # Check if required input files exist
    echo ""
    echo "Checking for required input files..."

    GENOME_DIR="../data/genome_data"
    GENOME_PATH="$GENOME_DIR/GRCh38.p7.genome.fa"

    required_files=(
        "$GENOME_PATH"
        "$DATASET_GTF"
        "$DATASET_BED"
        "$GENE_LIST"
    )

    echo ""
    echo "Checking required files..."
    missing_files=false
    for file in "${required_files[@]}"; do
        if [ -f "$file" ]; then
            echo "✓ $(basename $file)"
        else
            echo "✗ $(basename $file) missing"
            missing_files=true
        fi
    done

    if [ "$missing_files" = true ]; then
        echo ""
        echo "❌ Missing required files! Run 1_cleanup_files.sh first"
        exit 1
    fi

    # Create results directory structure
    mkdir -p "../results/${DATASET}/mutations"
    mkdir -p ../results/temp/chunks
else
    # Other tasks need to read the config too
    CONFIG_FILE="../data/ribosome_profiling/dataset_config.yaml"
    read -r DATASET_BED DATASET_GTF GENE_LIST <<< $(python3 -c "
import yaml
from pathlib import Path

with open('$CONFIG_FILE') as f:
    config = yaml.safe_load(f)

dataset_config = next(ds for ds in config['datasets'] if ds['name'] == '$DATASET')
bed_file = f\"../data/ribosome_profiling/$DATASET\_isoforms_with_transcripts.bed\"
gene_list = f\"../data/ribosome_profiling/$DATASET\_isoforms_gene_list.txt\"

source_gtf_path = Path(dataset_config['source_gtf_path'])
v47_gtf = source_gtf_path.parent / f\"{source_gtf_path.stem}.v47names.gtf\"
gtf_file = str(v47_gtf) if v47_gtf.exists() else dataset_config['source_gtf_path']

print(bed_file, gtf_file, gene_list)
")
fi

# Wait for setup to complete and add staggered delays to avoid API rate limits
sleep 2
STAGGER_DELAY=$(( (SLURM_ARRAY_TASK_ID - 1) * 3 ))
echo "Adding ${STAGGER_DELAY}s staggered delay for task ${SLURM_ARRAY_TASK_ID}"
sleep $STAGGER_DELAY

# Activate conda environment
echo ""
echo "Activating conda environment..."
source ~/.bashrc
conda activate swissisoform || {
    echo "❌ Failed to activate swissisoform conda environment"
    echo "Please run: conda env create --file=../environment.yml"
    exit 1
}

# Define common paths (using dataset-specific values from config)
GENOME_PATH="../data/genome_data/GRCh38.p7.genome.fa"
ANNOTATION_PATH="$DATASET_GTF"
TRUNCATIONS_PATH="$DATASET_BED"
OUTPUT_DIR="../results/${DATASET}/mutations/chunk_${SLURM_ARRAY_TASK_ID}"
CHUNK_ID=$SLURM_ARRAY_TASK_ID
TOTAL_CHUNKS=8

# Create task-specific output directory
mkdir -p "$OUTPUT_DIR"

# Create gene chunk for this task
CHUNK_FILE="../results/temp/chunks/${DATASET}_chunk_${CHUNK_ID}.txt"
split_gene_list "$GENE_LIST" "$TOTAL_CHUNKS" "$CHUNK_ID" "$CHUNK_FILE"

# Skip if no genes in chunk
if [ ! -s "$CHUNK_FILE" ]; then
    echo "No genes in chunk ${CHUNK_ID}, exiting"
    exit 0
fi

echo "Processing ${DATASET} dataset chunk ${CHUNK_ID}"
echo "Gene list: $(wc -l < "$CHUNK_FILE") genes"

# Run analysis on the gene chunk
python3 analyze_mutations.py "$CHUNK_FILE" "$OUTPUT_DIR" \
  --genome "$GENOME_PATH" \
  --annotation "$ANNOTATION_PATH" \
  --bed "$TRUNCATIONS_PATH" \
  --sources "clinvar" \
  --impact-types "missense variant" "nonsense variant" "frameshift variant" "synonymous variant" "inframe deletion" "inframe insertion" \
  --visualize

echo "Completed ${DATASET} dataset chunk ${CHUNK_ID}"

# Clean up chunk file
rm -f "$CHUNK_FILE"

# Verification (only last task)
if [ "$SLURM_ARRAY_TASK_ID" -eq 8 ]; then
    # Wait for other tasks to finish
    echo "Waiting for all tasks to complete..."
    sleep 30
    
    echo ""
    echo "Merging results from all chunks..."

    # Create final output directory
    FINAL_OUTPUT_DIR="../results/${DATASET}/mutations"
    mkdir -p "$FINAL_OUTPUT_DIR"

    # Merge gene level results
    echo "Merging gene level results..."
    first_file=true
    for i in {1..8}; do
        chunk_file="../results/${DATASET}/mutations/chunk_${i}/gene_level_results.csv"
        if [ -f "$chunk_file" ]; then
            if [ "$first_file" = true ]; then
                # Copy first file with header
                cp "$chunk_file" "$FINAL_OUTPUT_DIR/gene_level_results.csv"
                first_file=false
                echo "  ├─ Added chunk $i with header ($(wc -l < "$chunk_file") rows)"
            else
                # Append without header
                tail -n +2 "$chunk_file" >> "$FINAL_OUTPUT_DIR/gene_level_results.csv"
                rows=$(($(wc -l < "$chunk_file") - 1))
                echo "  ├─ Added chunk $i ($rows rows)"
            fi
        else
            echo "  ├─ Warning: chunk $i gene results missing"
        fi
    done
    
    # Merge isoform level results
    echo "Merging isoform level results..."
    first_file=true
    for i in {1..8}; do
        chunk_file="../results/${DATASET}/mutations/chunk_${i}/isoform_level_results.csv"
        if [ -f "$chunk_file" ]; then
            if [ "$first_file" = true ]; then
                # Copy first file with header
                cp "$chunk_file" "$FINAL_OUTPUT_DIR/isoform_level_results.csv"
                first_file=false
                echo "  ├─ Added chunk $i with header ($(wc -l < "$chunk_file") rows)"
            else
                # Append without header
                tail -n +2 "$chunk_file" >> "$FINAL_OUTPUT_DIR/isoform_level_results.csv"
                rows=$(($(wc -l < "$chunk_file") - 1))
                echo "  ├─ Added chunk $i ($rows rows)"
            fi
        else
            echo "  ├─ Warning: chunk $i truncation results missing"
        fi
    done
    
    # Merge visualization files (move all gene directories)
    echo "Merging visualization files..."
    for i in {1..8}; do
        chunk_viz_dir="../results/${DATASET}/mutations/chunk_${i}"
        if [ -d "$chunk_viz_dir" ]; then
            # Move all gene directories to final location
            find "$chunk_viz_dir" -maxdepth 1 -type d -not -name "chunk_*" -exec mv {} "$FINAL_OUTPUT_DIR/" \; 2>/dev/null || true
        fi
    done

    echo ""
    echo "Verifying merged mutation analysis outputs..."

    expected_files=(
        "../results/${DATASET}/mutations/gene_level_results.csv"
        "../results/${DATASET}/mutations/isoform_level_results.csv"
    )

    all_files_present=true
    for file in "${expected_files[@]}"; do
        if [ -f "$file" ]; then
            count=$(($(wc -l < "$file") - 1))  # Subtract header
            echo "✓ $(basename $file) ($count rows)"
        else
            echo "✗ $(basename $file) missing"
            all_files_present=false
        fi
    done

    # Check for visualization outputs
    echo ""
    echo "Checking visualization outputs..."
    vis_count=$(find "../results/${DATASET}/mutations" -name "*.pdf" -not -path "*/chunk_*/*" 2>/dev/null | wc -l)
    echo "✓ ${DATASET} dataset visualizations: $vis_count PDFs"

    if [ "$all_files_present" = true ]; then
        echo ""
        echo "🎉 ${DATASET} dataset mutation analysis completed successfully!"
        echo ""
        echo "Generated analysis results:"
        echo "  └─ ${DATASET}/mutations/              # All alternative isoform sites"
        echo "     ├─ gene_level_results.csv          # Gene-level mutation summary"
        echo "     ├─ isoform_level_results.csv       # Detailed transcript-isoform pairs"
        echo "     └─ [gene_name]/                    # Visualizations by gene"
        echo ""
        echo "Next step:"
        echo "  Run: sbatch --export=DATASET=${DATASET} 3_generate_proteins.sh"

        # Clean up chunk directories
        echo ""
        echo "Cleaning up chunk directories..."
        rm -rf "../results/${DATASET}/mutations/chunk_"*
    else
        echo ""
        echo "❌ ${DATASET} dataset mutation analysis failed. Some output files are missing."
        echo "Check the analyze_mutations.py script for errors."
        exit 1
    fi
fi

# Final cleanup
if [ "$SLURM_ARRAY_TASK_ID" -eq 8 ]; then
    sleep 30
    rm -rf ../results/temp/chunks
fi