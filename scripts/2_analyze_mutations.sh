#!/bin/bash

#SBATCH --job-name=mutations               # Job name
#SBATCH --partition=20                     # Partition name
#SBATCH --array=1-8                        # 8 chunks for full dataset only
#SBATCH --cpus-per-task=4                  # Reduced CPUs per task
#SBATCH --mem=16G                          # Reduced memory per task
#SBATCH --time=24:00:00                    # Time limit (hrs:min:sec)
#SBATCH --output=out/mutations-%A_%a.out   # %A = job ID, %a = array task ID

# 2_analyze_mutations_parallel.sh
# Analyzes mutations in alternative isoform truncation regions with gene-level parallelization

echo "========================================================"
echo "SwissIsoform Pipeline Step 2: Analyze Mutations (Parallel)"
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
echo "========================================================"

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
    # Check if required input files exist
    echo "Checking for required input files..."

    GENOME_DIR="../data/genome_data"
    RIBOPROF_DIR="../data/ribosome_profiling"

    required_files=(
        "$GENOME_DIR/GRCh38.p7.genome.fa"
        "$GENOME_DIR/gencode.v25.annotation.ensembl_cleaned.gtf"
        "$RIBOPROF_DIR/isoforms_with_transcripts.bed"
        "$RIBOPROF_DIR/isoforms_gene_list.txt"
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
    mkdir -p ../results/full/mutations
    mkdir -p ../results/temp/chunks
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

# Define common paths
GENOME_PATH="../data/genome_data/GRCh38.p7.genome.fa"
ANNOTATION_PATH="../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf"
TRUNCATIONS_PATH="../data/ribosome_profiling/isoforms_with_transcripts.bed"

# All tasks process full dataset chunks
DATASET="full"
GENE_LIST="../data/ribosome_profiling/isoforms_gene_list.txt"
OUTPUT_DIR="../results/full/mutations/chunk_${SLURM_ARRAY_TASK_ID}"
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

echo "Processing full dataset chunk ${CHUNK_ID}"
echo "Gene list: $(wc -l < "$CHUNK_FILE") genes"

# Run analysis on the gene chunk (same parameters as original)
python3 analyze_mutations.py "$CHUNK_FILE" "$OUTPUT_DIR" \
  --genome "$GENOME_PATH" \
  --annotation "$ANNOTATION_PATH" \
  --bed "$TRUNCATIONS_PATH" \
  --sources "clinvar" \
  --impact-types "missense variant" "nonsense variant" "frameshift variant" "synonymous variant" "inframe deletion" "inframe insertion" \
  --visualize

echo "Completed full dataset chunk ${CHUNK_ID}"

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
    FINAL_OUTPUT_DIR="../results/full/mutations"
    mkdir -p "$FINAL_OUTPUT_DIR"
    
    # Merge gene level results
    echo "Merging gene level results..."
    first_file=true
    for i in {1..8}; do
        chunk_file="../results/full/mutations/chunk_${i}/gene_level_results.csv"
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
    
    # Merge truncation level results
    echo "Merging truncation level results..."
    first_file=true
    for i in {1..8}; do
        chunk_file="../results/full/mutations/chunk_${i}/isoform_level_results.csv"
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
        chunk_viz_dir="../results/full/mutations/chunk_${i}"
        if [ -d "$chunk_viz_dir" ]; then
            # Move all gene directories to final location
            find "$chunk_viz_dir" -maxdepth 1 -type d -not -name "chunk_*" -exec mv {} "$FINAL_OUTPUT_DIR/" \; 2>/dev/null || true
        fi
    done
    
    echo ""
    echo "Verifying merged mutation analysis outputs..."

    expected_files=(
        "../results/full/mutations/gene_level_results.csv"
        "../results/full/mutations/isoform_level_results.csv"
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
    vis_count=$(find ../results/full/mutations -name "*.pdf" -not -path "*/chunk_*/*" 2>/dev/null | wc -l)
    echo "✓ Full dataset visualizations: $vis_count PDFs"

    if [ "$all_files_present" = true ]; then
        echo ""
        echo "🎉 Full dataset mutation analysis completed successfully!"
        echo ""
        echo "Generated analysis results:"
        echo "  └─ full/mutations/              # All truncation sites"
        echo "     ├─ gene_level_results.csv    # Gene-level mutation summary"
        echo "     ├─ isoform_level_results.csv  # Detailed transcript-truncation pairs"
        echo "     └─ [gene_name]/              # Visualizations by gene"
        echo ""
        echo "Next step:"
        echo "  Run: sbatch 3_generate_proteins.sh"
        
        # Clean up chunk directories
        echo ""
        echo "Cleaning up chunk directories..."
        rm -rf ../results/full/mutations/chunk_*
    else
        echo ""
        echo "❌ Full dataset mutation analysis failed. Some output files are missing."
        echo "Check the analyze_mutations.py script for errors."
        exit 1
    fi
fi

# Final cleanup
if [ "$SLURM_ARRAY_TASK_ID" -eq 8 ]; then
    sleep 30
    rm -rf ../results/temp/chunks
fi