#!/bin/bash

#SBATCH --job-name=mutations               # Job name
#SBATCH --partition=20                     # Partition name
#SBATCH --ntasks=2                         # Run 2 parallel tasks
#SBATCH --cpus-per-task=8                  # CPUs per task (2*8=16 total)
#SBATCH --mem=64G                          # Total memory for all tasks
#SBATCH --time=12:00:00                    # Time limit (hrs:min:sec)
#SBATCH --output=out/mutations-%j.out      # Standard output log

# 2_analyze_mutations.sh
# Analyzes mutations in alternative isoform truncation regions

echo "========================================================"
echo "SwissIsoform Pipeline Step 2: Analyze Mutations"
echo "========================================================"

# Check if required input files exist
echo "Checking for required input files..."

GENOME_DIR="../data/genome_data"
RIBOPROF_DIR="../data/ribosome_profiling"

required_files=(
    "$GENOME_DIR/GRCh38.p7.genome.fa"
    "$GENOME_DIR/gencode.v25.annotation.ensembl_cleaned.gtf"
    "$RIBOPROF_DIR/truncations_cleaned.bed"
    "$RIBOPROF_DIR/gene_list.txt"
    "$RIBOPROF_DIR/gene_list_reduced.txt"
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
TRUNCATIONS_PATH="../data/ribosome_profiling/truncations_cleaned.bed"
PREFERRED_TRANSCRIPTS="../data/genome_data/hela_top_transcript.txt"

# Create results directory structure
mkdir -p ../results/reduced/mutations
mkdir -p ../results/full/mutations

echo "Starting mutation analysis with 2 parallel tasks at $(date)"

# Function to run each mutation analysis task
run_task() {
    local task_id=$1
    local gene_list=$2
    local output_dir=$3
    local task_name=$4
    local extra_args=$5
    
    echo "Task $task_id: Starting $task_name at $(date)"
    
    python3 analyze_mutations.py "$gene_list" "$output_dir" \
      --genome "$GENOME_PATH" \
      --annotation "$ANNOTATION_PATH" \
      --bed "$TRUNCATIONS_PATH" \
      --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
      --sources "clinvar" \
      --impact-types "missense variant" "nonsense variant" "frameshift variant" \
      $extra_args
    
    echo "Task $task_id: Completed $task_name at $(date)"
}

# Launch 2 tasks in parallel using background processes
run_task 1 "../data/ribosome_profiling/gene_list_reduced.txt" "../results/reduced/mutations" "reduced mutations" "--visualize" &
TASK1_PID=$!

run_task 2 "../data/ribosome_profiling/gene_list.txt" "../results/full/mutations" "full mutations" "" &
TASK2_PID=$!

# Wait for all tasks to complete
echo "Waiting for all tasks to complete..."
wait $TASK1_PID
echo "Task 1 (reduced mutations) finished"

wait $TASK2_PID
echo "Task 2 (full mutations) finished"

echo "All mutation analysis completed at $(date)"

# Verify outputs
echo ""
echo "Verifying mutation analysis outputs..."

expected_files=(
    "../results/reduced/mutations/gene_level_results.csv"
    "../results/reduced/mutations/truncation_level_results.csv"
    "../results/full/mutations/gene_level_results.csv"
    "../results/full/mutations/truncation_level_results.csv"
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

# Check for visualization outputs (only for reduced dataset)
echo ""
echo "Checking visualization outputs..."
reduced_vis_count=$(find ../results/reduced/mutations -name "*.pdf" 2>/dev/null | wc -l)
full_vis_count=$(find ../results/full/mutations -name "*.pdf" 2>/dev/null | wc -l)

echo "✓ Reduced dataset visualizations: $reduced_vis_count PDFs"
echo "✓ Full dataset visualizations: $full_vis_count PDFs"

if [ "$all_files_present" = true ]; then
    echo ""
    echo "🎉 Mutation analysis completed successfully!"
    echo ""
    echo "Generated analysis results:"
    echo "  ├─ reduced/mutations/           # Curated truncation sites"
    echo "  │  ├─ gene_level_results.csv    # Gene-level mutation summary"
    echo "  │  ├─ truncation_level_results.csv  # Detailed transcript-truncation pairs"
    echo "  │  └─ [gene_name]/              # Visualizations by gene"
    echo "  └─ full/mutations/              # All truncation sites"
    echo "     ├─ gene_level_results.csv    # Gene-level mutation summary"
    echo "     ├─ truncation_level_results.csv  # Detailed transcript-truncation pairs"
    echo "     └─ [gene_name]/              # Visualizations by gene (if any)"
    echo ""
    echo "Next step:"
    echo "  Run: sbatch 3_generate_proteins.sh"
else
    echo ""
    echo "❌ Mutation analysis failed. Some output files are missing."
    echo "Check the analyze_mutations.py script for errors."
    exit 1
fi