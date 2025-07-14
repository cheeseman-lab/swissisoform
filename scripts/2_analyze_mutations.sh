#!/bin/bash

#SBATCH --job-name=swissisoform_mutations  # Job name
#SBATCH --partition=20                     # Partition name
#SBATCH --ntasks=2                         # Run 2 parallel tasks (reduced and full)
#SBATCH --cpus-per-task=8                  # CPUs per task (2*8=16 total)
#SBATCH --mem=24G                          # Total memory for all tasks
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
    "$RIBOPROF_DIR/full_truncations_JL_cleaned.bed"
    "$RIBOPROF_DIR/selected_truncations_JL_cleaned.bed"
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
FULL_TRUNCATIONS_PATH="../data/ribosome_profiling/full_truncations_JL_cleaned.bed"
SELECTED_TRUNCATIONS_PATH="../data/ribosome_profiling/selected_truncations_JL_cleaned.bed"
PREFERRED_TRANSCRIPTS="../data/genome_data/hela_top_transcript.txt"

# Create results directory
mkdir -p ../results/mutations/reduced
mkdir -p ../results/mutations/full

echo "Starting mutation analysis with 2 parallel tasks at $(date)"

# Function to run each mutation analysis task
run_mutation_analysis() {
    local task_id=$1
    local gene_list=$2
    local bed_file=$3
    local output_dir=$4
    local task_name=$5
    
    echo "Task $task_id: Starting $task_name at $(date)"
    
    python3 analyze_mutations.py "$gene_list" "$output_dir" \
      --genome "$GENOME_PATH" \
      --annotation "$ANNOTATION_PATH" \
      --bed "$bed_file" \
      --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
      --visualize \
      --sources "clinvar" \
      --impact-types "missense variant" "nonsense variant" "frameshift variant"
    
    echo "Task $task_id: Completed $task_name at $(date)"
}

# Launch 2 tasks in parallel using background processes
run_mutation_analysis 1 "../data/ribosome_profiling/gene_list_reduced.txt" "$SELECTED_TRUNCATIONS_PATH" "../results/mutations/reduced" "reduced dataset" &
TASK1_PID=$!

run_mutation_analysis 2 "../data/ribosome_profiling/gene_list.txt" "$FULL_TRUNCATIONS_PATH" "../results/mutations/full" "full dataset" &
TASK2_PID=$!

# Wait for all tasks to complete
echo "Waiting for all tasks to complete..."
wait $TASK1_PID
echo "Task 1 (reduced dataset) finished"

wait $TASK2_PID
echo "Task 2 (full dataset) finished"

echo "All mutation analysis completed at $(date)"

# Verify outputs
echo ""
echo "Verifying mutation analysis outputs..."

expected_files=(
    "../results/mutations/reduced/gene_level_results.csv"
    "../results/mutations/reduced/truncation_level_results.csv"
    "../results/mutations/full/gene_level_results.csv"
    "../results/mutations/full/truncation_level_results.csv"
)

all_files_present=true
for file in "${expected_files[@]}"; do
    if [ -f "$file" ]; then
        if [[ "$file" == *.csv ]]; then
            count=$(($(wc -l < "$file") - 1))  # Subtract header
            echo "✓ $(basename $file) ($count rows)"
        else
            echo "✓ $(basename $file)"
        fi
    else
        echo "✗ $(basename $file) missing"
        all_files_present=false
    fi
done

# Check for visualization outputs (optional)
echo ""
echo "Checking visualization outputs..."
reduced_vis_count=$(find ../results/mutations/reduced -name "*.pdf" 2>/dev/null | wc -l)
full_vis_count=$(find ../results/mutations/full -name "*.pdf" 2>/dev/null | wc -l)

echo "✓ Reduced dataset visualizations: $reduced_vis_count PDFs"
echo "✓ Full dataset visualizations: $full_vis_count PDFs"

if [ "$all_files_present" = true ]; then
    echo ""
    echo "🎉 Mutation analysis completed successfully!"
    echo ""
    echo "Generated analysis results:"
    echo "  ├─ reduced/                     # Curated truncation sites"
    echo "  │  ├─ gene_level_results.csv    # Gene-level mutation summary"
    echo "  │  ├─ truncation_level_results.csv  # Detailed transcript-truncation pairs"
    echo "  │  └─ [gene_name]/              # Visualizations by gene"
    echo "  └─ full/                        # All truncation sites"
    echo "     ├─ gene_level_results.csv    # Gene-level mutation summary"
    echo "     ├─ truncation_level_results.csv  # Detailed transcript-truncation pairs"
    echo "     └─ [gene_name]/              # Visualizations by gene"
    echo ""
    echo "Next step:"
    echo "  Run: sbatch 3_generate_proteins.sh"
else
    echo ""
    echo "❌ Mutation analysis failed. Some output files are missing."
    echo "Check the analyze_mutations.py script for errors."
    exit 1
fi