#!/bin/bash

#SBATCH --job-name=translate              # Job name
#SBATCH --partition=20                    # Partition name
#SBATCH --ntasks=4                        # Run 4 parallel tasks
#SBATCH --cpus-per-task=6                 # CPUs per task (4*6=24 total)
#SBATCH --mem=32G                         # Total memory for all tasks
#SBATCH --time=24:00:00                   # Time limit (hrs:min:sec)
#SBATCH --output=out/translate-%j.out     # Standard output log

# Activate conda environment (adjust path as needed)
source ~/.bashrc
conda activate swissisoform

# Define common paths
GENOME_PATH="../data/genome_data/GRCh38.p7.genome.fa"
ANNOTATION_PATH="../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf"
TRUNCATIONS_PATH="../data/ribosome_profiling/full_truncations_JL_cleaned.bed"
PREFERRED_TRANSCRIPTS="../data/genome_data/hela_top_transcript.txt"
MIN_LENGTH=10
MAX_LENGTH=100000
FORMAT="fasta,csv"

# Create results directory outside of scripts
mkdir -p ../results/reduced
mkdir -p ../results/full

echo "Starting protein sequence translation with 4 parallel tasks at $(date)"

# Function to run each task
run_task() {
    local task_id=$1
    local gene_list=$2
    local output_dir=$3
    local task_name=$4
    local extra_args=$5
    
    echo "Task $task_id: Starting $task_name at $(date)"
    
    python3 translate.py "$gene_list" "$output_dir" \
      --genome "$GENOME_PATH" \
      --annotation "$ANNOTATION_PATH" \
      --bed "$TRUNCATIONS_PATH" \
      --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
      --min-length "$MIN_LENGTH" \
      --max-length "$MAX_LENGTH" \
      --format "$FORMAT" \
      --include-canonical \
      $extra_args
    
    echo "Task $task_id: Completed $task_name at $(date)"
}

# Launch 4 tasks in parallel using background processes
run_task 1 "../data/ribosome_profiling/gene_list_reduced.txt" "../results/reduced" "reduced pairs" "--pairs-only" &
TASK1_PID=$!

run_task 2 "../data/ribosome_profiling/gene_list_reduced.txt" "../results/reduced" "reduced mutations" '--include-mutations --impact-types "missense variant" "nonsense variant" "frameshift variant"' &
TASK2_PID=$!

run_task 3 "../data/ribosome_profiling/gene_list.txt" "../results/full" "full pairs" "--pairs-only" &
TASK3_PID=$!

run_task 4 "../data/ribosome_profiling/gene_list.txt" "../results/full" "full mutations" '--include-mutations --impact-types "missense variant" "nonsense variant" "frameshift variant"' &
TASK4_PID=$!

# Wait for all tasks to complete
echo "Waiting for all tasks to complete..."
wait $TASK1_PID
echo "Task 1 (reduced pairs) finished"

wait $TASK2_PID
echo "Task 2 (reduced mutations) finished"

wait $TASK3_PID
echo "Task 3 (full pairs) finished"

wait $TASK4_PID
echo "Task 4 (full mutations) finished"

echo "All protein sequence translation completed at $(date)"

# Print summary of results
echo ""
echo "=== Results Summary ==="
echo "Results saved to ../results/ directory:"
echo "  ├─ reduced/                     # Reduced dataset (both standard and mutation variants)"
echo "  └─ full/                        # Full dataset (both standard and mutation variants)"

echo ""
echo "Output files in each directory:"
echo "  Standard pairs (canonical + truncated):"
echo "  ├─ protein_sequences.fasta"
echo "  ├─ protein_sequences.csv"
echo "  Mutation variants (canonical + truncated + mutated):"
echo "  ├─ protein_sequences_with_mutations.fasta"
echo "  └─ protein_sequences_with_mutations.csv"

echo ""
echo "Checking output files:"
for dir in "reduced" "full"; do
    echo "  $dir directory:"
    if [ -f "../results/$dir/protein_sequences.fasta" ]; then
        count=$(grep -c '^>' "../results/$dir/protein_sequences.fasta")
        echo "    ✓ protein_sequences.fasta ($count sequences)"
    else
        echo "    ✗ protein_sequences.fasta missing"
    fi
    
    if [ -f "../results/$dir/protein_sequences_with_mutations.fasta" ]; then
        count=$(grep -c '^>' "../results/$dir/protein_sequences_with_mutations.fasta")
        echo "    ✓ protein_sequences_with_mutations.fasta ($count sequences)"
    else
        echo "    ✗ protein_sequences_with_mutations.fasta missing"
    fi
done