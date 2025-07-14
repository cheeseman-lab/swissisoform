#!/bin/bash

#SBATCH --job-name=swissisoform_proteins   # Job name
#SBATCH --partition=20                     # Partition name
#SBATCH --ntasks=4                         # Run 4 parallel tasks
#SBATCH --cpus-per-task=6                  # CPUs per task (4*6=24 total)
#SBATCH --mem=32G                          # Total memory for all tasks
#SBATCH --time=24:00:00                    # Time limit (hrs:min:sec)
#SBATCH --output=out/proteins-%j.out       # Standard output log

# 3_generate_proteins.sh
# Generates protein sequence datasets from truncated transcripts

echo "======================================================="
echo "SwissIsoform Pipeline Step 3: Generate Protein Datasets"
echo "======================================================="

# Activate conda environment
source ~/.bashrc
conda activate swissisoform

# Define common paths
GENOME_PATH="../data/genome_data/GRCh38.p7.genome.fa"
ANNOTATION_PATH="../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf"
TRUNCATIONS_PATH="../data/ribosome_profiling/truncations_cleaned.bed"
PREFERRED_TRANSCRIPTS="../data/genome_data/hela_top_transcript.txt"
MIN_LENGTH=10
MAX_LENGTH=100000
FORMAT="fasta,csv"

# Create results directory structure
mkdir -p ../results/reduced/proteins
mkdir -p ../results/reduced/mutations
mkdir -p ../results/full/proteins
mkdir -p ../results/full/mutations

echo "Starting protein sequence generation with 4 parallel tasks at $(date)"

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
run_task 1 "../data/ribosome_profiling/gene_list_reduced.txt" "../results/reduced/proteins" "reduced pairs" "--pairs-only" &
TASK1_PID=$!

run_task 2 "../data/ribosome_profiling/gene_list_reduced.txt" "../results/reduced/mutations" "reduced mutations" '--include-mutations --impact-types "missense variant" "nonsense variant" "frameshift variant"' &
TASK2_PID=$!

run_task 3 "../data/ribosome_profiling/gene_list.txt" "../results/full/proteins" "full pairs" "--pairs-only" &
TASK3_PID=$!

run_task 4 "../data/ribosome_profiling/gene_list.txt" "../results/full/mutations" "full mutations" '--include-mutations --impact-types "missense variant" "nonsense variant" "frameshift variant"' &
TASK4_PID=$!

# Wait for all tasks to complete
echo "Waiting for all tasks to complete..."
wait $TASK1_PID
wait $TASK2_PID
wait $TASK3_PID
wait $TASK4_PID

echo "All protein sequence generation completed at $(date)"

# Verify outputs
echo ""
echo "Verifying generated datasets..."

expected_files=(
    "../results/reduced/proteins/protein_sequences.fasta"
    "../results/reduced/proteins/protein_sequences.csv"
    "../results/reduced/mutations/protein_sequences_with_mutations.fasta"
    "../results/reduced/mutations/protein_sequences_with_mutations.csv"
    "../results/full/proteins/protein_sequences.fasta"
    "../results/full/proteins/protein_sequences.csv"
    "../results/full/mutations/protein_sequences_with_mutations.fasta"
    "../results/full/mutations/protein_sequences_with_mutations.csv"
)

all_files_present=true
for file in "${expected_files[@]}"; do
    if [ -f "$file" ]; then
        if [[ "$file" == *.fasta ]]; then
            count=$(grep -c '^>' "$file")
            echo "✓ $(basename $file) ($count sequences)"
        elif [[ "$file" == *.csv ]]; then
            count=$(($(wc -l < "$file") - 1))
            echo "✓ $(basename $file) ($count rows)"
        fi
    else
        echo "✗ $(basename $file) missing"
        all_files_present=false
    fi
done

if [ "$all_files_present" = true ]; then
    echo ""
    echo "🎉 Protein sequence generation completed successfully!"
else
    echo ""
    echo "❌ Protein generation failed. Some output files are missing."
    exit 1
fi
