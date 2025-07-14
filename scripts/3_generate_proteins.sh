#!/bin/bash

#SBATCH --job-name=proteins                # Job name
#SBATCH --partition=20                     # Partition name
#SBATCH --ntasks=4                         # Run 4 parallel tasks
#SBATCH --cpus-per-task=6                  # CPUs per task (4*6=24 total)
#SBATCH --mem=64G                          # Total memory for all tasks
#SBATCH --time=24:00:00                    # Time limit (hrs:min:sec)
#SBATCH --output=out/proteins-%j.out       # Standard output log

# 3_generate_proteins.sh
# Generates protein sequence datasets from truncated transcripts

echo "======================================================="
echo "SwissIsoform Pipeline Step 3: Generate Protein Datasets"
echo "======================================================="

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
MIN_LENGTH=10
MAX_LENGTH=100000
FORMAT="fasta,csv"

# Create results directory structure
mkdir -p ../results/reduced/proteins
mkdir -p ../results/full/proteins

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

run_task 2 "../data/ribosome_profiling/gene_list_reduced.txt" "../results/reduced/proteins" "reduced mutations" '--include-mutations --impact-types "missense variant" "nonsense variant" "frameshift variant"' &
TASK2_PID=$!

run_task 3 "../data/ribosome_profiling/gene_list.txt" "../results/full/proteins" "full pairs" "--pairs-only" &
TASK3_PID=$!

run_task 4 "../data/ribosome_profiling/gene_list.txt" "../results/full/proteins" "full mutations" '--include-mutations --impact-types "missense variant" "nonsense variant" "frameshift variant"' &
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

echo "All protein sequence generation completed at $(date)"

# Verify outputs
echo ""
echo "Verifying generated datasets..."

expected_files=(
    "../results/reduced/proteins/protein_sequences.fasta"
    "../results/reduced/proteins/protein_sequences.csv"
    "../results/reduced/proteins/protein_sequences_with_mutations.fasta"
    "../results/reduced/proteins/protein_sequences_with_mutations.csv"
    "../results/full/proteins/protein_sequences.fasta"
    "../results/full/proteins/protein_sequences.csv"
    "../results/full/proteins/protein_sequences_with_mutations.fasta"
    "../results/full/proteins/protein_sequences_with_mutations.csv"
)

all_files_present=true
for file in "${expected_files[@]}"; do
    if [ -f "$file" ]; then
        if [[ "$file" == *.fasta ]]; then
            count=$(grep -c '^>' "$file")
            echo "✓ $(basename $file) ($count sequences)"
        elif [[ "$file" == *.csv ]]; then
            count=$(($(wc -l < "$file") - 1))  # Subtract header
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
    echo ""
    echo "Generated datasets:"
    echo "  ├─ reduced/proteins/            # Curated truncation sites"
    echo "  │  ├─ protein_sequences.*       # Canonical + truncated pairs"
    echo "  │  └─ protein_sequences_with_mutations.*  # With mutations"
    echo "  └─ full/proteins/               # All truncation sites"
    echo "     ├─ protein_sequences.*       # Canonical + truncated pairs"
    echo "     └─ protein_sequences_with_mutations.*  # With mutations"
    echo ""
    echo "Next step:"
    echo "  Run: sbatch 4_predict_localization.sh"
else
    echo ""
    echo "❌ Protein generation failed. Some output files are missing."
    exit 1
fi