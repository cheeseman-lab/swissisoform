#!/bin/bash

#SBATCH --job-name=swissisoform_mutations  # Job name
#SBATCH --partition=20                     # Partition
#SBATCH --ntasks=2                         # Two tasks in parallel
#SBATCH --cpus-per-task=8                  # 8 CPUs per task
#SBATCH --mem=24G                          # Total memory
#SBATCH --time=12:00:00                    # Time limit
#SBATCH --output=out/mutations-%j.out      # Output log

# 2_analyze_mutations.sh
# Mutation analysis of truncation regions in alternative isoforms

echo "========================================================"
echo "SwissIsoform Pipeline Step 2: Analyze Mutations"
echo "========================================================"

echo "Checking input files..."

GENOME_DIR="../data/genome_data"
RIBOPROF_DIR="../data/ribosome_profiling"

required_files=(
    "$GENOME_DIR/GRCh38.p7.genome.fa"
    "$GENOME_DIR/gencode.v25.annotation.ensembl_cleaned.gtf"
    "$RIBOPROF_DIR/truncations_cleaned.bed"
    "$RIBOPROF_DIR/gene_list.txt"
    "$RIBOPROF_DIR/gene_list_reduced.txt"
)

missing_files=false
for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "✗ $(basename $file) missing"
        missing_files=true
    else
        echo "✓ $(basename $file)"
    fi
done

if [ "$missing_files" = true ]; then
    echo "❌ Required files missing. Run 1_cleanup_files.sh first."
    exit 1
fi

echo ""
echo "Activating swissisoform conda environment..."
source ~/.bashrc
conda activate swissisoform || {
    echo "❌ Failed to activate environment. Create it with:"
    echo "conda env create --file=../environment.yml"
    exit 1
}

GENOME_PATH="$GENOME_DIR/GRCh38.p7.genome.fa"
ANNOTATION_PATH="$GENOME_DIR/gencode.v25.annotation.ensembl_cleaned.gtf"
TRUNCATIONS_PATH="$RIBOPROF_DIR/truncations_cleaned.bed"
PREFERRED_TRANSCRIPTS="$GENOME_DIR/hela_top_transcript.txt"

mkdir -p ../results/mutations/reduced ../results/mutations/full

echo "Starting parallel mutation analysis at $(date)"

run_mutation_analysis() {
    local gene_list=$1
    local output_dir=$2
    local label=$3

    echo "[$label] Analysis started at $(date)"

    python3 analyze_mutations.py "$gene_list" "$output_dir" \
        --genome "$GENOME_PATH" \
        --annotation "$ANNOTATION_PATH" \
        --bed "$TRUNCATIONS_PATH" \
        --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
        --visualize \
        --sources "clinvar" \
        --impact-types "missense variant" "nonsense variant" "frameshift variant"

    echo "[$label] Analysis completed at $(date)"
}

run_mutation_analysis "$RIBOPROF_DIR/gene_list_reduced.txt" "../results/reduced/mutations" "Reduced" &
PID1=$!

run_mutation_analysis "$RIBOPROF_DIR/gene_list.txt" "../results/full/mutations" "Full" &
PID2=$!

echo "Waiting for both tasks to complete..."
wait $PID1
echo "Reduced dataset completed"
wait $PID2
echo "Full dataset completed"

echo ""
echo "Verifying outputs..."

expected_outputs=(
    "../results/reduced/mutations/gene_level_results.csv"
    "../results/reduced/mutations/truncation_level_results.csv"
    "../results/full/mutations/gene_level_results.csv"
    "../results/full/mutations/truncation_level_results.csv"
)

all_outputs_present=true
for file in "${expected_outputs[@]}"; do
    if [ ! -f "$file" ]; then
        echo "✗ $(basename $file) missing"
        all_outputs_present=false
    else
        count=$(($(wc -l < "$file") - 1))
        echo "✓ $(basename $file) ($count rows)"
    fi
done

echo ""
echo "Checking visualization outputs..."
reduced_pdfs=$(find ../results/mutations/reduced -name "*.pdf" | wc -l)
full_pdfs=$(find ../results/mutations/full -name "*.pdf" | wc -l)

echo "✓ Reduced dataset visualizations: $reduced_pdfs PDFs"
echo "✓ Full dataset visualizations: $full_pdfs PDFs"

if [ "$all_outputs_present" = true ]; then
    echo ""
    echo "🎉 Mutation analysis complete!"
    echo "Results in ../results/mutations/"
    echo "Next step: sbatch 3_generate_proteins.sh"
else
    echo ""
    echo "❌ Some outputs are missing. Check logs and analyze_mutations.py."
    exit 1
fi
