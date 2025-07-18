#!/bin/bash

#SBATCH --job-name=proteins                # Job name
#SBATCH --partition=20                     # Partition name
#SBATCH --array=1-4                        # Job array with 4 tasks
#SBATCH --cpus-per-task=6                  # CPUs per task
#SBATCH --mem=16G                          # Memory per task (64G total = 16G × 4)
#SBATCH --time=24:00:00                    # Time limit (hrs:min:sec)
#SBATCH --output=out/proteins-%A_%a.out    # %A = job ID, %a = array task ID

# 3_generate_proteins.sh
# Generates protein sequence datasets from truncated transcripts

echo "======================================================="
echo "SwissIsoform Pipeline Step 3: Generate Protein Datasets"
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
echo "======================================================="

# Only run setup checks on the first task to avoid race conditions
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
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

    # Create results directory structure
    mkdir -p ../results/reduced/proteins
    mkdir -p ../results/full/proteins
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

# Run the appropriate task based on array task ID
case $SLURM_ARRAY_TASK_ID in
    1)
        # Task 1: Reduced dataset pairs
        echo "Array Task 1: Starting reduced pairs generation at $(date)"
        python3 generate_proteins.py "../data/ribosome_profiling/gene_list_reduced.txt" "../results/reduced/proteins" \
          --genome "$GENOME_PATH" \
          --annotation "$ANNOTATION_PATH" \
          --bed "$TRUNCATIONS_PATH" \
          --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
          --min-length "$MIN_LENGTH" \
          --max-length "$MAX_LENGTH" \
          --format "$FORMAT"
        echo "Array Task 1: Completed reduced pairs generation at $(date)"
        ;;
    2)
        # Task 2: Reduced dataset with mutations
        echo "Array Task 2: Starting reduced mutations generation at $(date)"
        python3 generate_proteins.py "../data/ribosome_profiling/gene_list_reduced.txt" "../results/reduced/proteins" \
          --genome "$GENOME_PATH" \
          --annotation "$ANNOTATION_PATH" \
          --bed "$TRUNCATIONS_PATH" \
          --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
          --min-length "$MIN_LENGTH" \
          --max-length "$MAX_LENGTH" \
          --format "$FORMAT" \
          --mutations \
          --impact-types "missense variant"
        echo "Array Task 2: Completed reduced mutations generation at $(date)"
        ;;
    3)
        # Task 3: Full dataset pairs
        echo "Array Task 3: Starting full pairs generation at $(date)"
        python3 generate_proteins.py "../data/ribosome_profiling/gene_list.txt" "../results/full/proteins" \
          --genome "$GENOME_PATH" \
          --annotation "$ANNOTATION_PATH" \
          --bed "$TRUNCATIONS_PATH" \
          --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
          --min-length "$MIN_LENGTH" \
          --max-length "$MAX_LENGTH" \
          --format "$FORMAT"
        echo "Array Task 3: Completed full pairs generation at $(date)"
        ;;
    4)
        # Task 4: Full dataset with mutations
        echo "Array Task 4: Starting full mutations generation at $(date)"
        python3 generate_proteins.py "../data/ribosome_profiling/gene_list.txt" "../results/full/proteins" \
          --genome "$GENOME_PATH" \
          --annotation "$ANNOTATION_PATH" \
          --bed "$TRUNCATIONS_PATH" \
          --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
          --min-length "$MIN_LENGTH" \
          --max-length "$MAX_LENGTH" \
          --format "$FORMAT" \
          --mutations \
          --impact-types "missense variant"
        echo "Array Task 4: Completed full mutations generation at $(date)"
        ;;
    *)
        echo "Unknown array task ID: $SLURM_ARRAY_TASK_ID"
        exit 1
        ;;
esac

# Only verify outputs on the last task to complete
if [ "$SLURM_ARRAY_TASK_ID" -eq 4 ]; then
    # Wait a moment for any file system sync
    sleep 5
    
    echo ""
    echo "Verifying generated datasets..."

    expected_files=(
        "../results/reduced/proteins/protein_sequences_pairs.fasta"
        "../results/reduced/proteins/protein_sequences_pairs.csv"
        "../results/reduced/proteins/protein_sequences_with_mutations.fasta"
        "../results/reduced/proteins/protein_sequences_with_mutations.csv"
        "../results/full/proteins/protein_sequences_pairs.fasta"
        "../results/full/proteins/protein_sequences_pairs.csv"
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
        echo "  │  ├─ protein_sequences_pairs.*       # Canonical + truncated pairs"
        echo "  │  └─ protein_sequences_with_mutations.*  # With mutations"
        echo "  └─ full/proteins/               # All truncation sites"
        echo "     ├─ protein_sequences_pairs.*       # Canonical + truncated pairs"
        echo "     └─ protein_sequences_with_mutations.*  # With mutations"
        echo ""
        echo "Next step:"
        echo "  Run: sbatch 4_predict_localization.sh"
    else
        echo ""
        echo "❌ Protein generation failed. Some output files are missing."
        exit 1
    fi
fi