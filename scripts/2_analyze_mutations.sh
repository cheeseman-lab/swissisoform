#!/bin/bash

#SBATCH --job-name=mutations               # Job name
#SBATCH --partition=20                     # Partition name
#SBATCH --array=1-2                        # Job array with 2 tasks
#SBATCH --cpus-per-task=8                  # CPUs per task
#SBATCH --mem=32G                          # Memory per task (64G total = 32G × 2)
#SBATCH --time=48:00:00                    # Time limit (hrs:min:sec)
#SBATCH --output=out/mutations-%A_%a.out   # %A = job ID, %a = array task ID

# 2_analyze_mutations.sh
# Analyzes mutations in alternative isoform truncation regions

echo "========================================================"
echo "SwissIsoform Pipeline Step 2: Analyze Mutations"
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
echo "========================================================"

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
        "$RIBOPROF_DIR/isoforms_gene_list_reduced.txt"
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
    mkdir -p ../results/reduced/mutations
    mkdir -p ../results/full/mutations
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
TRUNCATIONS_PATH="../data/ribosome_profiling/isoforms_with_transcripts.bed"

# Run the appropriate task based on array task ID
case $SLURM_ARRAY_TASK_ID in
    1)
        # Task 1: Reduced dataset with visualization
        echo "Array Task 1: Starting reduced mutations analysis at $(date)"
        python3 analyze_mutations.py "../data/ribosome_profiling/isoforms_gene_list_reduced.txt" "../results/reduced/mutations" \
          --genome "$GENOME_PATH" \
          --annotation "$ANNOTATION_PATH" \
          --bed "$TRUNCATIONS_PATH" \
          --sources "clinvar" \
          --impact-types "missense variant" "nonsense variant" "frameshift variant" "synonymous variant" "inframe deletion" "inframe insertion" \
          --visualize
        echo "Array Task 1: Completed reduced mutations analysis at $(date)"
        ;;
    2)
        # Task 2: Full dataset
        echo "Array Task 2: Starting full mutations analysis at $(date)"
        python3 analyze_mutations.py "../data/ribosome_profiling/isoforms_gene_list.txt" "../results/full/mutations" \
          --genome "$GENOME_PATH" \
          --annotation "$ANNOTATION_PATH" \
          --bed "$TRUNCATIONS_PATH" \
          --sources "clinvar" \
          --impact-types "missense variant" "nonsense variant" "frameshift variant" "synonymous variant" "inframe deletion" "inframe insertion" \
        echo "Array Task 2: Completed full mutations analysis at $(date)"
        ;;
    *)
        echo "Unknown array task ID: $SLURM_ARRAY_TASK_ID"
        exit 1
        ;;
esac

# Only verify outputs on the last task to complete
if [ "$SLURM_ARRAY_TASK_ID" -eq 2 ]; then
    # Wait a moment for any file system sync
    sleep 5
    
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
fi