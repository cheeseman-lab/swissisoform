#!/bin/bash

#SBATCH --job-name=summarize                  # Job name
#SBATCH --partition=20                        # Partition name
#SBATCH --ntasks=1                            # Run a single task
#SBATCH --cpus-per-task=4                     # CPUs per task
#SBATCH --mem=16G                             # Memory for the job
#SBATCH --time=2:00:00                        # Time limit
#SBATCH --output=out/summarize-%j.out         # Standard output log

# 5_summarize_results.sh
# Summarizes pipeline results: mutations and localization predictions

echo "=========================================================="
echo "SwissIsoform Pipeline Step 5: Summarize Results"
echo "=========================================================="

# Check if required input files exist
echo "Checking for required input files..."

DATASETS=("reduced" "full")
required_files=()
missing_files=()

for dataset in "${DATASETS[@]}"; do
    # Mutation analysis results
    mutation_gene_file="../results/$dataset/mutations/gene_level_results.csv"
    mutation_pair_file="../results/$dataset/mutations/truncation_level_results.csv"
    
    # Protein sequences
    proteins_file="../results/$dataset/proteins/protein_sequences_pairs.csv"
    proteins_mut_file="../results/$dataset/proteins/protein_sequences_with_mutations.csv"
    
    # Localization predictions (check for both pairs and mutations)
    loc_pairs_accurate="../results/$dataset/localization/protein_sequences_pairs_Accurate_results.csv"
    loc_pairs_fast="../results/$dataset/localization/protein_sequences_pairs_Fast_results.csv"
    loc_mut_accurate="../results/$dataset/localization/protein_sequences_mutations_Accurate_results.csv"
    loc_mut_fast="../results/$dataset/localization/protein_sequences_mutations_Fast_results.csv"
    
    # Check which files exist
    files_to_check=(
        "$mutation_gene_file"
        "$mutation_pair_file"
        "$proteins_file"
        "$proteins_mut_file"
        "$loc_pairs_accurate"
        "$loc_pairs_fast"
        "$loc_mut_accurate"
        "$loc_mut_fast"
    )
    
    echo ""
    echo "Checking $dataset dataset files..."
    for file in "${files_to_check[@]}"; do
        if [ -f "$file" ]; then
            echo "✓ $(basename $file)"
            required_files+=("$file")
        else
            echo "✗ $(basename $file) missing"
            missing_files+=("$file")
        fi
    done
done

# Check if we have the minimum required files
mutation_results_exist=false
localization_results_exist=false

for dataset in "${DATASETS[@]}"; do
    if [ -f "../results/$dataset/mutations/gene_level_results.csv" ]; then
        mutation_results_exist=true
    fi
    
    if [ -f "../results/$dataset/localization/protein_sequences_pairs_Accurate_results.csv" ] || \
       [ -f "../results/$dataset/localization/protein_sequences_pairs_Fast_results.csv" ]; then
        localization_results_exist=true
    fi
done

if [ "$mutation_results_exist" = false ]; then
    echo ""
    echo "❌ No mutation analysis results found!"
    echo "Please run: sbatch 2_analyze_mutations.sh"
    exit 1
fi

if [ "$localization_results_exist" = false ]; then
    echo ""
    echo "❌ No localization prediction results found!"
    echo "Please run: sbatch 4_predict_localization.sh"
    exit 1
fi

# Activate conda environment
echo ""
echo "Activating conda environment..."
source ~/.bashrc
conda activate swissisoform || {
    echo "❌ Failed to activate swissisoform conda environment"
    exit 1
}

# Create summary output directory
mkdir -p ../results/summary

echo ""
echo "Starting summary analysis..."
python3 summarize_results.py

# Verify outputs
echo ""
echo "Verifying summary outputs..."

expected_summary_files=(
    "../results/summary/mutation_summary.txt"
    "../results/summary/localization_summary.txt"
    "../results/summary/genes_with_localization_changes.csv"
    "../results/summary/detailed_localization_analysis.csv"
)

all_summary_files_present=true
for file in "${expected_summary_files[@]}"; do
    if [ -f "$file" ]; then
        if [[ "$file" == *.csv ]]; then
            count=$(($(wc -l < "$file") - 1))  # Subtract header
            echo "✓ $(basename $file) ($count rows)"
        else
            echo "✓ $(basename $file)"
        fi
    else
        echo "✗ $(basename $file) missing"
        all_summary_files_present=false
    fi
done

if [ "$all_summary_files_present" = true ]; then
    echo ""
    echo "🎉 Results summary completed successfully!"
    echo ""
    echo "Generated summary files:"
    echo "  ├─ summary/mutation_summary.txt           # Overview of mutation analysis"
    echo "  ├─ summary/localization_summary.txt       # Overview of localization predictions"
    echo "  ├─ summary/genes_with_localization_changes.csv  # Genes with interesting localization changes"
    echo "  └─ summary/detailed_localization_analysis.csv   # Detailed localization comparison"
    echo ""
    echo "Key findings:"
    echo ""
    
    # Show key findings from the text summaries
    if [ -f "../results/summary/mutation_summary.txt" ]; then
        echo "=== MUTATION ANALYSIS SUMMARY ==="
        cat "../results/summary/mutation_summary.txt"
        echo ""
    fi
    
    if [ -f "../results/summary/localization_summary.txt" ]; then
        echo "=== LOCALIZATION ANALYSIS SUMMARY ==="
        cat "../results/summary/localization_summary.txt"
        echo ""
    fi
    
    # Show preview of genes with localization changes
    if [ -f "../results/summary/genes_with_localization_changes.csv" ]; then
        echo "=== GENES WITH LOCALIZATION CHANGES (Preview) ==="
        head -n 10 "../results/summary/genes_with_localization_changes.csv"
        total_genes=$(tail -n +2 "../results/summary/genes_with_localization_changes.csv" | wc -l)
        echo "... (showing first 9 genes, $total_genes total genes with localization changes)"
        echo ""
    fi
    
    echo "Pipeline summary completed!"
    echo "All detailed results are available in ../results/summary/"
    
else
    echo ""
    echo "❌ Summary generation failed. Some output files are missing."
    exit 1
fi