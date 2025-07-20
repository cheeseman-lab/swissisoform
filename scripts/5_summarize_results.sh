#!/bin/bash

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

# Run summary analysis for both datasets
echo ""
echo "Starting summary analysis..."
python3 summarize_results.py

echo ""
echo "🎉 Pipeline summary analysis completed!"
echo ""
echo "Generated summary files by dataset and model:"

for dataset in "${DATASETS[@]}"; do
    summary_dir="../results/$dataset/summary"
    if [ -d "$summary_dir" ]; then
        echo ""
        echo "$dataset dataset summary:"
        echo "  ├─ $dataset/summary/mutation_summary.txt (shared)"
        
        # Check for model-specific subdirectories
        for model in "accurate" "fast"; do
            model_dir="$summary_dir/$model"
            if [ -d "$model_dir" ]; then
                echo "  ├─ $dataset/summary/$model/ (${model^} model results)"
                echo "  │  ├─ localization_summary.txt"
                echo "  │  ├─ genes_with_localization_changes.csv"
                echo "  │  ├─ detailed_localization_analysis.csv"
                echo "  │  └─ gene_level_summary.csv"
            fi
        done
        
        # Show key findings for each dataset and model
        echo ""
        echo "=== $dataset DATASET KEY FINDINGS ==="
        
        if [ -f "$summary_dir/mutation_summary.txt" ]; then
            echo ""
            echo "MUTATION ANALYSIS (shared across models):"
            cat "$summary_dir/mutation_summary.txt"
        fi
        
        # Show findings for each model
        for model in "accurate" "fast"; do
            model_dir="$summary_dir/$model"
            if [ -d "$model_dir" ]; then
                echo ""
                echo "LOCALIZATION ANALYSIS - ${model^} MODEL:"
                if [ -f "$model_dir/localization_summary.txt" ]; then
                    cat "$model_dir/localization_summary.txt"
                fi
                
                # Show preview of genes with localization changes for this model
                if [ -f "$model_dir/genes_with_localization_changes.csv" ]; then
                    echo ""
                    echo "GENES WITH LOCALIZATION CHANGES - ${model^} MODEL (Preview):"
                    head -n 6 "$model_dir/genes_with_localization_changes.csv"
                    total_genes=$(tail -n +2 "$model_dir/genes_with_localization_changes.csv" | wc -l)
                    if [ $total_genes -gt 5 ]; then
                        echo "... (showing first 5 genes, $total_genes total genes with localization changes in ${model^} model)"
                    fi
                fi
            fi
        done
        
        echo ""
        echo "============================================================"
    fi
done

echo ""
echo "Pipeline summary completed!"
echo "Detailed results are available in ../results/[dataset]/summary/[model]/"