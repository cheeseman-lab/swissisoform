#!/bin/bash
#
# SwissIsoform Pipeline Step 5: Summarize Results
#
# This script summarizes pipeline results including mutation analysis and
# localization predictions. It generates comprehensive reports for all
# processed datasets.
#
# Usage:
#   bash 5_summarize_results.sh
#
# Prerequisites:
#   - 2_analyze_mutations.sh must have been run
#   - 4_predict_localization.sh must have been run
#   - Mutation and localization results must exist
#

set -e  # Exit on error

# Source shared utilities (colors, helper functions)
# Use relative path from scripts directory
UTILS_PATH="$(dirname "$0")/utils.sh"
if [ -f "$UTILS_PATH" ]; then
    source "$UTILS_PATH"
elif [ -f "utils.sh" ]; then
    source "utils.sh"
elif [ -f "scripts/utils.sh" ]; then
    source "scripts/utils.sh"
else
    # Fallback: define colors inline
    if [ -t 1 ]; then
        RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
        BLUE='\033[0;34m'; CYAN='\033[0;36m'; NC='\033[0m'
    else
        RED=''; GREEN=''; YELLOW=''; BLUE=''; CYAN=''; NC=''
    fi
fi

# Start timing
START_TIME=$(date +%s)

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   SwissIsoform Pipeline Step 5: Summarize Results            ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "This script analyzes pipeline results and generates comprehensive summaries."
echo ""

# ============================================================================
# Validation
# ============================================================================

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Validation                                                  ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

DATASETS=("reduced" "full")
required_files=()
missing_files=()

echo -e "${YELLOW}→${NC} Checking for input files..."
echo ""

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

    echo -e "${CYAN}$dataset dataset:${NC}"
    for file in "${files_to_check[@]}"; do
        if [ -f "$file" ]; then
            SIZE=$(du -h "$file" | cut -f1)
            echo -e "  ${GREEN}✓${NC} $(basename $file) (${SIZE})"
            required_files+=("$file")
        else
            echo -e "  ${YELLOW}⚠${NC} $(basename $file) missing"
            missing_files+=("$file")
        fi
    done
    echo ""
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
    echo -e "${RED}✗ No mutation analysis results found!${NC}"
    echo "Please run: sbatch 2_analyze_mutations.sh"
    exit 1
fi

if [ "$localization_results_exist" = false ]; then
    echo ""
    echo -e "${RED}✗ No localization prediction results found!${NC}"
    echo "Please run: sbatch 4_predict_localization.sh"
    exit 1
fi

echo -e "${GREEN}✓ Minimum required files found${NC}"

# ============================================================================
# Environment Setup
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Environment Setup                                           ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${YELLOW}→${NC} Activating conda environment..."
source ~/.bashrc
conda activate swissisoform || {
    echo -e "${RED}✗${NC} Failed to activate swissisoform conda environment"
    echo ""
    echo "Please run: conda env create --file=../environment.yml"
    exit 1
}
echo -e "${GREEN}✓${NC} Environment activated"

# ============================================================================
# Summary Generation
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Generating Summary                                          ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${YELLOW}→${NC} Analyzing pipeline results..."
echo ""
python3 summarize_results.py

# Calculate duration
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
if [ $DURATION -lt 60 ]; then
    DURATION_STR="${DURATION}s"
elif [ $DURATION -lt 3600 ]; then
    MINUTES=$((DURATION / 60))
    SECS=$((DURATION % 60))
    DURATION_STR="${MINUTES}m ${SECS}s"
else
    HOURS=$((DURATION / 3600))
    MINUTES=$(((DURATION % 3600) / 60))
    SECS=$((DURATION % 60))
    DURATION_STR="${HOURS}h ${MINUTES}m ${SECS}s"
fi

# ============================================================================
# Results Display
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Summary Complete                                            ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${GREEN}✓ Pipeline summary analysis completed in $DURATION_STR${NC}"
echo ""

# Display generated summary files
echo "Generated summary files by dataset and model:"
echo ""

for dataset in "${DATASETS[@]}"; do
    summary_dir="../results/$dataset/summary"
    if [ -d "$summary_dir" ]; then
        echo -e "${CYAN}$dataset dataset summary:${NC}"

        # Check for mutation summary
        if [ -f "$summary_dir/mutation_summary.txt" ]; then
            echo -e "  ${GREEN}✓${NC} mutation_summary.txt (shared across models)"
        fi

        # Check for model-specific subdirectories
        for model in "accurate" "fast"; do
            model_dir="$summary_dir/$model"
            if [ -d "$model_dir" ]; then
                echo -e "  ${GREEN}✓${NC} $model/ (${model^} model results)"

                # List files in model directory
                if [ -f "$model_dir/localization_summary.txt" ]; then
                    echo -e "    ├─ localization_summary.txt"
                fi
                if [ -f "$model_dir/genes_with_localization_changes.csv" ]; then
                    echo -e "    ├─ genes_with_localization_changes.csv"
                fi
                if [ -f "$model_dir/detailed_localization_analysis.csv" ]; then
                    echo -e "    ├─ detailed_localization_analysis.csv"
                fi
                if [ -f "$model_dir/gene_level_summary.csv" ]; then
                    echo -e "    └─ gene_level_summary.csv"
                fi
            fi
        done
        echo ""
    fi
done

# Display key findings
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Key Findings                                                ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

for dataset in "${DATASETS[@]}"; do
    summary_dir="../results/$dataset/summary"
    if [ -d "$summary_dir" ]; then
        echo -e "${CYAN}═══ $dataset DATASET ═══${NC}"
        echo ""

        # Show mutation analysis findings
        if [ -f "$summary_dir/mutation_summary.txt" ]; then
            echo -e "${YELLOW}Mutation Analysis:${NC}"
            cat "$summary_dir/mutation_summary.txt"
            echo ""
        fi

        # Show findings for each model
        for model in "accurate" "fast"; do
            model_dir="$summary_dir/$model"
            if [ -d "$model_dir" ]; then
                echo -e "${YELLOW}Localization Analysis - ${model^} Model:${NC}"

                if [ -f "$model_dir/localization_summary.txt" ]; then
                    cat "$model_dir/localization_summary.txt"
                fi

                # Show preview of genes with localization changes
                if [ -f "$model_dir/genes_with_localization_changes.csv" ]; then
                    echo ""
                    echo -e "${CYAN}Genes with Localization Changes (Preview):${NC}"
                    head -n 6 "$model_dir/genes_with_localization_changes.csv"

                    total_genes=$(tail -n +2 "$model_dir/genes_with_localization_changes.csv" | wc -l)
                    if [ $total_genes -gt 5 ]; then
                        echo "... (showing first 5 genes, $total_genes total with localization changes)"
                    fi
                fi
                echo ""
            fi
        done

        echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
        echo ""
    fi
done

# ============================================================================
# Final Summary
# ============================================================================

echo -e "${GREEN}✓ Pipeline summary completed successfully!${NC}"
echo ""
echo "Detailed results available in:"
for dataset in "${DATASETS[@]}"; do
    if [ -d "../results/$dataset/summary" ]; then
        echo "  └─ ../results/$dataset/summary/"
        echo "     ├─ mutation_summary.txt"
        for model in "accurate" "fast"; do
            if [ -d "../results/$dataset/summary/$model" ]; then
                echo "     └─ $model/"
            fi
        done
    fi
done
echo ""
