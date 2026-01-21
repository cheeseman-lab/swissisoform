#!/bin/bash
#
# SwissIsoform Pipeline Step 5: Summarize Results
#
# This script summarizes pipeline results including mutation analysis and
# localization predictions. It generates comprehensive reports for the
# specified dataset.
#
# Usage:
#   bash 5_summarize_results.sh
#   DATASET=hela bash 5_summarize_results.sh
#   DATASET=hela_bch bash 5_summarize_results.sh
#
# Environment Variables:
#   DATASET - Dataset to process (default: hela)
#
# Prerequisites:
#   - 2_analyze_mutations.sh must have been run for the dataset
#   - 4_predict_localization.sh must have been run for the dataset
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

# Dataset selection (default: hela)
DATASET="${DATASET:-hela}"

# Start timing
START_TIME=$(date +%s)

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   SwissIsoform Pipeline Step 5: Summarize Results            ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Dataset: $DATASET"
echo "This script analyzes pipeline results and generates comprehensive summaries."
echo ""

# ============================================================================
# Validation
# ============================================================================

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Validation                                                  ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

required_files=()
missing_files=()

echo -e "${YELLOW}→${NC} Checking for input files..."
echo ""

# Mutation analysis results
mutation_gene_file="../results/$DATASET/mutations/gene_level_results.csv"
mutation_isoform_file="../results/$DATASET/mutations/isoform_level_results.csv"

# Protein sequences
proteins_file="../results/$DATASET/proteins/protein_sequences_pairs.csv"
proteins_mut_file="../results/$DATASET/proteins/protein_sequences_with_mutations.csv"

# Localization predictions (check for both pairs and mutations)
loc_pairs_accurate="../results/$DATASET/localization/protein_sequences_pairs_Accurate_results.csv"
loc_pairs_fast="../results/$DATASET/localization/protein_sequences_pairs_Fast_results.csv"
loc_mut_accurate="../results/$DATASET/localization/protein_sequences_mutations_Accurate_results.csv"
loc_mut_fast="../results/$DATASET/localization/protein_sequences_mutations_Fast_results.csv"

# Check which files exist
files_to_check=(
    "$mutation_gene_file"
    "$mutation_isoform_file"
    "$proteins_file"
    "$proteins_mut_file"
    "$loc_pairs_accurate"
    "$loc_pairs_fast"
    "$loc_mut_accurate"
    "$loc_mut_fast"
)

echo -e "${CYAN}$DATASET dataset:${NC}"
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

# Check if we have the minimum required files
mutation_results_exist=false
localization_results_exist=false

if [ -f "../results/$DATASET/mutations/gene_level_results.csv" ] || \
   [ -f "../results/$DATASET/mutations/isoform_level_results.csv" ]; then
    mutation_results_exist=true
fi

if [ -f "../results/$DATASET/localization/protein_sequences_pairs_Accurate_results.csv" ] || \
   [ -f "../results/$DATASET/localization/protein_sequences_pairs_Fast_results.csv" ]; then
    localization_results_exist=true
fi

if [ "$mutation_results_exist" = false ]; then
    echo ""
    echo -e "${RED}✗ No mutation analysis results found for $DATASET!${NC}"
    echo "Please run: sbatch --export=DATASET=$DATASET 2_analyze_mutations.sh"
    exit 1
fi

if [ "$localization_results_exist" = false ]; then
    echo ""
    echo -e "${RED}✗ No localization prediction results found for $DATASET!${NC}"
    echo "Please run: sbatch --export=DATASET=$DATASET 4_predict_localization.sh"
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

echo -e "${YELLOW}→${NC} Analyzing pipeline results and generating summaries for $DATASET..."
echo ""
python3 summarize_results.py "$DATASET"

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
echo "Generated summary files for $DATASET:"
echo ""

summary_dir="../results/$DATASET/summary"
if [ -d "$summary_dir" ]; then
    echo -e "${CYAN}$DATASET dataset summary:${NC}"

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
                echo -e "    ├─ gene_level_summary.csv"
            fi
            # Ranked summary TSVs for skill agent analysis
            if [ -f "$model_dir/summary_localization_changes.csv" ]; then
                echo -e "    ├─ summary_localization_changes.csv (ranked)"
            fi
            if [ -f "$model_dir/summary_deleterious_burden.csv" ]; then
                echo -e "    └─ summary_deleterious_burden.csv (ranked)"
            fi
        fi
    done
    echo ""
fi

# Display key findings
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Key Findings                                                ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

summary_dir="../results/$DATASET/summary"
if [ -d "$summary_dir" ]; then
    echo -e "${CYAN}═══ $DATASET DATASET ═══${NC}"
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

    # Show ranked summaries preview
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Ranked Summaries (for Skill Agent Analysis)                 ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""

    for model in "accurate" "fast"; do
        model_dir="$summary_dir/$model"
        if [ -d "$model_dir" ]; then
            # Show top localization changes
            if [ -f "$model_dir/summary_localization_changes.csv" ]; then
                echo -e "${YELLOW}Top Localization Changes - ${model^} Model:${NC}"
                head -n 6 "$model_dir/summary_localization_changes.csv" | column -t -s ','
                total_changes=$(tail -n +2 "$model_dir/summary_localization_changes.csv" | wc -l)
                echo "... ($total_changes total features)"
                echo ""
            fi

            # Show top deleterious burden
            if [ -f "$model_dir/summary_deleterious_burden.csv" ]; then
                echo -e "${YELLOW}Top Deleterious Burden - ${model^} Model:${NC}"
                head -n 6 "$model_dir/summary_deleterious_burden.csv" | cut -d',' -f1,5,12,13,14 | column -t -s ','
                total_burden=$(tail -n +2 "$model_dir/summary_deleterious_burden.csv" | wc -l)
                echo "... ($total_burden total features with mutations)"
                echo ""
            fi
        fi
    done

    echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
    echo ""
fi

# ============================================================================
# Final Summary
# ============================================================================

echo -e "${GREEN}✓ Pipeline summary completed successfully for $DATASET!${NC}"
echo ""
echo "Detailed results available in:"
if [ -d "../results/$DATASET/summary" ]; then
    echo "  └─ ../results/$DATASET/summary/"
    echo "     ├─ mutation_summary.txt"
    for model in "accurate" "fast"; do
        if [ -d "../results/$DATASET/summary/$model" ]; then
            echo "     └─ $model/"
        fi
    done
fi
echo ""
