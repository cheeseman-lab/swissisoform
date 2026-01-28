#!/bin/bash
#
# SwissIsoform Pipeline Step 5: Summarize Results
#
# This script summarizes pipeline results including mutation analysis and
# localization predictions. It generates comprehensive reports for the
# specified dataset and mutation source(s).
#
# Usage:
#   DATASET=hela SOURCE=gnomad bash 5_summarize_results.sh
#   DATASET=hela SOURCES="gnomad|clinvar|cosmic" bash 5_summarize_results.sh
#
# Environment Variables:
#   DATASET - Dataset to process (default: hela)
#   SOURCE  - Single mutation source (e.g., gnomad, clinvar, cosmic, custom_*)
#   SOURCES - Multiple sources pipe-separated (e.g., "gnomad|clinvar|cosmic")
#
# Prerequisites:
#   - 2_analyze_mutations.sh must have been run for the dataset+source
#   - 4_predict_localization.sh must have been run for the dataset+source
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

# Source selection
# - If SOURCES is set (pipe-separated), process multiple sources
# - Otherwise use SOURCE (default: gnomad)
if [ -n "$SOURCES" ]; then
    IFS='|' read -ra SOURCES_ARRAY <<< "$SOURCES"
    echo "Multi-source mode: ${SOURCES_ARRAY[@]}"
else
    SOURCE="${SOURCE:-gnomad}"
    SOURCES_ARRAY=("$SOURCE")
    echo "Single source mode: $SOURCE"
fi

# Start timing
START_TIME=$(date +%s)

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   SwissIsoform Pipeline Step 5: Summarize Results            ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Dataset: $DATASET"
echo "Sources: ${SOURCES_ARRAY[@]}"
echo "This script analyzes pipeline results and generates comprehensive summaries."
echo ""

# ============================================================================
# Validation
# ============================================================================

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Validation                                                  ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${YELLOW}→${NC} Checking for input files..."
echo ""

# Check files for each source
ALL_SOURCES_VALID=true

for SOURCE in "${SOURCES_ARRAY[@]}"; do
    echo -e "${CYAN}$DATASET / $SOURCE:${NC}"

    # Check if this source has required files
    SOURCE_VALID=true

    # Special handling for "default" source (structural analysis, no mutations)
    if [ "$SOURCE" = "default" ]; then
        # Default source only needs pairs localization files
        loc_pairs_accurate="../results/$DATASET/default/localization/protein_sequences_pairs_Accurate_results.csv"
        loc_pairs_fast="../results/$DATASET/default/localization/protein_sequences_pairs_Fast_results.csv"

        LOC_EXISTS=false
        if [ -f "$loc_pairs_accurate" ]; then
            SIZE=$(du -h "$loc_pairs_accurate" | cut -f1)
            echo -e "  ${GREEN}✓${NC} protein_sequences_pairs_Accurate_results.csv (${SIZE})"
            LOC_EXISTS=true
        else
            echo -e "  ${YELLOW}⚠${NC} Accurate pairs localization missing"
        fi

        if [ -f "$loc_pairs_fast" ]; then
            SIZE=$(du -h "$loc_pairs_fast" | cut -f1)
            echo -e "  ${GREEN}✓${NC} protein_sequences_pairs_Fast_results.csv (${SIZE})"
            LOC_EXISTS=true
        else
            echo -e "  ${YELLOW}⚠${NC} Fast pairs localization missing"
        fi

        if [ "$LOC_EXISTS" = false ]; then
            echo -e "  ${RED}✗${NC} No pairs localization results found"
            SOURCE_VALID=false
        fi
    else
        # Mutation sources need mutation files + localization
        # Mutation analysis results
        mutation_isoform_file="../results/$DATASET/$SOURCE/mutations/isoform_level_results.csv"

        # Protein sequences
        proteins_mut_file="../results/$DATASET/$SOURCE/proteins/protein_sequences_with_mutations.csv"

        # Localization predictions
        loc_mut_accurate="../results/$DATASET/$SOURCE/localization/protein_sequences_mutations_Accurate_results.csv"
        loc_mut_fast="../results/$DATASET/$SOURCE/localization/protein_sequences_mutations_Fast_results.csv"

        # Check mutation results
        if [ -f "$mutation_isoform_file" ]; then
            SIZE=$(du -h "$mutation_isoform_file" | cut -f1)
            echo -e "  ${GREEN}✓${NC} isoform_level_results.csv (${SIZE})"
        else
            echo -e "  ${RED}✗${NC} isoform_level_results.csv missing"
            SOURCE_VALID=false
        fi

        # Check protein results
        if [ -f "$proteins_mut_file" ]; then
            SIZE=$(du -h "$proteins_mut_file" | cut -f1)
            echo -e "  ${GREEN}✓${NC} protein_sequences_with_mutations.csv (${SIZE})"
        else
            echo -e "  ${YELLOW}⚠${NC} protein_sequences_with_mutations.csv missing (optional)"
        fi

        # Check localization results (at least one mode required)
        LOC_EXISTS=false
        if [ -f "$loc_mut_accurate" ]; then
            SIZE=$(du -h "$loc_mut_accurate" | cut -f1)
            echo -e "  ${GREEN}✓${NC} protein_sequences_mutations_Accurate_results.csv (${SIZE})"
            LOC_EXISTS=true
        else
            echo -e "  ${YELLOW}⚠${NC} Accurate localization results missing"
        fi

        if [ -f "$loc_mut_fast" ]; then
            SIZE=$(du -h "$loc_mut_fast" | cut -f1)
            echo -e "  ${GREEN}✓${NC} protein_sequences_mutations_Fast_results.csv (${SIZE})"
            LOC_EXISTS=true
        else
            echo -e "  ${YELLOW}⚠${NC} Fast localization results missing"
        fi

        if [ "$LOC_EXISTS" = false ]; then
            echo -e "  ${RED}✗${NC} No localization results found"
            SOURCE_VALID=false
        fi
    fi

    if [ "$SOURCE_VALID" = false ]; then
        ALL_SOURCES_VALID=false
        echo -e "  ${RED}→ Source $SOURCE validation failed${NC}"
    fi

    echo ""
done

# Also check for base (default) localization results
echo -e "${CYAN}$DATASET / default (base):${NC}"
BASE_PAIRS_EXIST=false

base_loc_accurate="../results/$DATASET/default/localization/protein_sequences_pairs_Accurate_results.csv"
base_loc_fast="../results/$DATASET/default/localization/protein_sequences_pairs_Fast_results.csv"

if [ -f "$base_loc_accurate" ]; then
    SIZE=$(du -h "$base_loc_accurate" | cut -f1)
    echo -e "  ${GREEN}✓${NC} protein_sequences_pairs_Accurate_results.csv (${SIZE})"
    BASE_PAIRS_EXIST=true
fi

if [ -f "$base_loc_fast" ]; then
    SIZE=$(du -h "$base_loc_fast" | cut -f1)
    echo -e "  ${GREEN}✓${NC} protein_sequences_pairs_Fast_results.csv (${SIZE})"
    BASE_PAIRS_EXIST=true
fi

if [ "$BASE_PAIRS_EXIST" = false ]; then
    echo -e "  ${YELLOW}⚠${NC} Base pair localization results missing (optional for summary)"
fi
echo ""

if [ "$ALL_SOURCES_VALID" = false ]; then
    echo -e "${RED}✗ Some sources are missing required files!${NC}"
    echo ""
    echo "Please ensure you have run:"
    echo "  1. sbatch --export=DATASET=$DATASET,SOURCES=\"${SOURCES}\" 2_analyze_mutations.sh"
    echo "  2. sbatch --export=DATASET=$DATASET,SOURCES=\"${SOURCES}\" 3_generate_proteins.sh"
    echo "  3. sbatch --export=DATASET=$DATASET,SOURCES=\"${SOURCES}\" 4_predict_localization.sh"
    exit 1
fi

echo -e "${GREEN}✓ All sources have minimum required files${NC}"

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

# Process each source
for SOURCE in "${SOURCES_ARRAY[@]}"; do
    echo -e "${YELLOW}→${NC} Analyzing pipeline results for $DATASET / $SOURCE..."
    echo ""

    # Call summarize_results.py with both dataset and source
    python3 summarize_results.py "$DATASET" "$SOURCE"

    echo ""
done

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

# Display summary for each source
for SOURCE in "${SOURCES_ARRAY[@]}"; do
    summary_dir="../results/$DATASET/$SOURCE/summary"
    if [ -d "$summary_dir" ]; then
        echo -e "${CYAN}═══ $DATASET / $SOURCE ═══${NC}"
        echo ""

        # Show the unified summary text
        if [ -f "$summary_dir/summary.txt" ]; then
            cat "$summary_dir/summary.txt"
            echo ""
        fi

        # List generated files
        echo -e "${YELLOW}Generated files:${NC}"
        for model in "accurate" "fast"; do
            model_dir="$summary_dir/$model"
            if [ -d "$model_dir" ] && [ -f "$model_dir/feature_analysis.csv" ]; then
                rows=$(tail -n +2 "$model_dir/feature_analysis.csv" | wc -l)
                echo -e "  ${GREEN}✓${NC} $model/feature_analysis.csv ($rows features)"
            fi
        done
        if [ -d "$summary_dir/model_comparison" ]; then
            echo -e "  ${GREEN}✓${NC} model_comparison/"
        fi
        echo ""

        echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
        echo ""
    fi
done

# ============================================================================
# Final Summary
# ============================================================================

echo -e "${GREEN}✓ Pipeline summary completed successfully for $DATASET!${NC}"
echo ""
echo "Detailed results available in:"

for SOURCE in "${SOURCES_ARRAY[@]}"; do
    if [ -d "../results/$DATASET/$SOURCE/summary" ]; then
        echo "  └─ ../results/$DATASET/$SOURCE/summary/"
        echo "     ├─ summary.txt"
        for model in "accurate" "fast"; do
            if [ -d "../results/$DATASET/$SOURCE/summary/$model" ]; then
                echo "     ├─ $model/feature_analysis.csv"
            fi
        done
        if [ -d "../results/$DATASET/$SOURCE/summary/model_comparison" ]; then
            echo "     └─ model_comparison/"
        fi
        echo ""
    fi
done
