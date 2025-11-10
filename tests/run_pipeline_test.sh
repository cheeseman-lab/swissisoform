#!/bin/bash
#
# SwissIsoform Pipeline Test Suite
#
# This script runs the complete SwissIsoform pipeline on a test gene list,
# covering both major aspects:
#   1. Mutation analysis (analyze_mutations.py)
#   2. Protein generation (generate_proteins.py)
#
# Usage:
#   bash tests/run_pipeline_test.sh [--debug] [--visualize]
#
# Options:
#   --debug      Enable debug logging for detailed output
#   --visualize  Generate visualizations for each gene
#   --clean      Remove previous test results before running
#

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default settings
VERBOSITY=0  # 0=WARNING, 1=INFO, 2=DEBUG
VISUALIZE=false
CLEAN=false
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -v)
            VERBOSITY=1
            shift
            ;;
        -vv)
            VERBOSITY=2
            shift
            ;;
        --debug)
            # Legacy flag for backwards compatibility
            VERBOSITY=2
            shift
            ;;
        --visualize)
            VISUALIZE=true
            shift
            ;;
        --clean)
            CLEAN=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [-v|-vv] [--visualize] [--clean]"
            echo ""
            echo "Options:"
            echo "  -v           Enable info logging (shows INFO messages)"
            echo "  -vv          Enable debug logging (shows DEBUG messages)"
            echo "  --debug      Same as -vv (legacy flag)"
            echo "  --visualize  Generate visualizations"
            echo "  --clean      Remove previous test results"
            echo "  -h, --help   Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
GENE_LIST="$SCRIPT_DIR/test_genes.txt"
OUTPUT_DIR="$SCRIPT_DIR/test_run_$TIMESTAMP"
GENOME_DATA="$PROJECT_ROOT/data/genome_data"

# Use default hela dataset
DATASET="hela"
BED_FILE="$PROJECT_ROOT/data/ribosome_profiling/${DATASET}_isoforms_with_transcripts.bed"

# Configuration
GENOME_FA="$GENOME_DATA/GRCh38.p7.genome.fa"

# Use v47names GTF if available (hela uses gencode v25)
ANNOTATION_GTF_V47="$GENOME_DATA/gencode.v25.annotation.v47names.gtf"
ANNOTATION_GTF_ORIG="$GENOME_DATA/gencode.v25.annotation.gtf"

if [ -f "$ANNOTATION_GTF_V47" ]; then
    ANNOTATION_GTF="$ANNOTATION_GTF_V47"
    echo -e "${GREEN}→${NC} Using v47-updated GTF annotation"
else
    ANNOTATION_GTF="$ANNOTATION_GTF_ORIG"
    echo -e "${YELLOW}→${NC} Using original GTF annotation (v47 update not found)"
fi

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║           SwissIsoform Pipeline Test Suite                   ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Validation
echo -e "${YELLOW}→${NC} Validating setup..."

if [ ! -f "$GENE_LIST" ]; then
    echo -e "${RED}✗ Gene list not found: $GENE_LIST${NC}"
    exit 1
fi

# Count non-comment, non-empty lines
GENE_COUNT=$(grep -v '^#' "$GENE_LIST" | grep -v '^$' | wc -l)
if [ "$GENE_COUNT" -eq 0 ]; then
    echo -e "${RED}✗ No genes found in $GENE_LIST${NC}"
    echo -e "${YELLOW}  Add gene names (one per line) to the test gene list.${NC}"
    exit 1
fi

if [ ! -f "$GENOME_FA" ]; then
    echo -e "${RED}✗ Genome FASTA not found: $GENOME_FA${NC}"
    exit 1
fi

if [ ! -f "$ANNOTATION_GTF" ]; then
    echo -e "${RED}✗ Annotation GTF not found: $ANNOTATION_GTF${NC}"
    exit 1
fi

if [ ! -f "$BED_FILE" ]; then
    echo -e "${RED}✗ BED file not found: $BED_FILE${NC}"
    exit 1
fi

echo -e "${GREEN}✓${NC} All required files found"
echo ""

# Display configuration
echo -e "${BLUE}Configuration:${NC}"
echo "  Dataset:      $DATASET"
echo "  Gene list:    $GENE_LIST"
echo "  Gene count:   $GENE_COUNT"
echo "  Output dir:   $OUTPUT_DIR"
echo "  Genome:       $GENOME_FA"
echo "  Annotation:   $ANNOTATION_GTF"
echo "  BED file:     $BED_FILE"
echo "  Log level:    $([ $VERBOSITY -eq 0 ] && echo 'WARNING' || [ $VERBOSITY -eq 1 ] && echo 'INFO' || echo 'DEBUG')"
echo "  Visualize:    $VISUALIZE"
echo ""

# Clean previous results if requested
if [ "$CLEAN" = true ]; then
    echo -e "${YELLOW}→${NC} Cleaning previous test results..."
    rm -rf "$SCRIPT_DIR/test_run_"*
    echo -e "${GREEN}✓${NC} Previous results removed"
    echo ""
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Set up log file
LOG_FILE="$OUTPUT_DIR/pipeline_test.log"

# Function to log to both console and file
# Strips ANSI color codes when writing to log file
log_output() {
    # Use tee to write to both stdout and log file
    # sed removes ANSI escape codes from the log file copy
    tee >(sed 's/\x1b\[[0-9;]*m//g' >> "$LOG_FILE")
}

# Start logging
echo "Pipeline Test Started: $(date)" | log_output
echo "" | log_output

# Build command arguments
MUTATION_ARGS="$GENE_LIST $OUTPUT_DIR/mutations"
MUTATION_ARGS="$MUTATION_ARGS --genome $GENOME_FA"
MUTATION_ARGS="$MUTATION_ARGS --annotation $ANNOTATION_GTF"
MUTATION_ARGS="$MUTATION_ARGS --bed $BED_FILE"
MUTATION_ARGS="$MUTATION_ARGS --sources clinvar gnomad cosmic"

# Add verbosity flags
if [ $VERBOSITY -eq 1 ]; then
    MUTATION_ARGS="$MUTATION_ARGS -v"
elif [ $VERBOSITY -eq 2 ]; then
    MUTATION_ARGS="$MUTATION_ARGS -vv"
fi

if [ "$VISUALIZE" = true ]; then
    MUTATION_ARGS="$MUTATION_ARGS --visualize"
fi

PROTEIN_ARGS="$GENE_LIST $OUTPUT_DIR/proteins"
PROTEIN_ARGS="$PROTEIN_ARGS --genome $GENOME_FA"
PROTEIN_ARGS="$PROTEIN_ARGS --annotation $ANNOTATION_GTF"
PROTEIN_ARGS="$PROTEIN_ARGS --bed $BED_FILE"
PROTEIN_ARGS="$PROTEIN_ARGS --mutations-file $OUTPUT_DIR/mutations/isoform_level_results.csv"

# Add verbosity flags
if [ $VERBOSITY -eq 1 ]; then
    PROTEIN_ARGS="$PROTEIN_ARGS -v"
elif [ $VERBOSITY -eq 2 ]; then
    PROTEIN_ARGS="$PROTEIN_ARGS -vv"
fi

# Save test configuration
cat > "$OUTPUT_DIR/test_config.txt" << EOF
SwissIsoform Pipeline Test Run
==============================
Timestamp: $TIMESTAMP
Gene count: $GENE_COUNT
Log level: $([ $VERBOSITY -eq 0 ] && echo 'WARNING' || [ $VERBOSITY -eq 1 ] && echo 'INFO' || echo 'DEBUG')
Visualize: $VISUALIZE

Paths:
------
Gene list: $GENE_LIST
Genome: $GENOME_FA
Annotation: $ANNOTATION_GTF
BED file: $BED_FILE

Test Genes:
-----------
$(cat "$GENE_LIST" | grep -v '^#' | grep -v '^$')
EOF

echo -e "${GREEN}✓${NC} Configuration saved to $OUTPUT_DIR/test_config.txt"
echo ""

# ============================================================================
# STEP 1: Mutation Analysis
# ============================================================================

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  STEP 1: Mutation Analysis                                   ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${YELLOW}→${NC} Running analyze_mutations.py..."
echo ""

START_TIME=$(date +%s)

# Run analysis - use 'tee' to capture output while preserving TTY for progress bars
# Only redirect stdout to log, let stderr go to terminal for progress bars
if python -u "$PROJECT_ROOT/scripts/analyze_mutations.py" $MUTATION_ARGS 2> >(tee -a "$LOG_FILE" >&2); then
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo ""
    echo -e "${GREEN}✓ Mutation analysis completed in ${DURATION}s${NC}" | log_output

    # Check output files
    if [ -f "$OUTPUT_DIR/mutations/gene_level_results.csv" ]; then
        GENE_SUCCESS=$(tail -n +2 "$OUTPUT_DIR/mutations/gene_level_results.csv" | wc -l)
        echo -e "${GREEN}  → Generated results for $GENE_SUCCESS genes${NC}" | log_output
    fi

    if [ -f "$OUTPUT_DIR/mutations/isoform_level_results.csv" ]; then
        PAIR_COUNT=$(tail -n +2 "$OUTPUT_DIR/mutations/isoform_level_results.csv" | wc -l)
        echo -e "${GREEN}  → Generated $PAIR_COUNT transcript-isoform pairs${NC}" | log_output
    fi
else
    echo ""
    echo -e "${RED}✗ Mutation analysis failed${NC}" | log_output
    echo -e "${YELLOW}  Check logs above for errors${NC}" | log_output
    exit 1
fi

echo ""

# ============================================================================
# STEP 2: Protein Generation
# ============================================================================

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  STEP 2: Protein Generation                                  ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${YELLOW}→${NC} Running generate_proteins.py..."
echo ""

START_TIME=$(date +%s)

# Run protein generation - use 'tee' to capture output while preserving TTY for progress bars
# Only redirect stdout to log, let stderr go to terminal for progress bars
if python -u "$PROJECT_ROOT/scripts/generate_proteins.py" $PROTEIN_ARGS 2> >(tee -a "$LOG_FILE" >&2); then
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo ""
    echo -e "${GREEN}✓ Protein generation completed in ${DURATION}s${NC}" | log_output

    # Check output files
    if [ -f "$OUTPUT_DIR/proteins/protein_sequences_with_mutations.csv" ]; then
        PROTEIN_COUNT=$(tail -n +2 "$OUTPUT_DIR/proteins/protein_sequences_with_mutations.csv" | wc -l)
        echo -e "${GREEN}  → Generated $PROTEIN_COUNT protein sequences${NC}" | log_output
    fi

    if [ -f "$OUTPUT_DIR/proteins/protein_sequences_with_mutations.fasta" ]; then
        FASTA_COUNT=$(grep -c '^>' "$OUTPUT_DIR/proteins/protein_sequences_with_mutations.fasta" || true)
        echo -e "${GREEN}  → Generated FASTA with $FASTA_COUNT sequences${NC}" | log_output
    fi
else
    echo ""
    echo -e "${RED}✗ Protein generation failed${NC}" | log_output
    echo -e "${YELLOW}  Check logs above for errors${NC}" | log_output
    exit 1
fi

echo ""

# ============================================================================
# SUMMARY
# ============================================================================

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Test Summary                                                ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${GREEN}✓ Pipeline test completed successfully!${NC}"
echo ""
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Output structure:"
echo "  $OUTPUT_DIR/"
echo "  ├── test_config.txt              # Test configuration"
echo "  ├── mutations/"
echo "  │   ├── gene_level_results.csv   # Gene-level mutation summary"
echo "  │   └── isoform_level_results.csv  # Transcript-isoform pairs"
echo "  └── proteins/"
echo "      ├── protein_sequences_with_mutations.csv   # Protein sequences (CSV)"
echo "      └── protein_sequences_with_mutations.fasta # Protein sequences (FASTA)"
echo ""

# Create a summary file
cat > "$OUTPUT_DIR/SUMMARY.md" << EOF
# SwissIsoform Pipeline Test Summary

**Test Run:** $TIMESTAMP
**Genes Tested:** $GENE_COUNT

## Results

### Mutation Analysis
- Gene-level results: \`mutations/gene_level_results.csv\`
- Transcript-isoform pairs: \`mutations/isoform_level_results.csv\`

### Protein Generation
- Protein sequences (CSV): \`proteins/protein_sequences_with_mutations.csv\`
- Protein sequences (FASTA): \`proteins/protein_sequences_with_mutations.fasta\`

## Test Genes

\`\`\`
$(cat "$GENE_LIST" | grep -v '^#' | grep -v '^$')
\`\`\`

## Configuration

See \`test_config.txt\` for full configuration details.

## Next Steps

1. Review the results in the output files
2. Check for expected truncations/extensions
3. Validate protein sequences
4. Look for edge cases and failures

EOF

echo -e "${GREEN}✓ Summary saved to $OUTPUT_DIR/SUMMARY.md${NC}" | log_output
echo "" | log_output
echo -e "${BLUE}Next steps:${NC}" | log_output
echo "  1. Review results: cd $OUTPUT_DIR" | log_output
echo "  2. Check gene-level results: less mutations/gene_level_results.csv" | log_output
echo "  3. Check protein sequences: less proteins/protein_sequences_with_mutations.csv" | log_output
echo "" | log_output
echo -e "${GREEN}Test complete! 🎉${NC}" | log_output
echo "" | log_output
echo "Pipeline Test Completed: $(date)" | log_output
echo "Log saved to: $LOG_FILE" | log_output
