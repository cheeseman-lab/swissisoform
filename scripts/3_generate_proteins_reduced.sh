#!/bin/bash

#SBATCH --job-name=proteins_reduced        # Job name
#SBATCH --partition=20                     # Partition name
#SBATCH --cpus-per-task=8                  # CPUs per task
#SBATCH --mem=32G                          # Memory per task
#SBATCH --time=12:00:00                    # Time limit (reduced since using fast mode)
#SBATCH --output=out/proteins_reduced-%j.out

# 3_generate_proteins_reduced.sh
# Fast protein sequence generation using pre-validated mutations from step 2

echo "======================================================="
echo "SwissIsoform Pipeline Step 3: Generate Proteins (Reduced - Fast Mode)"
echo "======================================================="

# Check if required input files exist
echo "Checking for required input files..."

GENOME_DIR="../data/genome_data"
RIBOPROF_DIR="../data/ribosome_profiling"
MUTATIONS_DIR="../results/reduced/mutations"

required_files=(
    "$GENOME_DIR/GRCh38.p7.genome.fa"
    "$GENOME_DIR/gencode.v25.annotation.ensembl_cleaned.gtf"
    "$RIBOPROF_DIR/isoforms_with_transcripts.bed"
    "$RIBOPROF_DIR/isoforms_gene_list_reduced.txt"
    "$MUTATIONS_DIR/isoform_level_results.csv"
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
    echo "❌ Missing required files! Run 1_cleanup_files.sh and 2_analyze_mutations_reduced.sh first"
    exit 1
fi

# Check that we have mutation results
mutation_count=$(tail -n +2 "$MUTATIONS_DIR/isoform_level_results.csv" | wc -l)
if [ "$mutation_count" -eq 0 ]; then
    echo "❌ No mutation results found! Run 2_analyze_mutations_reduced.sh first"
    exit 1
fi

echo "✓ Found $mutation_count validated mutation results"

# Create results directory structure
mkdir -p ../results/reduced/proteins

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
MUTATIONS_FILE="../results/reduced/mutations/isoform_level_results.csv"
GENE_LIST="../data/ribosome_profiling/isoforms_gene_list_reduced.txt"
OUTPUT_DIR="../results/reduced/proteins"
MIN_LENGTH=10
MAX_LENGTH=100000
FORMAT="fasta,csv"

echo ""
echo "Starting protein sequence generation..."
echo "Configuration:"
echo "  ├─ Gene list: $(basename $GENE_LIST)"
echo "  ├─ Pre-validated mutations: $(basename $MUTATIONS_FILE)"
echo "  ├─ Source: ClinVar"
echo "  ├─ Impact types: missense variant"
echo "  ├─ Length range: $MIN_LENGTH-$MAX_LENGTH amino acids"
echo "  └─ Output format: $FORMAT"

# Run the single script that generates both datasets
echo ""
echo "Generating both pairs and mutations datasets..."
python3 generate_proteins.py "$GENE_LIST" "$OUTPUT_DIR" \
  --mutations-file "$MUTATIONS_FILE" \
  --genome "$GENOME_PATH" \
  --annotation "$ANNOTATION_PATH" \
  --bed "$TRUNCATIONS_PATH" \
  --sources clinvar \
  --impact-types "missense variant" \
  --min-length "$MIN_LENGTH" \
  --max-length "$MAX_LENGTH" \
  --format "$FORMAT" \
  --fast-mode

exit_code=$?

if [ $exit_code -ne 0 ]; then
    echo "❌ Protein generation failed with exit code $exit_code"
    exit 1
fi

# Verify outputs
echo ""
echo "Verifying generated datasets..."

expected_files=(
    "../results/reduced/proteins/protein_sequences_pairs.fasta"
    "../results/reduced/proteins/protein_sequences_pairs.csv"
    "../results/reduced/proteins/protein_sequences_with_mutations.fasta"
    "../results/reduced/proteins/protein_sequences_with_mutations.csv"
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
    echo "🎉 Reduced dataset protein generation completed successfully!"
    echo ""
    echo "Generated datasets:"
    echo "  └─ reduced/proteins/               # Fast generation with pre-validated mutations"
    echo "     ├─ protein_sequences_pairs.*            # Canonical + truncated/extended pairs"
    echo "     └─ protein_sequences_with_mutations.*   # With pre-validated mutations applied"
    echo ""
    echo "Performance benefits:"
    echo "  ├─ ⚡ No mutation re-fetching (used cached results)"
    echo "  ├─ ⚡ No mutation re-validation (used pre-validated impacts)"
    echo "  └─ ⚡ Direct mutation application from step 2 results"
    echo ""
    echo "Dataset composition:"
    echo "  ├─ Pairs: canonical + alternative (truncated/extended) proteins"
    echo "  ├─ Mutations: canonical + alternative + mutated variants (fast mode)"
    echo "  └─ Mutations source: ClinVar missense variants only"
    echo ""
    echo "Next step:"
    echo "  Run: sbatch 4_predict_localization.sh"
else
    echo ""
    echo "❌ Protein generation failed. Some output files are missing."
    echo "Check the generate_proteins.py script and logs for errors."
    exit 1
fi