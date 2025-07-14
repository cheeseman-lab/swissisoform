#!/bin/bash

#SBATCH --job-name=translate              # Job name
#SBATCH --partition=20                    # Partition name
#SBATCH --ntasks=1                        # Run a single task
#SBATCH --cpus-per-task=12                # Single CPU for the controller job
#SBATCH --mem=24G                         # Memory for the controller job
#SBATCH --time=24:00:00                   # Time limit (hrs:min:sec)
#SBATCH --output=out/translate-%j.out     # Standard output log

# Activate conda environment (adjust path as needed)
source ~/.bashrc
conda activate swissisoform

# Define common paths
GENOME_PATH="../data/genome_data/GRCh38.p7.genome.fa"
ANNOTATION_PATH="../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf"
PREFERRED_TRANSCRIPTS="../data/genome_data/hela_top_transcript.txt"
MIN_LENGTH=10
MAX_LENGTH=100000
FORMAT="fasta,csv"

# Create results directory outside of scripts
mkdir -p ../results

echo "Starting protein sequence translation at $(date)"

# Run reduced dataset (selected truncations)
echo "=== Processing reduced dataset ==="
GENE_LIST="../data/ribosome_profiling/gene_list_reduced.txt"
BED_PATH="../data/ribosome_profiling/selected_truncations_JL_cleaned.bed"
OUTPUT_DIR="../results/reduced"

python3 translate.py "$GENE_LIST" "$OUTPUT_DIR" \
  --genome "$GENOME_PATH" \
  --annotation "$ANNOTATION_PATH" \
  --bed "$BED_PATH" \
  --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
  --min-length "$MIN_LENGTH" \
  --max-length "$MAX_LENGTH" \
  --format "$FORMAT" \
  --include-canonical \
  --pairs-only

echo "Reduced dataset completed at $(date)"

# Run reduced dataset with mutations
echo "=== Processing reduced dataset with mutations ==="
OUTPUT_DIR="../results/reduced_mutations"

python3 translate.py "$GENE_LIST" "$OUTPUT_DIR" \
  --genome "$GENOME_PATH" \
  --annotation "$ANNOTATION_PATH" \
  --bed "$BED_PATH" \
  --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
  --min-length "$MIN_LENGTH" \
  --max-length "$MAX_LENGTH" \
  --format "$FORMAT" \
  --include-canonical \
  --include-mutations \
  --impact-types "missense variant"

echo "Reduced dataset with mutations completed at $(date)"

# Run full dataset (all truncations) - uncomment when ready
echo "=== Processing full dataset ==="
GENE_LIST="../data/ribosome_profiling/gene_list.txt"
BED_PATH="../data/ribosome_profiling/full_truncations_JL_cleaned.bed"
OUTPUT_DIR="../results/full"

python3 translate.py "$GENE_LIST" "$OUTPUT_DIR" \
  --genome "$GENOME_PATH" \
  --annotation "$ANNOTATION_PATH" \
  --bed "$BED_PATH" \
  --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
  --min-length "$MIN_LENGTH" \
  --max-length "$MAX_LENGTH" \
  --format "$FORMAT" \
  --include-canonical \
  --pairs-only

echo "Full dataset completed at $(date)"

Run full dataset with mutations - uncomment when ready
echo "=== Processing full dataset with mutations ==="
OUTPUT_DIR="../results/full_mutations"

python3 translate.py "$GENE_LIST" "$OUTPUT_DIR" \
  --genome "$GENOME_PATH" \
  --annotation "$ANNOTATION_PATH" \
  --bed "$BED_PATH" \
  --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
  --min-length "$MIN_LENGTH" \
  --max-length "$MAX_LENGTH" \
  --format "$FORMAT" \
  --include-canonical \
  --include-mutations \
  --impact-types "missense variant"

echo "Full dataset with mutations completed at $(date)"

echo "All protein sequence translation completed at $(date)"

# Print summary of results
echo ""
echo "=== Results Summary ==="
echo "Results saved to ../results/ directory:"
echo "  ├─ reduced/                     # Reduced dataset (pairs only)"
echo "  ├─ reduced_mutations/           # Reduced dataset with mutations"
echo "  ├─ full/                        # Full dataset (pairs only)"
echo "  └─ full_mutations/              # Full dataset with mutations"

echo ""
echo "Output files in each directory:"
echo "  ├─ protein_sequence_dataset_pairs.fasta"
echo "  ├─ protein_sequence_dataset_pairs.csv"
echo "  ├─ protein_sequences_with_mutations.fasta (mutations only)"
echo "  └─ protein_sequences_with_mutations.csv (mutations only)"