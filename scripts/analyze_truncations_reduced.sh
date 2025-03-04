#!/bin/bash

#SBATCH --job-name=truncations          # Job name
#SBATCH --partition=20                  # Partition name
#SBATCH --ntasks=1                      # Run a single task
#SBATCH --cpus-per-task=12              # Single CPU for the controller job
#SBATCH --mem=24G                       # Memory for the controller job
#SBATCH --time=1:00:00                  # Time limit (hrs:min:sec)
#SBATCH --output=out/truncations-%j.out # Standard output log

# Activate conda environment (adjust path as needed)
source ~/.bashrc
conda activate swissisoform

# Create output directory
mkdir -p results_reduced

# Define paths for command-line arguments
GENE_LIST="../data/ribosome_profiling/gene_list_reduced.txt"
OUTPUT_DIR="results_reduced"
GENOME_PATH="../data/genome_data/GRCh38.p7.genome.fa"
ANNOTATION_PATH="../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf"
BED_PATH="../data/ribosome_profiling/selected_truncations_JL_cleaned.bed"
PREFERRED_TRANSCRIPTS="../data/genome_data/hela_top_transcript.txt"

# Run the analysis script with all arguments
python3 analyze_truncations.py "$GENE_LIST" "$OUTPUT_DIR" \
  --genome "$GENOME_PATH" \
  --annotation "$ANNOTATION_PATH" \
  --bed "$BED_PATH" \
  --preferred-transcripts "$PREFERRED_TRANSCRIPTS" \
  --visualize

# Print completion message
echo "Analysis completed at $(date)"