#!/bin/bash

#SBATCH --job-name=deeploc          # Job name
#SBATCH --partition=20                  # Partition name
#SBATCH --ntasks=1                      # Run a single task
#SBATCH --cpus-per-task=12              # Single CPU for the controller job
#SBATCH --mem=32G                       # Memory for the controller job
#SBATCH --time=36:00:00                  # Time limit (hrs:min:sec)
#SBATCH --output=out/deeploc-%j.out # Standard output log

# Activate conda environment (adjust path as needed)
source ~/.bashrc
# Switch to DeepLoc environment
conda activate deeploc

# Run DeepLoc predictions
echo "Starting DeepLoc predictions for results reduced at $(date)"
deeploc2 -f results_reduced/protein_sequence_dataset_pairs.fasta -o results_reduced/

echo "Starting DeepLoc predictions for results at $(date)"
deeploc2 -f results/protein_sequence_dataset_pairs.fasta -o results/

# Print completion message
echo "Analysis completed at $(date)"