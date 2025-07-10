#!/bin/bash

#SBATCH --job-name=deeploc_gpu          # Job name
#SBATCH --partition=nvidia-t4-20        # GPU partition (change based on availability)
#SBATCH --ntasks=1                      # Run a single task
#SBATCH --cpus-per-task=8               # CPUs per task (adjust based on GPU type)
#SBATCH --mem=24G                       # Memory for the job
#SBATCH --gres=gpu:1                    # Request 1 GPU
#SBATCH --time=12:00:00                 # Time limit (should be much faster with GPU)
#SBATCH --output=out/deeploc_gpu-%j.out # Standard output log

# Activate conda environment (adjust path as needed)
source ~/.bashrc
# Switch to DeepLoc environment
conda activate deeploc

# Check GPU availability
echo "Checking GPU availability..."
nvidia-smi

# Run DeepLoc predictions with CUDA
echo "Starting DeepLoc predictions for results_reduced at $(date)"
deeploc2 -f results_reduced/protein_sequence_dataset_pairs.fasta -o results_reduced/ -d cuda

echo "Starting DeepLoc predictions for results at $(date)"
deeploc2 -f results/protein_sequence_dataset_pairs.fasta -o results/ -d cuda

# Print completion message
echo "Analysis completed at $(date)"