#!/bin/bash

#SBATCH --job-name=generate_protein_sequences          # Job name
#SBATCH --partition=20                  # Partition name
#SBATCH --ntasks=1                      # Run a single task
#SBATCH --cpus-per-task=12              # Single CPU for the controller job
#SBATCH --mem=12G                       # Memory for the controller job
#SBATCH --time=24:00:00                 # Time limit (hrs:min:sec)
#SBATCH --output=out/generate_protein_sequences-%j.out     # Standard output log

# Activate conda environment (adjust path as needed)
source ~/.bashrc
conda activate swissisoform

# Run the truncations script
mkdir -p results

# Run the analysis script
python3 generate_protein_sequences.py '../data/ribosome_profiling/gene_list.txt' 'results'