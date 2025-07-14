#!/bin/bash

#SBATCH --job-name=swissisoform_deeploc    # Job name
#SBATCH --partition=nvidia-t4-20           # GPU partition
#SBATCH --ntasks=1                         # Run a single task
#SBATCH --cpus-per-task=8                  # CPUs per task
#SBATCH --mem=48G                          # Memory for the job
#SBATCH --gres=gpu:1                       # Request 1 GPU
#SBATCH --time=12:00:00                    # Time limit
#SBATCH --output=out/deeploc-%j.out        # Standard output log

# 4_generate_deeploc.sh
# Generates subcellular localization predictions using DeepLoc

echo "=========================================================="
echo "SwissIsoform Pipeline Step 4: Generate DeepLoc Predictions"
echo "=========================================================="

# Check GPU availability
echo "Checking GPU availability..."
nvidia-smi

# Set huggingface cache directory
export HF_HOME="../.cache/huggingface"
mkdir -p "$HF_HOME"

# Check if protein sequence files exist
echo ""
echo "Checking for input files..."

input_files=(
    "../results/reduced/protein_sequences.fasta"
    "../results/full/protein_sequences.fasta"
)

missing_files=false
for file in "${input_files[@]}"; do
    if [ -f "$file" ]; then
        count=$(grep -c '^>' "$file")
        echo "✓ $(basename $file) ($count sequences)"
    else
        echo "✗ $(basename $file) missing"
        missing_files=true
    fi
done

if [ "$missing_files" = true ]; then
    echo ""
    echo "❌ Missing input files! Run 3_generate_proteins.sh first"
    exit 1
fi

# Activate DeepLoc environment
echo ""
echo "Activating DeepLoc environment..."
source ~/.bashrc
conda activate deeploc || {
    echo "❌ Failed to activate deeploc environment"
    echo ""
    echo "To install DeepLoc:"
    echo "  conda create -n deeploc python=3.8"
    echo "  conda activate deeploc"
    echo "  pip install DeepLoc-2.1.0.tar.gz"
    exit 1
}

# Run DeepLoc predictions
echo ""
echo "Running DeepLoc predictions..."

datasets=("reduced" "full")

for dataset in "${datasets[@]}"; do
    echo ""
    echo "Processing $dataset dataset..."
    
    input_file="../results/$dataset/protein_sequences.fasta"
    output_dir="../results/$dataset"
    
    if [ -f "$input_file" ]; then
        echo "  Starting DeepLoc Fast mode for $dataset at $(date)"
        deeploc2 -f "$input_file" -m Fast -o "$output_dir/" -d cuda
        
        echo "  Starting DeepLoc Accurate mode for $dataset at $(date)"
        deeploc2 -f "$input_file" -m Accurate -o "$output_dir/" -d cuda
        
        echo "  Completed $dataset dataset at $(date)"
    else
        echo "  ⚠ Skipping $dataset - input file not found"
    fi
done

# Verify outputs
echo ""
echo "Verifying DeepLoc outputs..."

expected_outputs=(
    "../results/reduced/protein_sequences_results.csv"
    "../results/reduced/protein_sequences_Accurate_results.csv"
    "../results/full/protein_sequences_results.csv"
    "../results/full/protein_sequences_Accurate_results.csv"
)

all_outputs_present=true
for file in "${expected_outputs[@]}"; do
    if [ -f "$file" ]; then
        count=$(($(wc -l < "$file") - 1))  # Subtract header
        echo "✓ $(basename $file) ($count predictions)"
    else
        echo "✗ $(basename $file) missing"
        all_outputs_present=false
    fi
done

if [ "$all_outputs_present" = true ]; then
    echo ""
    echo "🎉 DeepLoc predictions completed successfully!"
    echo ""
    echo "Generated predictions:"
    echo "  ├─ reduced/"
    echo "  │  ├─ protein_sequences_results.csv         # Fast mode"
    echo "  │  └─ protein_sequences_Accurate_results.csv # Accurate mode"
    echo "  └─ full/"
    echo "     ├─ protein_sequences_results.csv         # Fast mode"
    echo "     └─ protein_sequences_Accurate_results.csv # Accurate mode"
    echo ""
    echo "Next step:"
    echo "  Run analysis scripts for your specific research questions"
else
    echo ""
    echo "❌ DeepLoc prediction failed. Some output files are missing."
    echo "Check the SLURM output file for detailed error messages."
    exit 1
fi