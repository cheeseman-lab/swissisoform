#!/bin/bash

#SBATCH --job-name=swissisoform_deeploc    # Job name
#SBATCH --partition=nvidia-t4-20           # GPU partition
#SBATCH --ntasks=1                         # Run a single task
#SBATCH --cpus-per-task=8                  # CPUs per task
#SBATCH --mem=48G                          # Memory for the job
#SBATCH --gres=gpu:1                       # Request 1 GPU
#SBATCH --time=12:00:00                    # Time limit
#SBATCH --output=out/deeploc-%j.out        # Standard output log

# 4_predict_localization.sh
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

# Check what protein sequence files are available
echo ""
echo "Checking for input files..."

datasets=("reduced" "full")
available_files=()
missing_files=()

for dataset in "${datasets[@]}"; do
    pairs_file="../results/$dataset/proteins/protein_sequences_pairs.fasta"
    mutations_file="../results/$dataset/proteins/protein_sequences_with_mutations.fasta"
    
    if [ -f "$pairs_file" ]; then
        count=$(grep -c '^>' "$pairs_file")
        echo "✓ $dataset/protein_sequences_pairs.fasta ($count sequences)"
        available_files+=("$dataset:pairs:$pairs_file")
    else
        echo "✗ $dataset/protein_sequences_pairs.fasta missing"
        missing_files+=("$dataset:pairs")
    fi
    
    if [ -f "$mutations_file" ]; then
        count=$(grep -c '^>' "$mutations_file")
        echo "✓ $dataset/protein_sequences_with_mutations.fasta ($count sequences)"
        available_files+=("$dataset:mutations:$mutations_file")
    else
        echo "✗ $dataset/protein_sequences_with_mutations.fasta missing"
        missing_files+=("$dataset:mutations")
    fi
done

if [ ${#available_files[@]} -eq 0 ]; then
    echo ""
    echo "❌ No input files found! Run 3_generate_proteins.sh first"
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
    echo "  1. Download DeepLoc 2.1 from: https://services.healthtech.dtu.dk/services/DeepLoc-2.0/"
    echo "  2. Place DeepLoc-2.1.0.tar.gz in the main swissisoform/ directory"
    echo "  3. conda create -n deeploc python=3.8"
    echo "  4. conda activate deeploc"
    echo "  5. pip install DeepLoc-2.1.0.tar.gz"
    exit 1
}

# Create localization output directories
for dataset in "${datasets[@]}"; do
    mkdir -p ../results/$dataset/localization
done

# Run DeepLoc predictions for each available file
echo ""
echo "Running DeepLoc predictions..."

for file_info in "${available_files[@]}"; do
    IFS=':' read -r dataset file_type input_file <<< "$file_info"
    
    echo ""
    echo "Processing $dataset dataset ($file_type sequences)..."
    echo "  Input file: $(basename $input_file)"
    
    if [ -f "$input_file" ]; then
        # Create descriptive temporary subfolder for this specific run
        temp_subdir="../results/$dataset/localization/${dataset}_${file_type}_temp"
        mkdir -p "$temp_subdir"
        
        echo "  Starting DeepLoc Fast mode for $dataset ($file_type) at $(date)"
        deeploc2 -f "$input_file" -m Fast -o "$temp_subdir/" -d cuda
        
        # Find and move the Fast results file
        fast_result=$(find "$temp_subdir" -name "results_*.csv" | head -n 1)
        if [ -n "$fast_result" ] && [ -f "$fast_result" ]; then
            mv "$fast_result" "../results/$dataset/localization/protein_sequences_${file_type}_Fast_results.csv"
            echo "  ✓ Moved Fast results to protein_sequences_${file_type}_Fast_results.csv"
        else
            echo "  ✗ Fast results not found"
        fi
        
        echo "  Starting DeepLoc Accurate mode for $dataset ($file_type) at $(date)"
        deeploc2 -f "$input_file" -m Accurate -o "$temp_subdir/" -d cuda
        
        # Find and move the Accurate results file
        accurate_result=$(find "$temp_subdir" -name "results_*.csv" | head -n 1)
        if [ -n "$accurate_result" ] && [ -f "$accurate_result" ]; then
            mv "$accurate_result" "../results/$dataset/localization/protein_sequences_${file_type}_Accurate_results.csv"
            echo "  ✓ Moved Accurate results to protein_sequences_${file_type}_Accurate_results.csv"
        else
            echo "  ✗ Accurate results not found"
        fi
        
        # Clean up temp subfolder
        rm -rf "$temp_subdir"
        
        echo "  Completed $dataset ($file_type) at $(date)"
    else
        echo "  ⚠ Skipping $dataset ($file_type) - input file not found"
    fi
done

# Verify outputs
echo ""
echo "Verifying DeepLoc outputs..."

# Check for outputs based on what files were actually processed
found_outputs=()
processed_files=()

# Build list of files that should have been processed
for file_info in "${available_files[@]}"; do
    IFS=':' read -r dataset file_type input_file <<< "$file_info"
    base_name="protein_sequences_${file_type}"
    
    fast_output="../results/$dataset/localization/${base_name}_Fast_results.csv"
    accurate_output="../results/$dataset/localization/${base_name}_Accurate_results.csv"
    
    if [ -f "$fast_output" ]; then
        count=$(($(wc -l < "$fast_output") - 1))  # Subtract header
        echo "✓ $(basename $fast_output) ($count predictions)"
        found_outputs+=("$fast_output")
    else
        echo "✗ $(basename $fast_output) missing"
    fi
    
    if [ -f "$accurate_output" ]; then
        count=$(($(wc -l < "$accurate_output") - 1))  # Subtract header
        echo "✓ $(basename $accurate_output) ($count predictions)"
        found_outputs+=("$accurate_output")
    else
        echo "✗ $(basename $accurate_output) missing"
    fi
    
    processed_files+=("$dataset:$file_type")
done

if [ ${#found_outputs[@]} -eq 0 ]; then
    echo ""
    echo "❌ No DeepLoc outputs were generated!"
    echo "Check the SLURM log for DeepLoc error messages."
    
    # Debug: List what files exist in localization directories
    echo ""
    echo "Debug: Files in localization directories:"
    for dataset in "${datasets[@]}"; do
        localization_dir="../results/$dataset/localization"
        if [ -d "$localization_dir" ]; then
            echo "  $dataset/localization/:"
            ls -la "$localization_dir" | sed 's/^/    /'
        fi
    done
    
    exit 1
else
    echo ""
    echo "🎉 DeepLoc predictions completed successfully!"
    echo ""
    echo "Generated predictions:"
    
    # Group outputs by dataset
    for dataset in "${datasets[@]}"; do
        dataset_outputs=($(printf '%s\n' "${found_outputs[@]}" | grep "$dataset/localization/"))
        
        if [ ${#dataset_outputs[@]} -gt 0 ]; then
            echo "  ├─ $dataset/localization/"
            for output in "${dataset_outputs[@]}"; do
                echo "  │  ├─ $(basename $output)"
            done
        fi
    done
    
    echo ""
    echo "Summary: Generated ${#found_outputs[@]} prediction files from ${#processed_files[@]} input files"
    echo ""
    echo "Next step:"
    echo "  Analyze localization predictions for your specific research questions"
fi