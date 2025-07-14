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
        # Create temporary working directory for this prediction
        temp_dir="../results/$dataset/localization/temp_$(basename $input_file .fasta)"
        mkdir -p "$temp_dir"
        
        # Copy input file to temp directory (DeepLoc outputs in same dir as input)
        temp_input="$temp_dir/$(basename $input_file)"
        cp "$input_file" "$temp_input"
        
        echo "  Starting DeepLoc Fast mode for $dataset ($file_type) at $(date)"
        deeploc2 -f "$temp_input" -m Fast -d cuda
        
        echo "  Starting DeepLoc Accurate mode for $dataset ($file_type) at $(date)"
        deeploc2 -f "$temp_input" -m Accurate -d cuda
        
        # Move results to organized location and rename appropriately
        base_name=$(basename $input_file .fasta)
        final_dir="../results/$dataset/localization"
        
        if [ -f "$temp_dir/${base_name}_results.csv" ]; then
            mv "$temp_dir/${base_name}_results.csv" "$final_dir/${base_name}_Fast_results.csv"
            echo "  ✓ Moved Fast results to $final_dir/${base_name}_Fast_results.csv"
        fi
        
        if [ -f "$temp_dir/${base_name}_Accurate_results.csv" ]; then
            mv "$temp_dir/${base_name}_Accurate_results.csv" "$final_dir/${base_name}_Accurate_results.csv"
            echo "  ✓ Moved Accurate results to $final_dir/${base_name}_Accurate_results.csv"
        fi
        
        # Clean up temp directory
        rm -rf "$temp_dir"
        
        echo "  Completed $dataset ($file_type) at $(date)"
    else
        echo "  ⚠ Skipping $dataset ($file_type) - input file not found"
    fi
done

# Verify outputs
echo ""
echo "Verifying DeepLoc outputs..."

# Check for all possible output combinations
expected_outputs=()
for dataset in "${datasets[@]}"; do
    for file_type in "pairs" "mutations"; do
        base_name="protein_sequences_${file_type}"
        expected_outputs+=("../results/$dataset/localization/${base_name}_Fast_results.csv")
        expected_outputs+=("../results/$dataset/localization/${base_name}_Accurate_results.csv")
    done
done

all_outputs_present=true
found_outputs=()

for file in "${expected_outputs[@]}"; do
    if [ -f "$file" ]; then
        count=$(($(wc -l < "$file") - 1))  # Subtract header
        echo "✓ $(basename $file) ($count predictions)"
        found_outputs+=("$file")
    fi
done

if [ ${#found_outputs[@]} -eq 0 ]; then
    echo "❌ No DeepLoc outputs were generated!"
    all_outputs_present=false
else
    echo ""
    echo "🎉 DeepLoc predictions completed successfully!"
    echo ""
    echo "Generated predictions:"
    
    for dataset in "${datasets[@]}"; do
        echo "  ├─ $dataset/localization/"
        
        # Check pairs files
        pairs_fast="../results/$dataset/localization/protein_sequences_pairs_Fast_results.csv"
        pairs_accurate="../results/$dataset/localization/protein_sequences_pairs_Accurate_results.csv"
        mutations_fast="../results/$dataset/localization/protein_sequences_with_mutations_Fast_results.csv"
        mutations_accurate="../results/$dataset/localization/protein_sequences_with_mutations_Accurate_results.csv"
        
        if [ -f "$pairs_fast" ]; then
            echo "  │  ├─ protein_sequences_pairs_Fast_results.csv"
        fi
        if [ -f "$pairs_accurate" ]; then
            echo "  │  ├─ protein_sequences_pairs_Accurate_results.csv"
        fi
        if [ -f "$mutations_fast" ]; then
            echo "  │  ├─ protein_sequences_with_mutations_Fast_results.csv"
        fi
        if [ -f "$mutations_accurate" ]; then
            echo "  │  └─ protein_sequences_with_mutations_Accurate_results.csv"
        fi
    done
fi

if [ "$all_outputs_present" = true ] || [ ${#found_outputs[@]} -gt 0 ]; then
    echo ""
    echo "Next step:"
    echo "  Analyze localization predictions for your specific research questions"
else
    echo ""
    echo "❌ DeepLoc prediction failed. Check the SLURM output file for detailed error messages."
    exit 1
fi