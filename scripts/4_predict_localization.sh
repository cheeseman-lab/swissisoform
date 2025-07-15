#!/bin/bash

#SBATCH --job-name=deeploc    # Job name
#SBATCH --partition=nvidia-t4-20           # GPU partition
#SBATCH --array=1-4                        # Job array with 4 tasks
#SBATCH --cpus-per-task=8                  # CPUs per task
#SBATCH --mem=48G                          # Memory per task (48G total = 12G × 4)
#SBATCH --gres=gpu:1                       # Request 1 GPU per task
#SBATCH --time=24:00:00                    # Time limit
#SBATCH --output=out/deeploc-%A_%a.out     # %A = job ID, %a = array task ID

# 4_predict_localization.sh
# Generates subcellular localization predictions using DeepLoc

echo "=========================================================="
echo "SwissIsoform Pipeline Step 4: Generate DeepLoc Predictions"
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
echo "=========================================================="

# Check GPU availability
echo "Checking GPU availability..."
nvidia-smi

# Set huggingface cache directory
export HF_HOME="../.cache/huggingface"
mkdir -p "$HF_HOME"

# Only run setup checks on the first task to avoid race conditions
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
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
        else
            echo "✗ $dataset/protein_sequences_pairs.fasta missing"
            missing_files+=("$dataset:pairs")
        fi
        
        if [ -f "$mutations_file" ]; then
            count=$(grep -c '^>' "$mutations_file")
            echo "✓ $dataset/protein_sequences_with_mutations.fasta ($count sequences)"
        else
            echo "✗ $dataset/protein_sequences_with_mutations.fasta missing"
            missing_files+=("$dataset:mutations")
        fi
    done

    # Create localization output directories
    for dataset in "${datasets[@]}"; do
        mkdir -p ../results/$dataset/localization
    done
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

# Define which file each task processes
case $SLURM_ARRAY_TASK_ID in
    1)
        # Task 1: Reduced pairs
        dataset="reduced"
        file_type="pairs"
        input_file="../results/$dataset/proteins/protein_sequences_pairs.fasta"
        ;;
    2)
        # Task 2: Reduced mutations
        dataset="reduced"
        file_type="mutations"
        input_file="../results/$dataset/proteins/protein_sequences_with_mutations.fasta"
        ;;
    3)
        # Task 3: Full pairs
        dataset="full"
        file_type="pairs"
        input_file="../results/$dataset/proteins/protein_sequences_pairs.fasta"
        ;;
    4)
        # Task 4: Full mutations
        dataset="full"
        file_type="mutations"
        input_file="../results/$dataset/proteins/protein_sequences_with_mutations.fasta"
        ;;
    *)
        echo "Unknown array task ID: $SLURM_ARRAY_TASK_ID"
        exit 1
        ;;
esac

# Process the assigned file
echo ""
echo "Array Task ${SLURM_ARRAY_TASK_ID}: Processing $dataset dataset ($file_type sequences)..."
echo "  Input file: $(basename $input_file)"

if [ -f "$input_file" ]; then
    # Remove outputs directory if it exists
    if [ -d "outputs/" ]; then
        echo "Removing existing outputs/ directory..."
        rm -rf outputs/
    fi
    
    # Create descriptive temporary subfolder for this specific run
    temp_subdir="../results/$dataset/localization/${dataset}_${file_type}_temp_$$"
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
    
    # Remove outputs directory if it exists
    if [ -d "outputs/" ]; then
        echo "Removing outputs/ directory..."
        rm -rf outputs/
    fi
    
    echo "  Completed $dataset ($file_type) at $(date)"
else
    echo "  ⚠ Skipping $dataset ($file_type) - input file not found"
    exit 1
fi

# Only verify outputs on the last task to complete
if [ "$SLURM_ARRAY_TASK_ID" -eq 4 ]; then
    # Wait a moment for any file system sync
    sleep 10
    
    echo ""
    echo "Verifying DeepLoc outputs..."

    # Check for outputs from all tasks
    datasets=("reduced" "full")
    file_types=("pairs" "mutations")
    found_outputs=()
    
    for dataset in "${datasets[@]}"; do
        for file_type in "${file_types[@]}"; do
            base_name="protein_sequences_${file_type}"
            
            fast_output="../results/$dataset/localization/${base_name}_Fast_results.csv"
            accurate_output="../results/$dataset/localization/${base_name}_Accurate_results.csv"
            
            if [ -f "$fast_output" ]; then
                count=$(($(wc -l < "$fast_output") - 1))  # Subtract header
                echo "✓ $dataset/$(basename $fast_output) ($count predictions)"
                found_outputs+=("$fast_output")
            else
                echo "✗ $dataset/$(basename $fast_output) missing"
            fi
            
            if [ -f "$accurate_output" ]; then
                count=$(($(wc -l < "$accurate_output") - 1))  # Subtract header
                echo "✓ $dataset/$(basename $accurate_output) ($count predictions)"
                found_outputs+=("$accurate_output")
            else
                echo "✗ $dataset/$(basename $accurate_output) missing"
            fi
        done
    done

    if [ ${#found_outputs[@]} -eq 0 ]; then
        echo ""
        echo "❌ No DeepLoc outputs were generated!"
        echo "Check the SLURM logs for DeepLoc error messages."
        
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
        echo "Summary: Generated ${#found_outputs[@]} prediction files"
        echo ""
        echo "Next step:"
        echo "  Analyze localization predictions for your specific research questions"
    fi
fi