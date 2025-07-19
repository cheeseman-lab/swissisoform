#!/bin/bash
#SBATCH --job-name=deeploc                 # Job name
#SBATCH --partition=nvidia-A4000-20        # Faster GPU partition
#SBATCH --array=1-8                        # Split into 8 tasks for better parallelization
#SBATCH --cpus-per-task=4                  # CPUs 
#SBATCH --mem=36G                          # Memory per task
#SBATCH --gres=gpu:1                       # Request 1 GPU per task
#SBATCH --time=36:00:00                    # Time limit
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

# Define which file and mode each task processes (8 tasks total)
case $SLURM_ARRAY_TASK_ID in
    1)
        dataset="reduced"; file_type="pairs"; mode="Fast"
        input_file="../results/$dataset/proteins/protein_sequences_pairs.fasta"
        ;;
    2)
        dataset="reduced"; file_type="pairs"; mode="Accurate"
        input_file="../results/$dataset/proteins/protein_sequences_pairs.fasta"
        ;;
    3)
        dataset="reduced"; file_type="mutations"; mode="Fast"
        input_file="../results/$dataset/proteins/protein_sequences_with_mutations.fasta"
        ;;
    4)
        dataset="reduced"; file_type="mutations"; mode="Accurate"
        input_file="../results/$dataset/proteins/protein_sequences_with_mutations.fasta"
        ;;
    5)
        dataset="full"; file_type="pairs"; mode="Fast"
        input_file="../results/$dataset/proteins/protein_sequences_pairs.fasta"
        ;;
    6)
        dataset="full"; file_type="pairs"; mode="Accurate"
        input_file="../results/$dataset/proteins/protein_sequences_pairs.fasta"
        ;;
    7)
        dataset="full"; file_type="mutations"; mode="Fast"
        input_file="../results/$dataset/proteins/protein_sequences_with_mutations.fasta"
        ;;
    8)
        dataset="full"; file_type="mutations"; mode="Accurate"
        input_file="../results/$dataset/proteins/protein_sequences_with_mutations.fasta"
        ;;
    *)
        echo "Unknown array task ID: $SLURM_ARRAY_TASK_ID"
        exit 1
        ;;
esac

# Process the assigned file with specified mode
echo ""
echo "Array Task ${SLURM_ARRAY_TASK_ID}: Processing $dataset dataset ($file_type sequences, $mode mode)..."
echo "  Input file: $(basename $input_file)"

if [ -f "$input_file" ]; then
    # Remove outputs directory if it exists
    if [ -d "outputs/" ]; then
        echo "Removing existing outputs/ directory..."
        rm -rf outputs/
    fi
    
    # Create descriptive temporary subfolder for this specific run
    temp_subdir="../results/$dataset/localization/${dataset}_${file_type}_${mode}_temp_$$"
    mkdir -p "$temp_subdir"
    
    echo "  Starting DeepLoc $mode mode for $dataset ($file_type) at $(date)"
    
    # Set GPU memory growth to avoid OOM errors
    export TF_FORCE_GPU_ALLOW_GROWTH=true
    export CUDA_VISIBLE_DEVICES=0
    
    deeploc2 -f "$input_file" -m "$mode" -o "$temp_subdir/" -d cuda
    
    # Find and move the results file
    result_file=$(find "$temp_subdir" -name "results_*.csv" | head -n 1)
    if [ -n "$result_file" ] && [ -f "$result_file" ]; then
        output_file="../results/$dataset/localization/protein_sequences_${file_type}_${mode}_results.csv"
        mv "$result_file" "$output_file"
        echo "  ✓ Moved $mode results to protein_sequences_${file_type}_${mode}_results.csv"
    else
        echo "  ✗ $mode results not found"
        exit 1
    fi
    
    # Clean up temp subfolder
    rm -rf "$temp_subdir"
    
    # Remove outputs directory if it exists
    if [ -d "outputs/" ]; then
        echo "Removing outputs/ directory..."
        rm -rf outputs/
    fi
    
    echo "  Completed $dataset ($file_type, $mode) at $(date)"
else
    echo "  ⚠ Skipping $dataset ($file_type, $mode) - input file not found"
    exit 1
fi

# Only verify outputs on the last task to complete
if [ "$SLURM_ARRAY_TASK_ID" -eq 8 ]; then
    # Wait a moment for any file system sync
    sleep 10
    
    echo ""
    echo "Verifying DeepLoc outputs..."
    # Check for outputs from all tasks
    datasets=("reduced" "full")
    file_types=("pairs" "mutations")
    modes=("Fast" "Accurate")
    found_outputs=()
    
    for dataset in "${datasets[@]}"; do
        for file_type in "${file_types[@]}"; do
            for mode in "${modes[@]}"; do
                output_file="../results/$dataset/localization/protein_sequences_${file_type}_${mode}_results.csv"
                
                if [ -f "$output_file" ]; then
                    count=$(($(wc -l < "$output_file") - 1))  # Subtract header
                    echo "✓ $dataset/$(basename $output_file) ($count predictions)"
                    found_outputs+=("$output_file")
                else
                    echo "✗ $dataset/$(basename $output_file) missing"
                fi
            done
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