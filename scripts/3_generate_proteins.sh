#!/bin/bash

#SBATCH --job-name=proteins_chunked        # Job name
#SBATCH --partition=20                     # Partition name
#SBATCH --array=1-8                        # 8 chunks for full dataset
#SBATCH --cpus-per-task=4                  # CPUs per task
#SBATCH --mem=16G                          # Memory per task
#SBATCH --time=24:00:00                    # Time limit (hrs:min:sec)
#SBATCH --output=out/proteins-%A_%a.out    # %A = job ID, %a = array task ID

# 3_generate_proteins_chunked.sh
# Generates protein sequences using chunked parallel processing with pre-validated mutations

echo "======================================================="
echo "SwissIsoform Pipeline Step 3: Generate Proteins (Chunked)"
echo "Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
echo "======================================================="

# Function to split gene list into chunks
split_gene_list() {
    local input_file=$1
    local total_chunks=$2
    local chunk_id=$3
    local output_file=$4
    
    total_genes=$(wc -l < "$input_file")
    genes_per_chunk=$(( (total_genes + total_chunks - 1) / total_chunks ))
    start_line=$(( (chunk_id - 1) * genes_per_chunk + 1 ))
    end_line=$(( chunk_id * genes_per_chunk ))
    
    sed -n "${start_line},${end_line}p" "$input_file" > "$output_file"
}

# Only run setup checks on the first task to avoid race conditions
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    # Check if required input files exist
    echo "Checking for required input files..."

    GENOME_DIR="../data/genome_data"
    RIBOPROF_DIR="../data/ribosome_profiling"
    MUTATIONS_DIR="../results/full/mutations"

    required_files=(
        "$GENOME_DIR/GRCh38.p7.genome.fa"
        "$GENOME_DIR/gencode.v25.annotation.ensembl_cleaned.gtf"
        "$RIBOPROF_DIR/isoforms_with_transcripts.bed"
        "$RIBOPROF_DIR/isoforms_gene_list.txt"
        "$MUTATIONS_DIR/isoform_level_results.csv"
    )

    echo ""
    echo "Checking required files..."
    missing_files=false
    for file in "${required_files[@]}"; do
        if [ -f "$file" ]; then
            echo "✓ $(basename $file)"
        else
            echo "✗ $(basename $file) missing"
            missing_files=true
        fi
    done

    if [ "$missing_files" = true ]; then
        echo ""
        echo "❌ Missing required files! Run 1_cleanup_files.sh and 2_analyze_mutations.sh first"
        exit 1
    fi

    # Check that we have mutation results
    mutation_count=$(tail -n +2 "$MUTATIONS_DIR/isoform_level_results.csv" | wc -l)
    if [ "$mutation_count" -eq 0 ]; then
        echo "❌ No mutation results found! Run 2_analyze_mutations.sh first"
        exit 1
    fi

    echo "✓ Found $mutation_count validated mutation results"

    # Create results directory structure
    mkdir -p ../results/full/proteins
    mkdir -p ../results/temp/protein_chunks
fi

# Wait for setup to complete and add staggered delays
sleep 2
STAGGER_DELAY=$(( (SLURM_ARRAY_TASK_ID - 1) * 3 ))
echo "Adding ${STAGGER_DELAY}s staggered delay for task ${SLURM_ARRAY_TASK_ID}"
sleep $STAGGER_DELAY

# Activate conda environment
echo ""
echo "Activating conda environment..."
source ~/.bashrc
conda activate swissisoform || {
    echo "❌ Failed to activate swissisoform conda environment"
    echo "Please run: conda env create --file=../environment.yml"
    exit 1
}

# Define common paths
GENOME_PATH="../data/genome_data/GRCh38.p7.genome.fa"
ANNOTATION_PATH="../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf"
TRUNCATIONS_PATH="../data/ribosome_profiling/isoforms_with_transcripts.bed"
MUTATIONS_FILE="../results/full/mutations/isoform_level_results.csv"
GENE_LIST="../data/ribosome_profiling/isoforms_gene_list.txt"
OUTPUT_DIR="../results/full/proteins/chunk_${SLURM_ARRAY_TASK_ID}"
CHUNK_ID=$SLURM_ARRAY_TASK_ID
TOTAL_CHUNKS=8
MIN_LENGTH=10
MAX_LENGTH=100000
FORMAT="fasta,csv"

# Create task-specific output directory
mkdir -p "$OUTPUT_DIR"

# Create gene chunk for this task
CHUNK_FILE="../results/temp/protein_chunks/chunk_${CHUNK_ID}.txt"
split_gene_list "$GENE_LIST" "$TOTAL_CHUNKS" "$CHUNK_ID" "$CHUNK_FILE"

# Skip if no genes in chunk
if [ ! -s "$CHUNK_FILE" ]; then
    echo "No genes in chunk ${CHUNK_ID}, exiting"
    exit 0
fi

echo "Processing chunk ${CHUNK_ID}"
echo "Gene list: $(wc -l < "$CHUNK_FILE") genes"
echo "Configuration:"
echo "  ├─ Pre-validated mutations: $(basename $MUTATIONS_FILE)"
echo "  ├─ Sources: clinvar"
echo "  ├─ Impact types: missense variant"
echo "  ├─ Length range: $MIN_LENGTH-$MAX_LENGTH amino acids"
echo "  └─ Output format: $FORMAT"

# Generate both pairs and mutations datasets for this chunk
echo ""
echo "Generating protein sequences for chunk ${CHUNK_ID}..."
python3 generate_proteins.py "$CHUNK_FILE" "$OUTPUT_DIR" \
  --mutations-file "$MUTATIONS_FILE" \
  --genome "$GENOME_PATH" \
  --annotation "$ANNOTATION_PATH" \
  --bed "$TRUNCATIONS_PATH" \
  --sources clinvar \
  --impact-types "missense variant" \
  --min-length "$MIN_LENGTH" \
  --max-length "$MAX_LENGTH" \
  --format "$FORMAT" \
  --fast-mode

exit_code=$?

if [ $exit_code -ne 0 ]; then
    echo "❌ Protein generation failed for chunk ${CHUNK_ID} with exit code $exit_code"
    exit 1
fi

echo "Completed chunk ${CHUNK_ID}"

# Clean up chunk file
rm -f "$CHUNK_FILE"

# Verification and merging (only last task)
if [ "$SLURM_ARRAY_TASK_ID" -eq 8 ]; then
    # Wait for other tasks to finish
    echo "Waiting for all tasks to complete..."
    sleep 30
    
    echo ""
    echo "Merging results from all chunks..."
    
    # Create final output directory
    FINAL_OUTPUT_DIR="../results/full/proteins"
    mkdir -p "$FINAL_OUTPUT_DIR"
    
    # Merge pairs datasets
    echo "Merging pairs datasets..."
    
    # Merge FASTA files (pairs)
    echo "  ├─ Merging protein_sequences_pairs.fasta..."
    > "$FINAL_OUTPUT_DIR/protein_sequences_pairs.fasta"
    for i in {1..8}; do
        chunk_fasta="../results/full/proteins/chunk_${i}/protein_sequences_pairs.fasta"
        if [ -f "$chunk_fasta" ]; then
            cat "$chunk_fasta" >> "$FINAL_OUTPUT_DIR/protein_sequences_pairs.fasta"
            count=$(grep -c '^>' "$chunk_fasta" 2>/dev/null || echo 0)
            echo "    ├─ Added chunk $i ($count sequences)"
        else
            echo "    ├─ Warning: chunk $i pairs FASTA missing"
        fi
    done
    
    # Merge CSV files (pairs)
    echo "  ├─ Merging protein_sequences_pairs.csv..."
    first_file=true
    for i in {1..8}; do
        chunk_csv="../results/full/proteins/chunk_${i}/protein_sequences_pairs.csv"
        if [ -f "$chunk_csv" ]; then
            if [ "$first_file" = true ]; then
                # Copy first file with header
                cp "$chunk_csv" "$FINAL_OUTPUT_DIR/protein_sequences_pairs.csv"
                first_file=false
                rows=$(wc -l < "$chunk_csv")
                echo "    ├─ Added chunk $i with header ($rows rows)"
            else
                # Append without header
                tail -n +2 "$chunk_csv" >> "$FINAL_OUTPUT_DIR/protein_sequences_pairs.csv"
                rows=$(($(wc -l < "$chunk_csv") - 1))
                echo "    ├─ Added chunk $i ($rows rows)"
            fi
        else
            echo "    ├─ Warning: chunk $i pairs CSV missing"
        fi
    done
    
    # Merge mutations datasets
    echo "Merging mutations datasets..."
    
    # Merge FASTA files (mutations)
    echo "  ├─ Merging protein_sequences_with_mutations.fasta..."
    > "$FINAL_OUTPUT_DIR/protein_sequences_with_mutations.fasta"
    for i in {1..8}; do
        chunk_fasta="../results/full/proteins/chunk_${i}/protein_sequences_with_mutations.fasta"
        if [ -f "$chunk_fasta" ]; then
            cat "$chunk_fasta" >> "$FINAL_OUTPUT_DIR/protein_sequences_with_mutations.fasta"
            count=$(grep -c '^>' "$chunk_fasta" 2>/dev/null || echo 0)
            echo "    ├─ Added chunk $i ($count sequences)"
        else
            echo "    ├─ Warning: chunk $i mutations FASTA missing"
        fi
    done
    
    # Merge CSV files (mutations)
    echo "  ├─ Merging protein_sequences_with_mutations.csv..."
    first_file=true
    for i in {1..8}; do
        chunk_csv="../results/full/proteins/chunk_${i}/protein_sequences_with_mutations.csv"
        if [ -f "$chunk_csv" ]; then
            if [ "$first_file" = true ]; then
                # Copy first file with header
                cp "$chunk_csv" "$FINAL_OUTPUT_DIR/protein_sequences_with_mutations.csv"
                first_file=false
                rows=$(wc -l < "$chunk_csv")
                echo "    ├─ Added chunk $i with header ($rows rows)"
            else
                # Append without header
                tail -n +2 "$chunk_csv" >> "$FINAL_OUTPUT_DIR/protein_sequences_with_mutations.csv"
                rows=$(($(wc -l < "$chunk_csv") - 1))
                echo "    ├─ Added chunk $i ($rows rows)"
            fi
        else
            echo "    ├─ Warning: chunk $i mutations CSV missing"
        fi
    done
    
    # Verify final outputs
    echo ""
    echo "Verifying merged protein datasets..."

    expected_files=(
        "../results/full/proteins/protein_sequences_pairs.fasta"
        "../results/full/proteins/protein_sequences_pairs.csv"
        "../results/full/proteins/protein_sequences_with_mutations.fasta"
        "../results/full/proteins/protein_sequences_with_mutations.csv"
    )

    all_files_present=true
    for file in "${expected_files[@]}"; do
        if [ -f "$file" ]; then
            if [[ "$file" == *.fasta ]]; then
                count=$(grep -c '^>' "$file" 2>/dev/null || echo 0)
                echo "✓ $(basename $file) ($count sequences)"
            elif [[ "$file" == *.csv ]]; then
                count=$(($(wc -l < "$file") - 1))  # Subtract header
                echo "✓ $(basename $file) ($count rows)"
            fi
        else
            echo "✗ $(basename $file) missing"
            all_files_present=false
        fi
    done

    if [ "$all_files_present" = true ]; then
        echo ""
        echo "🎉 Chunked protein sequence generation completed successfully!"
        echo ""
        echo "Generated datasets:"
        echo "  └─ full/proteins/                        # Fast generation with pre-validated mutations"
        echo "     ├─ protein_sequences_pairs.*                # Canonical + truncated/extended pairs"
        echo "     └─ protein_sequences_with_mutations.*       # With pre-validated mutations applied"
        echo ""
        echo "Performance benefits:"
        echo "  ├─ ⚡ Parallel processing with 8 chunks"
        echo "  ├─ ⚡ No mutation re-fetching (used cached results)"
        echo "  ├─ ⚡ No mutation re-validation (used pre-validated impacts)"
        echo "  └─ ⚡ Direct mutation application from step 2 results"
        echo ""
        echo "Dataset composition:"
        echo "  ├─ Pairs: canonical + alternative (truncated/extended) proteins"
        echo "  ├─ Mutations: canonical + alternative + mutated variants (fast mode)"
        echo "  └─ Mutations source: ClinVar missense variants only"
        echo ""
        echo "Next step:"
        echo "  Run: sbatch 4_predict_localization.sh"
        
        # Clean up chunk directories
        echo ""
        echo "Cleaning up chunk directories..."
        rm -rf ../results/full/proteins/chunk_*
    else
        echo ""
        echo "❌ Chunked protein generation failed. Some output files are missing."
        echo "Check the generate_proteins.py script and logs for errors."
        exit 1
    fi
fi

# Final cleanup
if [ "$SLURM_ARRAY_TASK_ID" -eq 8 ]; then
    sleep 30
    rm -rf ../results/temp/protein_chunks
fi