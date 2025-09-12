#!/bin/bash
# 0_download_genome.sh
# Downloads required reference genome and annotation files

set -e  # Exit on any error

echo "================================================"
echo "SwissIsoform Pipeline Step 0: Download Genome Data"
echo "================================================"

# Create data directories
echo "Creating data directories..."
mkdir -p ../data/genome_data
mkdir -p ../data/mutation_data
cd ../data/genome_data

# Download reference genome (GRCh38.p7)
echo ""
echo "Downloading reference genome (GRCh38.p7)..."
if [ ! -f "GRCh38.p7.genome.fa" ]; then
    wget -O GRCh38.p7.genome.fa.gz \
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/GRCh38.p7.genome.fa.gz"
    gunzip GRCh38.p7.genome.fa.gz
    echo "✓ Downloaded and extracted GRCh38.p7.genome.fa"
else
    echo "✓ GRCh38.p7.genome.fa already exists"
fi

# Download GENCODE v25 annotations (primary)
echo ""
echo "Downloading GENCODE v25 annotations..."
if [ ! -f "gencode.v25.annotation.gtf" ]; then
    wget -O gencode.v25.annotation.gtf.gz \
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz"
    gunzip gencode.v25.annotation.gtf.gz
    echo "✓ Downloaded and extracted gencode.v25.annotation.gtf"
else
    echo "✓ gencode.v25.annotation.gtf already exists"
fi

# Download GENCODE v47 annotations (for gene name updates)
echo ""
echo "Downloading GENCODE v47 annotations (for gene name updates)..."
if [ ! -f "gencode.v47.annotation.gtf" ]; then
    wget -O gencode.v47.annotation.gtf.gz \
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz"
    gunzip gencode.v47.annotation.gtf.gz
    echo "✓ Downloaded and extracted gencode.v47.annotation.gtf"
else
    echo "✓ gencode.v47.annotation.gtf already exists"
fi

# Download HeLa preferred transcripts (if available)
echo ""
echo "Checking for HeLa preferred transcripts..."
if [ ! -f "hela_top_transcript.txt" ]; then
    echo "⚠ hela_top_transcript.txt not found"
    echo "  Please place your preferred transcript file here manually"
    echo "  Format: one transcript ID per line (e.g., ENST00000000000.1)"
    # Create empty file as placeholder
    touch hela_top_transcript.txt
else
    echo "✓ hela_top_transcript.txt already exists"
fi

# COSMIC database download
echo ""
echo "================================================"
echo "COSMIC Database Download"
echo "================================================"
cd ../mutation_data

# Check if any COSMIC database already exists
cosmic_exists=false
if ls *CancerMutationCensus_AllData*.tsv 1> /dev/null 2>&1; then
    cosmic_file=$(ls *CancerMutationCensus_AllData*.tsv | head -n 1)
    size=$(du -h "$cosmic_file" | cut -f1)
    echo "✓ COSMIC Cancer Mutation Census database already exists"
    echo "  Found: $(basename $cosmic_file) ($size)"
    cosmic_exists=true
fi

if ls *Cosmic_NonCodingVariants*.tsv 1> /dev/null 2>&1; then
    noncoding_file=$(ls *Cosmic_NonCodingVariants*.tsv | head -n 1)
    size=$(du -h "$noncoding_file" | cut -f1)
    echo "✓ COSMIC NonCoding Variants database already exists"
    echo "  Found: $(basename $noncoding_file) ($size)"
    cosmic_exists=true
fi

if [ "$cosmic_exists" = false ]; then
    echo ""
    echo "Downloading COSMIC databases (Cancer Mutation Census + NonCoding variants)..."
    echo ""
    echo "You will need your COSMIC credentials (register free at: https://cancer.sanger.ac.uk/cosmic/register)"
    echo ""
    
    # Get user credentials once for both downloads
    read -p "Enter your COSMIC email: " cosmic_email
    read -s -p "Enter your COSMIC password: " cosmic_password
    echo ""
    
    # Generate auth string
    auth_string=$(echo "${cosmic_email}:${cosmic_password}" | base64)
    echo "Generated authentication string..."
    
    # Step 1: Download Cancer Mutation Census (All Data)
    echo ""
    echo "================================================"
    echo "Step 1: Downloading Cancer Mutation Census (All Data)"
    echo "================================================"
    echo ""
    echo "Getting download URL for Cancer Mutation Census..."
    
    census_response=$(curl -s -H "Authorization: Basic ${auth_string}" \
        "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch37/cmc/v102/CancerMutationCensus_AllData_Tsv_v102_GRCh37.tar&bucket=downloads")
    
    census_url=$(echo "$census_response" | sed -n 's/.*"url":"\([^"]*\)".*/\1/p')
    
    if [ -n "$census_url" ] && [[ "$census_url" == http* ]]; then
        echo "Downloading Cancer Mutation Census..."
        curl -o "CancerMutationCensus_AllData_Tsv_v102_GRCh37.tar" "$census_url"
        
        if [ $? -eq 0 ] && [ -f "CancerMutationCensus_AllData_Tsv_v102_GRCh37.tar" ]; then
            echo "Extracting Cancer Mutation Census..."
            tar -xf CancerMutationCensus_AllData_Tsv_v102_GRCh37.tar
            
            # Find and extract any compressed files
            find . -name "*.tsv.gz" -exec gunzip {} \;
            
            # Clean up tar file
            rm -f CancerMutationCensus_AllData_Tsv_v102_GRCh37.tar
            
            echo "✓ Cancer Mutation Census downloaded and extracted successfully"
        else
            echo "❌ Cancer Mutation Census download failed"
        fi
    else
        echo "❌ Failed to get Cancer Mutation Census download URL"
    fi
    
    # Step 2: Download NonCoding variants
    echo ""
    echo "================================================"
    echo "Step 2: Downloading NonCoding Variants"
    echo "================================================"
    echo ""
    echo "Getting download URL for NonCoding variants..."
    
    noncoding_response=$(curl -s -H "Authorization: Basic ${auth_string}" \
        "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v102/Cosmic_NonCodingVariants_Tsv_v102_GRCh38.tar&bucket=downloads")
    
    noncoding_url=$(echo "$noncoding_response" | sed -n 's/.*"url":"\([^"]*\)".*/\1/p')
    
    if [ -n "$noncoding_url" ] && [[ "$noncoding_url" == http* ]]; then
        echo "Downloading NonCoding variants..."
        curl -o "Cosmic_NonCodingVariants_v102_GRCh38.tar" "$noncoding_url"
        
        if [ $? -eq 0 ] && [ -f "Cosmic_NonCodingVariants_v102_GRCh38.tar" ]; then
            echo "Extracting NonCoding variants..."
            tar -xf Cosmic_NonCodingVariants_v102_GRCh38.tar
            
            # Handle compressed files
            if [ -f "Cosmic_NonCodingVariants_v102_GRCh38.tsv.gz" ]; then
                gunzip Cosmic_NonCodingVariants_v102_GRCh38.tsv.gz
            fi
            
            # Move any NonCoding files to current directory
            find . -maxdepth 2 -name "*NonCodingVariant*" -type f -exec mv {} . \; 2>/dev/null || true
            
            # Clean up
            rm -f Cosmic_NonCodingVariants_v102_GRCh38.tar
            
            echo "✓ NonCoding variants downloaded and extracted successfully"
        else
            echo "❌ NonCoding variants download failed"
        fi
    else
        echo "❌ Failed to get NonCoding variants download URL"
    fi
    
    # Step 3: Convert to Parquet for faster querying
    echo ""
    echo "================================================"
    echo "Step 3: Converting COSMIC files to Parquet format"
    echo "================================================"
    echo ""
    
    # Activate conda environment for Python access
    if command -v conda &> /dev/null; then
        echo "Activating conda environment..."
        source ~/.bashrc 2>/dev/null || true
        conda activate swissisoform 2>/dev/null || {
            echo "⚠ Could not activate swissisoform environment"
            echo "  Parquet conversion may fail if pandas/pyarrow not installed"
        }
    fi
    
    python3 -c "
import pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
import sys

print('Converting COSMIC files to Parquet for faster querying...')

converted_files = []

# Convert Cancer Mutation Census
census_files = list(Path('.').glob('*CancerMutationCensus_AllData*.tsv'))
if census_files:
    census_file = census_files[0]
    print(f'Converting {census_file} to Parquet...')
    try:
        df = pd.read_csv(census_file, sep='\t', low_memory=False)
        parquet_path = str(census_file).replace('.tsv', '.parquet')
        df.to_parquet(parquet_path, compression='snappy', index=False)
        
        original_size = census_file.stat().st_size / 1024**2
        parquet_size = Path(parquet_path).stat().st_size / 1024**2
        print(f'✓ Created {parquet_path}')
        print(f'  Size: {original_size:.1f}MB → {parquet_size:.1f}MB ({parquet_size/original_size*100:.1f}%)')
        converted_files.append('Cancer Mutation Census')
    except Exception as e:
        print(f'❌ Failed to convert Cancer Mutation Census: {e}')

# Convert NonCoding variants
noncoding_files = list(Path('.').glob('*Cosmic_NonCodingVariants*.tsv'))
if noncoding_files:
    noncoding_file = noncoding_files[0]
    print(f'Converting {noncoding_file} to Parquet...')
    try:
        df = pd.read_csv(noncoding_file, sep='\t', low_memory=False)
        parquet_path = str(noncoding_file).replace('.tsv', '.parquet')
        df.to_parquet(parquet_path, compression='snappy', index=False)
        
        original_size = noncoding_file.stat().st_size / 1024**2
        parquet_size = Path(parquet_path).stat().st_size / 1024**2
        print(f'✓ Created {parquet_path}')
        print(f'  Size: {original_size:.1f}MB → {parquet_size:.1f}MB ({parquet_size/original_size*100:.1f}%)')
        converted_files.append('NonCoding variants')
    except Exception as e:
        print(f'❌ Failed to convert NonCoding variants: {e}')

if converted_files:
    print(f'\\n✓ Successfully converted {len(converted_files)} COSMIC datasets to Parquet')
    print('  Future gene queries will be much faster!')
else:
    print('\\n⚠ No COSMIC files were converted to Parquet')
    print('  TSV files will still work but queries will be slower')
"
    
fi

cd ../genome_data

# Verify downloads
echo ""
echo "================================================"
echo "Verifying Downloaded Files"
echo "================================================"

required_files=(
    "GRCh38.p7.genome.fa"
    "gencode.v25.annotation.gtf"
    "gencode.v47.annotation.gtf"
)

all_present=true
for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        size=$(du -h "$file" | cut -f1)
        echo "✓ $file ($size)"
    else
        echo "✗ $file missing"
        all_present=false
    fi
done

# Check COSMIC database
echo ""
echo "Checking COSMIC database..."
cosmic_available=false

if ls ../mutation_data/*CancerMutationCensus_AllData*.tsv 1> /dev/null 2>&1; then
    cosmic_file=$(ls ../mutation_data/*CancerMutationCensus_AllData*.tsv | head -n 1)
    size=$(du -h "$cosmic_file" | cut -f1)
    echo "✓ COSMIC Cancer Mutation Census: $(basename $cosmic_file) ($size)"
    cosmic_available=true
fi

if ls ../mutation_data/*CancerMutationCensus_AllData*.parquet 1> /dev/null 2>&1; then
    parquet_file=$(ls ../mutation_data/*CancerMutationCensus_AllData*.parquet | head -n 1)
    size=$(du -h "$parquet_file" | cut -f1)
    echo "✓ COSMIC Cancer Mutation Census (Parquet): $(basename $parquet_file) ($size)"
fi

if ls ../mutation_data/*Cosmic_NonCodingVariants*.tsv 1> /dev/null 2>&1; then
    noncoding_file=$(ls ../mutation_data/*Cosmic_NonCodingVariants*.tsv | head -n 1)
    size=$(du -h "$noncoding_file" | cut -f1)
    echo "✓ COSMIC NonCoding Variants: $(basename $noncoding_file) ($size)"
    cosmic_available=true
fi

if ls ../mutation_data/*Cosmic_NonCodingVariants*.parquet 1> /dev/null 2>&1; then
    parquet_file=$(ls ../mutation_data/*Cosmic_NonCodingVariants*.parquet | head -n 1)
    size=$(du -h "$parquet_file" | cut -f1)
    echo "✓ COSMIC NonCoding Variants (Parquet): $(basename $parquet_file) ($size)"
fi

if [ "$cosmic_available" = false ]; then
    echo "⚠ No COSMIC database found"
    echo "  Pipeline will work without COSMIC, but mutation analysis will be limited"
    echo "  Available sources: ClinVar, gnomAD"
fi

# Final summary
echo ""
echo "================================================"
echo "Download Summary"
echo "================================================"

if [ "$all_present" = true ]; then
    echo "🎉 Required genome files downloaded successfully!"
    
    if [ "$cosmic_available" = true ]; then
        echo "📊 COSMIC database is ready for mutation analysis"
        echo "🔬 Available mutation sources: ClinVar, gnomAD, COSMIC"
        echo "⚡ Parquet files created for fast gene queries"
    else
        echo "📊 Available mutation sources: ClinVar, gnomAD (COSMIC not available)"
    fi
    
    echo ""
    echo "Next steps:"
    echo "1. Place your ribosome profiling BED files in ../ribosome_profiling/"
    echo "2. Add preferred transcripts to hela_top_transcript.txt (optional)"
    echo "3. Run: bash 1_cleanup_files.sh"
    
else
    echo "❌ Some required files are missing. Please check the download errors above."
    exit 1
fi

cd ../../scripts