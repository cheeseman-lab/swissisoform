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

# Check if COSMIC MutantCensus database already exists
cosmic_exists=false
if ls *Cosmic_MutantCensus*.tsv 1> /dev/null 2>&1; then
    cosmic_file=$(ls *Cosmic_MutantCensus*.tsv | head -n 1)
    size=$(du -h "$cosmic_file" | cut -f1)
    echo "✓ COSMIC Mutant Census database already exists"
    echo "  Found: $(basename $cosmic_file) ($size)"
    cosmic_exists=true
fi

if [ "$cosmic_exists" = false ]; then
    echo ""
    echo "Downloading COSMIC Mutant Census database (GRCh38)..."
    echo ""
    echo "You will need your COSMIC credentials (register free at: https://cancer.sanger.ac.uk/cosmic/register)"
    echo ""
    
    # Get user credentials
    read -p "Enter your COSMIC email: " cosmic_email
    read -s -p "Enter your COSMIC password: " cosmic_password
    echo ""
    
    # Generate auth string
    auth_string=$(echo "${cosmic_email}:${cosmic_password}" | base64)
    echo "Generated authentication string..."
    
    # Get download URL for Cosmic Mutant Census
    echo ""
    echo "Getting download URL for COSMIC Mutant Census..."
    
    census_response=$(curl -s -H "Authorization: Basic ${auth_string}" \
        "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v102/Cosmic_MutantCensus_Tsv_v102_GRCh38.tar&bucket=downloads")
    
    census_url=$(echo "$census_response" | sed -n 's/.*"url":"\([^"]*\)".*/\1/p')
    
    if [ -n "$census_url" ] && [[ "$census_url" == http* ]]; then
        echo "Downloading COSMIC Mutant Census..."
        curl -o "Cosmic_MutantCensus_Tsv_v102_GRCh38.tar" "$census_url"
        
        if [ $? -eq 0 ] && [ -f "Cosmic_MutantCensus_Tsv_v102_GRCh38.tar" ]; then
            echo "Extracting COSMIC Mutant Census..."
            tar -xf Cosmic_MutantCensus_Tsv_v102_GRCh38.tar
            
            # Find and extract any compressed files
            find . -name "*.tsv.gz" -exec gunzip {} \;
            
            # Clean up tar file
            rm -f Cosmic_MutantCensus_Tsv_v102_GRCh38.tar
            
            echo "✓ COSMIC Mutant Census downloaded and extracted successfully"
            
            # Convert to Parquet for faster querying
            echo ""
            echo "Converting COSMIC file to Parquet format for faster querying..."
            
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

print('Converting COSMIC Mutant Census to Parquet for faster querying...')

# Convert Mutant Census
census_files = list(Path('.').glob('*Cosmic_MutantCensus*.tsv'))
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
        print('✓ Successfully converted COSMIC Mutant Census to Parquet')
        print('  Future gene queries will be much faster!')
    except Exception as e:
        print(f'❌ Failed to convert COSMIC Mutant Census: {e}')
        print('  TSV file will still work but queries will be slower')
else:
    print('❌ No COSMIC Mutant Census TSV file found for conversion')
"
            
        else
            echo "❌ COSMIC Mutant Census download failed"
        fi
    else
        echo "❌ Failed to get COSMIC Mutant Census download URL"
        echo "Response: $census_response"
    fi
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

if ls ../mutation_data/*Cosmic_MutantCensus*.tsv 1> /dev/null 2>&1; then
    cosmic_file=$(ls ../mutation_data/*Cosmic_MutantCensus*.tsv | head -n 1)
    size=$(du -h "$cosmic_file" | cut -f1)
    echo "✓ COSMIC Mutant Census: $(basename $cosmic_file) ($size)"
    cosmic_available=true
fi

if ls ../mutation_data/*Cosmic_MutantCensus*.parquet 1> /dev/null 2>&1; then
    parquet_file=$(ls ../mutation_data/*Cosmic_MutantCensus*.parquet | head -n 1)
    size=$(du -h "$parquet_file" | cut -f1)
    echo "✓ COSMIC Mutant Census (Parquet): $(basename $parquet_file) ($size)"
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
        echo "📊 COSMIC Mutant Census database is ready for mutation analysis"
        echo "🔬 Available mutation sources: ClinVar, gnomAD, COSMIC"
        echo "⚡ Parquet file created for fast gene queries"
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