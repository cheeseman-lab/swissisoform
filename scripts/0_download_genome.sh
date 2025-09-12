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
if ls CancerMutationCensus*/*CancerMutationCensus*.tsv 1> /dev/null 2>&1; then
    cosmic_file=$(ls CancerMutationCensus*/*CancerMutationCensus*.tsv | head -n 1)
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
    
    # Step 1: Download Cancer Mutation Census via gget
    echo "================================================"
    echo "Step 1: Downloading Cancer Mutation Census via gget"
    echo "================================================"
    echo ""
    echo "You will be prompted for your COSMIC email and password."
    echo "Register for free at: https://cancer.sanger.ac.uk/cosmic/register"
    echo ""
    
    # Activate conda environment for Python access
    if command -v conda &> /dev/null; then
        echo "Activating conda environment..."
        source ~/.bashrc 2>/dev/null || true
        conda activate swissisoform 2>/dev/null || {
            echo "⚠ Could not activate swissisoform environment"
            echo "  Please ensure the environment exists and gget is installed"
        }
    fi
    
    # Download COSMIC database using Python
    gget_success=false
    python3 -c "
import gget
import sys
try:
    print('Starting COSMIC Cancer Mutation Census download...')
    print('')
    
    # Download the cancer mutation census database
    result = gget.cosmic(
        searchterm=None, 
        download_cosmic=True, 
        cosmic_project='cancer',
        grch_version=37
    )
    print('✓ Cancer Mutation Census downloaded successfully')
    
except Exception as e:
    print(f'❌ Cancer Mutation Census download failed: {str(e)}')
    print('You can continue without COSMIC or try the manual download instructions below.')
    sys.exit(1)
" && gget_success=true
    
    # Step 2: Download NonCoding variants
    echo ""
    echo "================================================"
    echo "Step 2: Download NonCoding Variants (Interactive)"
    echo "================================================"
    echo ""
    echo "To get NonCoding variants, we need your COSMIC credentials and a download URL."
    echo ""
    
    # Get user credentials
    read -p "Enter your COSMIC email: " cosmic_email
    read -s -p "Enter your COSMIC password: " cosmic_password
    echo ""
    
    # Generate auth string
    auth_string=$(echo "${cosmic_email}:${cosmic_password}" | base64)
    echo "Generated authentication string..."
    
    # Get download URL
    echo "Getting download URL from COSMIC API..."
    api_response=$(curl -s -H "Authorization: Basic ${auth_string}" \
        "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v102/Cosmic_NonCodingVariants_Tsv_v102_GRCh38.tar&bucket=downloads")
    
    # Extract URL from JSON response
    download_url=$(echo "$api_response" | sed -n 's/.*"url":"\([^"]*\)".*/\1/p')
    
    if [ -n "$download_url" ] && [[ "$download_url" == http* ]]; then
        echo "Download URL obtained successfully."
        echo ""
        echo "Downloading NonCoding variants..."
        
        # Download the file
        curl -o "Cosmic_NonCodingVariants_v102_GRCh38.tar" "$download_url"
        
        if [ $? -eq 0 ] && [ -f "Cosmic_NonCodingVariants_v102_GRCh38.tar" ]; then
            echo "Download completed. Extracting and organizing files..."
            
            # Extract and organize
            tar -xf Cosmic_NonCodingVariants_v102_GRCh38.tar
            
            # Handle compressed files
            if [ -f "Cosmic_NonCodingVariants_v102_GRCh38.tsv.gz" ]; then
                gunzip Cosmic_NonCodingVariants_v102_GRCh38.tsv.gz
            fi
            
            # Move any remaining NonCodingVariant files to current directory
            find . -maxdepth 2 -name "*NonCodingVariant*" -not -path "./*" -exec mv {} . \; 2>/dev/null || true
            
            # Clean up
            rm -f Cosmic_NonCodingVariants_v102_GRCh38.tar
            
            echo "✓ NonCoding variants downloaded and organized successfully"
        else
            echo "❌ Download failed. Please check your connection and try again."
        fi
    else
        echo "❌ Failed to get download URL. Please check your credentials."
        echo "You can skip NonCoding variants and continue with just the Cancer Mutation Census."
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

if ls ../mutation_data/CancerMutationCensus*/*CancerMutationCensus*.tsv 1> /dev/null 2>&1; then
    cosmic_file=$(ls ../mutation_data/CancerMutationCensus*/*CancerMutationCensus*.tsv | head -n 1)
    size=$(du -h "$cosmic_file" | cut -f1)
    echo "✓ COSMIC Cancer Mutation Census: $(basename $cosmic_file) ($size)"
    cosmic_available=true
fi

if ls ../mutation_data/*Cosmic_NonCodingVariants*.tsv 1> /dev/null 2>&1; then
    noncoding_file=$(ls ../mutation_data/*Cosmic_NonCodingVariants*.tsv | head -n 1)
    size=$(du -h "$noncoding_file" | cut -f1)
    echo "✓ COSMIC NonCoding Variants: $(basename $noncoding_file) ($size)"
    cosmic_available=true
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