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

# Verify downloads
echo ""
echo "Verifying downloaded files..."
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

if [ "$all_present" = true ]; then
    echo ""
    echo "🎉 All required files downloaded successfully!"
    echo ""
    echo "Next steps:"
    echo "1. Place your ribosome profiling BED files in ../ribosome_profiling/"
    echo "2. Add preferred transcripts to hela_top_transcript.txt (optional)"
    echo "3. Run: bash ../scripts/1_cleanup_files.sh"
else
    echo ""
    echo "❌ Some files are missing. Please check the download errors above."
    exit 1
fi

cd ../../scripts