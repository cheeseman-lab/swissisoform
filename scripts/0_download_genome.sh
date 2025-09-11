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

# Download COSMIC database
echo ""
echo "Downloading COSMIC database..."
cd ../mutation_data

# Check if COSMIC database already exists
if ls *CancerMutationCensus*.tsv 1> /dev/null 2>&1; then
    echo "✓ COSMIC database already exists"
    cosmic_file=$(ls *CancerMutationCensus*.tsv | head -n 1)
    size=$(du -h "$cosmic_file" | cut -f1)
    echo "  Found: $cosmic_file ($size)"
else
    echo "Downloading COSMIC database (this may take several minutes and requires COSMIC login)..."
    echo "You will be prompted for your COSMIC email and password."
    echo ""
    
    # Activate conda environment for Python access
    if command -v conda &> /dev/null; then
        # Try to activate the swissisoform environment
        echo "Activating conda environment..."
        source ~/.bashrc 2>/dev/null || true
        conda activate swissisoform 2>/dev/null || {
            echo "⚠ Could not activate swissisoform environment"
            echo "  Please ensure the environment exists and gget is installed"
        }
    fi
    
    # Download COSMIC database using Python
    python3 -c "
import gget
import sys
try:
    print('Starting COSMIC database download...')
    print('Note: This requires a COSMIC account (free for academic use)')
    print('Register at: https://cancer.sanger.ac.uk/cosmic/register')
    print('')
    
    # Download the cancer mutation census database
    result = gget.cosmic(
        searchterm=None, 
        download_cosmic=True, 
        cosmic_project='cancer',
        grch_version=37
    )
    print('✓ COSMIC database downloaded successfully')
    
except Exception as e:
    print(f'❌ COSMIC download failed: {str(e)}')
    print('')
    print('Troubleshooting:')
    print('1. Ensure you have a COSMIC account (free for academic use)')
    print('2. Register at: https://cancer.sanger.ac.uk/cosmic/register')
    print('3. Verify your credentials are correct')
    print('4. Check your internet connection')
    print('')
    print('You can also download manually later using:')
    print('  python3 -c \"import gget; gget.cosmic(None, download_cosmic=True, cosmic_project=\\\"cancer\\\")\"')
    sys.exit(1)
"
    
    # Verify download
    if ls CancerMutationCensus*/*CancerMutationCensus*.tsv 1> /dev/null 2>&1; then
        cosmic_file=$(ls CancerMutationCensus*/*CancerMutationCensus*.tsv | head -n 1)
        size=$(du -h "$cosmic_file" | cut -f1)
        echo "✓ Downloaded: $cosmic_file ($size)"
    else
        echo "⚠ COSMIC download appears to have failed"
        echo "  You can download it manually later if needed"
    fi
fi

cd ../genome_data

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

# Check COSMIC database
echo ""
echo "Checking COSMIC database..."
if ls ../mutation_data/CancerMutationCensus*/*CancerMutationCensus*.tsv 1> /dev/null 2>&1; then
    cosmic_file=$(ls ../mutation_data/CancerMutationCensus*/*CancerMutationCensus*.tsv | head -n 1)
    size=$(du -h "$cosmic_file" | cut -f1)
    echo "✓ COSMIC database: $(basename $cosmic_file) ($size)"
    cosmic_available=true
else
    echo "⚠ COSMIC database not available"
    echo "  Pipeline will work without COSMIC, but mutation analysis will be limited"
    echo "  You can download it later using:"
    echo "    cd ../data/mutation_data"
    echo "    python3 -c \"import gget; gget.cosmic(None, download_cosmic=True, cosmic_project='cancer')\""
    cosmic_available=false
fi

if [ "$all_present" = true ]; then
    echo ""
    echo "🎉 Required files downloaded successfully!"
    
    if [ "$cosmic_available" = true ]; then
        echo "📊 COSMIC database is ready for mutation analysis"
    fi
    
    echo ""
    echo "Next steps:"
    echo "1. Place your ribosome profiling BED files in ../ribosome_profiling/"
    echo "2. Add preferred transcripts to hela_top_transcript.txt (optional)"
    
    if [ "$cosmic_available" = true ]; then
        echo "3. Run: bash 1_cleanup_files.sh"
        echo ""
        echo "Available mutation sources: ClinVar, gnomAD, COSMIC"
    else
        echo "3. Download COSMIC database (optional, for comprehensive mutation analysis)"
        echo "4. Run: bash 1_cleanup_files.sh"
        echo ""
        echo "Available mutation sources: ClinVar, gnomAD (COSMIC will be skipped)"
    fi
else
    echo ""
    echo "❌ Some files are missing. Please check the download errors above."
    exit 1
fi

cd ../../scripts