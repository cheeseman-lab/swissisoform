#!/bin/bash
#
# SwissIsoform Pipeline Step 0: Interactive Genome Data Setup
#
# This script guides you through setting up reference files for your datasets.
# It asks about your BED files, downloads only what you need, and creates
# a configuration file for downstream processing.
#
# Usage:
#   bash scripts/0_download_genome.sh [--force] [--skip-cosmic]
#
# Options:
#   --force        Force re-download even if files exist
#   --skip-cosmic  Skip COSMIC database download
#

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Default settings
FORCE_DOWNLOAD=false
SKIP_COSMIC=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --force)
            FORCE_DOWNLOAD=true
            shift
            ;;
        --skip-cosmic)
            SKIP_COSMIC=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [--force] [--skip-cosmic]"
            echo ""
            echo "Options:"
            echo "  --force        Force re-download even if files exist"
            echo "  --skip-cosmic  Skip COSMIC database download"
            echo "  -h, --help     Show this help message"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   SwissIsoform Pipeline Step 0: Interactive Genome Setup     ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "This script will help you set up reference files for your datasets."
echo "It will ask about your BED files and download only what you need."
echo ""

# Create data directories
echo -e "${YELLOW}→${NC} Setting up data directories..."
mkdir -p ../data/genome_data
mkdir -p ../data/mutation_data/cosmic_data
mkdir -p ../data/ribosome_profiling
echo -e "${GREEN}✓${NC} Directories created"
echo ""

# ============================================================================
# Interactive Dataset Configuration
# ============================================================================

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Dataset Configuration                                       ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Check for existing BED files
echo -e "${YELLOW}→${NC} Scanning for BED files in ../data/ribosome_profiling/..."
bed_files=($(find ../data/ribosome_profiling -maxdepth 1 -name "*.bed" -type f 2>/dev/null | sort))

if [ ${#bed_files[@]} -eq 0 ]; then
    echo -e "${YELLOW}⚠${NC} No BED files found in ../data/ribosome_profiling/"
    echo ""
    echo "Please add your BED files to ../data/ribosome_profiling/ first."
    echo "Then run this script again."
    exit 1
fi

echo -e "${GREEN}✓${NC} Found ${#bed_files[@]} BED file(s):"
for bed_file in "${bed_files[@]}"; do
    basename_file=$(basename "$bed_file")
    echo "  → $basename_file"
done
echo ""

# Ask how many datasets to configure
read -p "How many datasets would you like to configure? [${#bed_files[@]}]: " num_datasets
num_datasets=${num_datasets:-${#bed_files[@]}}

# Validate input
if ! [[ "$num_datasets" =~ ^[0-9]+$ ]] || [ "$num_datasets" -lt 1 ]; then
    echo -e "${RED}✗${NC} Invalid number of datasets"
    exit 1
fi

echo ""
echo -e "${CYAN}Configuring ${num_datasets} dataset(s)...${NC}"
echo ""

# Arrays to store dataset info
declare -a dataset_names
declare -a dataset_bed_files
declare -a dataset_id_types
declare -a dataset_source_gtfs
declare -a needed_gencode_versions

# Configure each dataset
for ((i=1; i<=num_datasets; i++)); do
    echo -e "${BLUE}═══════════════════════════════════════════════════════════${NC}"
    echo -e "${BLUE}Dataset $i of $num_datasets${NC}"
    echo -e "${BLUE}═══════════════════════════════════════════════════════════${NC}"
    echo ""

    # Dataset name
    read -p "Dataset name (e.g., ribosome_profiling_2015): " dataset_name
    while [ -z "$dataset_name" ]; do
        echo -e "${YELLOW}⚠ Dataset name cannot be empty${NC}"
        read -p "Dataset name: " dataset_name
    done
    dataset_names+=("$dataset_name")

    # BED file selection
    echo ""
    echo "Available BED files:"
    for idx in "${!bed_files[@]}"; do
        basename_file=$(basename "${bed_files[$idx]}")
        echo "  $((idx+1)). $basename_file"
    done
    read -p "Select BED file [1-${#bed_files[@]}]: " bed_selection

    # Validate BED selection
    if ! [[ "$bed_selection" =~ ^[0-9]+$ ]] || [ "$bed_selection" -lt 1 ] || [ "$bed_selection" -gt ${#bed_files[@]} ]; then
        echo -e "${RED}✗${NC} Invalid selection, using first file"
        bed_selection=1
    fi

    selected_bed=$(basename "${bed_files[$((bed_selection-1))]}")
    dataset_bed_files+=("$selected_bed")
    echo -e "${GREEN}✓${NC} Selected: $selected_bed"

    # All datasets use ENST (Ensembl/GENCODE)
    dataset_id_types+=("ensembl")
    echo ""
    echo -e "${CYAN}Transcript IDs: ENST (Ensembl/GENCODE) - e.g., ENST00000344517.4${NC}"

    # Ask for GENCODE version
    echo ""
    echo "Source GENCODE version:"
    echo "  1. GENCODE v24 (2015, GRCh38.p5)"
    echo "  2. GENCODE v25 (2015, GRCh38.p7)"
    read -p "Select version [1-2]: " gencode_selection

    case $gencode_selection in
        1)
            dataset_source_gtfs+=("gencode_v24")
            needed_gencode_versions+=("v24")
            echo -e "${GREEN}✓${NC} Source: GENCODE v24"
            ;;
        2)
            dataset_source_gtfs+=("gencode_v25")
            needed_gencode_versions+=("v25")
            echo -e "${GREEN}✓${NC} Source: GENCODE v25"
            ;;
        *)
            echo -e "${YELLOW}⚠ Invalid selection, defaulting to v25${NC}"
            dataset_source_gtfs+=("gencode_v25")
            needed_gencode_versions+=("v25")
            ;;
    esac

    echo ""
done

# ============================================================================
# Generate Dataset Configuration YAML
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Generating Configuration File                               ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

config_file="../data/ribosome_profiling/dataset_config.yaml"

cat > "$config_file" << 'EOF'
# SwissIsoform Dataset Configuration
# Generated by 0_download_genome.sh
#
# This file describes your datasets and reference files.
# Used by downstream scripts for standardization and cleanup.

genome:
  reference: "GRCh38.p7.genome.fa"
  reference_path: "../genome_data/GRCh38.p7.genome.fa"
  target_standard: "gencode_v25"
  target_gtf_path: "../genome_data/gencode.v25.annotation.gtf"

datasets:
EOF

# Add each dataset to config
for ((i=0; i<num_datasets; i++)); do
    version_num=$(echo "${dataset_source_gtfs[$i]}" | sed 's/gencode_v//')
    cat >> "$config_file" << EOF
  - name: "${dataset_names[$i]}"
    bed_file: "${dataset_bed_files[$i]}"
    id_type: "${dataset_id_types[$i]}"
    source_gtf: "${dataset_source_gtfs[$i]}"
    source_gtf_path: "../genome_data/gencode.v${version_num}.annotation.gtf"

EOF
done

echo -e "${GREEN}✓${NC} Configuration saved to: $config_file"
echo ""

# Show configuration summary
echo -e "${CYAN}Configuration Summary:${NC}"
echo "  Target standard: GENCODE v25"
echo "  Datasets configured: $num_datasets"
for ((i=0; i<num_datasets; i++)); do
    echo "    $((i+1)). ${dataset_names[$i]}"
    echo "       BED: ${dataset_bed_files[$i]}"
    echo "       Type: ${dataset_id_types[$i]} (ENST)"
    echo "       Source: ${dataset_source_gtfs[$i]}"
done
echo ""

# ============================================================================
# Download Reference Files
# ============================================================================

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Downloading Reference Files                                 ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

cd ../data/genome_data

# Download reference genome
echo -e "${BLUE}Reference Genome (GRCh38.p7)${NC}"
echo "────────────────────────────────────────────────────────────────"
if [ -f "GRCh38.p7.genome.fa" ] && [ "$FORCE_DOWNLOAD" = false ]; then
    SIZE=$(du -h "GRCh38.p7.genome.fa" | cut -f1)
    echo -e "${GREEN}✓${NC} GRCh38.p7.genome.fa already exists (${SIZE})"
    echo -e "${CYAN}  Skipping download...${NC}"
else
    echo -e "${YELLOW}→${NC} Downloading reference genome..."
    if wget -q --show-progress -O GRCh38.p7.genome.fa.gz \
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/GRCh38.p7.genome.fa.gz"; then
        echo -e "${YELLOW}→${NC} Extracting..."
        gunzip GRCh38.p7.genome.fa.gz
        SIZE=$(du -h "GRCh38.p7.genome.fa" | cut -f1)
        echo -e "${GREEN}✓${NC} Downloaded and extracted (${SIZE})"
    else
        echo -e "${RED}✗ Download failed${NC}"
        exit 1
    fi
fi
echo ""

# Download GENCODE versions (v24 and v25)
unique_versions=($(printf '%s\n' "${needed_gencode_versions[@]}" | sort -u))

for version in "${unique_versions[@]}"; do
    version_num=$(echo "$version" | sed 's/v//')

    echo -e "${BLUE}GENCODE $version (Source Version)${NC}"
    echo "────────────────────────────────────────────────────────────────"

    if [ -f "gencode.v${version_num}.annotation.gtf" ] && [ "$FORCE_DOWNLOAD" = false ]; then
        SIZE=$(du -h "gencode.v${version_num}.annotation.gtf" | cut -f1)
        echo -e "${GREEN}✓${NC} gencode.v${version_num}.annotation.gtf already exists (${SIZE})"
        echo -e "${CYAN}  Skipping download...${NC}"
    else
        echo -e "${YELLOW}→${NC} Downloading GENCODE v${version_num}..."
        if wget -q --show-progress -O gencode.v${version_num}.annotation.gtf.gz \
            "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${version_num}/gencode.v${version_num}.annotation.gtf.gz"; then
            echo -e "${YELLOW}→${NC} Extracting..."
            gunzip gencode.v${version_num}.annotation.gtf.gz
            SIZE=$(du -h "gencode.v${version_num}.annotation.gtf" | cut -f1)
            echo -e "${GREEN}✓${NC} Downloaded and extracted (${SIZE})"
        else
            echo -e "${RED}✗ Download failed${NC}"
            exit 1
        fi
    fi
    echo ""
done

# Download COSMIC (optional)
if [ "$SKIP_COSMIC" = false ]; then
    echo -e "${BLUE}COSMIC Database (Optional)${NC}"
    echo "────────────────────────────────────────────────────────────────"

    cd ../mutation_data

    if [ -f "cosmic_variants_combined.parquet" ] && [ "$FORCE_DOWNLOAD" = false ]; then
        SIZE=$(du -h "cosmic_variants_combined.parquet" | cut -f1)
        echo -e "${GREEN}✓${NC} COSMIC database already exists (${SIZE})"
        echo -e "${CYAN}  Skipping download...${NC}"
    else
        cd cosmic_data
        echo ""
        echo "COSMIC database requires authentication."
        echo "Register for free at: https://cancer.sanger.ac.uk/cosmic/register"
        echo ""
        read -p "Do you want to download COSMIC data? (y/N): " download_cosmic

        if [[ "$download_cosmic" =~ ^[Yy]$ ]]; then
            echo ""
            read -p "Enter your COSMIC email: " cosmic_email
            read -s -p "Enter your COSMIC password: " cosmic_password
            echo ""
            echo ""

            # Generate auth string
            auth_string=$(echo -n "${cosmic_email}:${cosmic_password}" | base64)

            # COSMIC v102 GRCh38 VCF tar files to download
            cosmic_files=(
                "Cosmic_GenomeScreensMutant_Vcf_v102_GRCh38.tar"
                "Cosmic_NonCodingVariants_Vcf_v102_GRCh38.tar"
                "Cosmic_CompleteTargetedScreensMutant_Vcf_v102_GRCh38.tar"
            )

            cosmic_paths=(
                "grch38/cosmic/v102/VCF/Cosmic_GenomeScreensMutant_Vcf_v102_GRCh38.tar"
                "grch38/cosmic/v102/VCF/Cosmic_NonCodingVariants_Vcf_v102_GRCh38.tar"
                "grch38/cosmic/v102/VCF/Cosmic_CompleteTargetedScreensMutant_Vcf_v102_GRCh38.tar"
            )

            cosmic_descriptions=(
                "Genome Screens Mutant (VCF)"
                "Non-Coding Variants (VCF)"
                "Complete Targeted Screens Mutant (VCF)"
            )

            # Download each COSMIC tar file
            all_downloaded=true
            for idx in "${!cosmic_files[@]}"; do
                file="${cosmic_files[$idx]}"
                path="${cosmic_paths[$idx]}"
                desc="${cosmic_descriptions[$idx]}"

                echo -e "${BLUE}${desc}${NC}"
                echo "────────────────────────────────────────────────────────────────"

                if [ -f "$file" ] && [ "$FORCE_DOWNLOAD" = false ]; then
                    SIZE=$(du -h "$file" | cut -f1)
                    echo -e "${GREEN}✓${NC} $file already exists (${SIZE})"
                    echo -e "${CYAN}  Skipping download...${NC}"
                else
                    echo -e "${YELLOW}→${NC} Downloading $file..."

                    # Get download URL from COSMIC API
                    url_response=$(curl -s -H "Authorization: Basic ${auth_string}" \
                        "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=${path}&bucket=downloads")

                    download_url=$(echo "$url_response" | grep -o '"url":"[^"]*"' | sed 's/"url":"//;s/"$//' | sed 's/\\//g')

                    if [ -z "$download_url" ]; then
                        echo -e "${RED}✗${NC} Failed to get download URL"
                        echo -e "${YELLOW}  Response: $url_response${NC}"
                        echo -e "${YELLOW}  Check your COSMIC email/password${NC}"
                        all_downloaded=false
                        break
                    fi

                    # Download the file
                    if wget -q --show-progress -O "$file" "$download_url"; then
                        SIZE=$(du -h "$file" | cut -f1)
                        echo -e "${GREEN}✓${NC} Downloaded $file (${SIZE})"
                    else
                        echo -e "${RED}✗${NC} Download failed"
                        all_downloaded=false
                        break
                    fi
                fi
                echo ""
            done

            # Extract VCF files from tar archives
            if [ "$all_downloaded" = true ]; then
                echo -e "${BLUE}Extracting VCF Files${NC}"
                echo "────────────────────────────────────────────────────────────────"
                for file in "${cosmic_files[@]}"; do
                    if [ -f "$file" ]; then
                        echo -e "${YELLOW}→${NC} Extracting $file..."
                        tar -xf "$file"
                        echo -e "${GREEN}✓${NC} Extracted VCF files from $file"
                    fi
                done
                echo ""
            fi

            # Convert VCF files to Parquet if all downloaded successfully
            if [ "$all_downloaded" = true ]; then
                echo -e "${BLUE}Converting VCF to Parquet Format${NC}"
                echo "────────────────────────────────────────────────────────────────"
                echo -e "${YELLOW}→${NC} Converting COSMIC VCF files to Parquet..."
                echo ""

                # Check if Python and required packages are available
                if ! python3 -c "import pandas, pyarrow" 2>/dev/null; then
                    echo -e "${YELLOW}⚠${NC} Python packages required: pandas, pyarrow"
                    echo -e "${CYAN}  Install with: pip install pandas pyarrow${NC}"
                    echo -e "${YELLOW}  Skipping Parquet conversion...${NC}"
                else
                    # Create Python script for VCF to Parquet conversion with INFO field parsing
                    python3 << 'PYTHON_SCRIPT'
import pandas as pd
import gzip
import os
import shutil
import time

def parse_info_field(info_string):
    """Parse VCF INFO field into a dictionary."""
    info_dict = {}
    for item in info_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True  # Flag fields without values
    return info_dict

def process_vcf_to_structured_df(vcf_file):
    """Read VCF and parse INFO field into separate columns."""
    vcf_size = os.path.getsize(vcf_file) / (1024**2)
    print(f"  → Parsing {vcf_file} ({vcf_size:.1f}MB compressed)...")
    print(f"    Reading and expanding INFO fields (this may take a few minutes)...")

    start_time = time.time()

    all_rows = []
    total_rows = 0
    chunk_size = 50000

    with gzip.open(vcf_file, 'rt') as f:
        # Skip header lines starting with ##
        for line in f:
            if line.startswith('#CHROM'):
                break

        # Process data lines
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom, pos, variant_id, ref, alt, qual, filt, info = fields[:8]

            # Parse INFO field
            info_dict = parse_info_field(info)

            # Create row with expanded columns matching mutations.py expectations
            row = {
                # Standard VCF fields
                'CHROMOSOME': chrom,
                'GENOME_START': int(pos),
                'GENOMIC_MUTATION_ID': variant_id,
                'GENOMIC_WT_ALLELE': ref,
                'GENOMIC_MUT_ALLELE': alt,

                # INFO field expansions (map to mutations.py expected names)
                'GENE_SYMBOL': info_dict.get('GENE', ''),
                'TRANSCRIPT_ACCESSION': info_dict.get('TRANSCRIPT', ''),
                'MUTATION_CDS': info_dict.get('CDS', ''),
                'MUTATION_AA': info_dict.get('AA', ''),
                'MUTATION_DESCRIPTION': info_dict.get('SO_TERM', ''),
                'HGVSG': info_dict.get('HGVSG', ''),
                'HGVSC': info_dict.get('HGVSC', ''),
                'HGVSP': info_dict.get('HGVSP', ''),
                'STRAND': info_dict.get('STRAND', ''),
                'LEGACY_MUTATION_ID': info_dict.get('LEGACY_ID', ''),

                # VCF-specific fields not in TSV format
                'GENOME_SCREEN_SAMPLE_COUNT': info_dict.get('GENOME_SCREEN_SAMPLE_COUNT', ''),
                'IS_CANONICAL': info_dict.get('IS_CANONICAL', ''),

                # Fields missing in VCF but expected by mutations.py (set to None)
                'COSMIC_SAMPLE_ID': None,
                'COSMIC_STUDY_ID': None,
                'MUTATION_SOMATIC_STATUS': None,
                'MUTATION_ZYGOSITY': None,
                'PUBMED_PMID': None,
                'SAMPLE_NAME': None,
            }

            all_rows.append(row)
            total_rows += 1

            # Print progress every 100k rows
            if total_rows % 100000 == 0:
                elapsed = time.time() - start_time
                print(f"    ... processed {total_rows:,} rows ({elapsed:.1f}s elapsed)")

    df = pd.DataFrame(all_rows)
    elapsed = time.time() - start_time

    print(f"    ✓ Parsed {total_rows:,} variants in {elapsed:.1f}s")
    print(f"    ✓ Expanded INFO into {len(df.columns)} columns")

    return df, vcf_size

# VCF filenames for COSMIC v102 (after extraction from tar)
vcf_files = [
    "Cosmic_GenomeScreensMutant_v102_GRCh38.vcf.gz",
    "Cosmic_NonCodingVariants_v102_GRCh38.vcf.gz",
    "Cosmic_CompleteTargetedScreensMutant_v102_GRCh38.vcf.gz"
]

parquet_files = []

for vcf_file in vcf_files:
    if not os.path.exists(vcf_file):
        print(f"  ⚠ {vcf_file} not found, skipping...")
        continue

    # Process VCF to structured DataFrame
    df, vcf_size = process_vcf_to_structured_df(vcf_file)

    # Write to parquet
    parquet_file = vcf_file.replace('.vcf.gz', '.parquet')
    print(f"    Writing to parquet...")
    df.to_parquet(parquet_file, compression='snappy', engine='pyarrow')

    parquet_size = os.path.getsize(parquet_file) / (1024**2)
    compression_ratio = (parquet_size / vcf_size) * 100

    print(f"    ✓ Created {parquet_file}")
    print(f"      Size: {vcf_size:.1f}MB → {parquet_size:.1f}MB ({compression_ratio:.1f}%)")
    print()

    parquet_files.append(parquet_file)

if parquet_files:
    print("  → Creating combined COSMIC variants file...")

    # Read all parquet files and combine
    all_data = [pd.read_parquet(pf) for pf in parquet_files]
    combined_df = pd.concat(all_data, ignore_index=True)

    print(f"    Sorting by chromosome and position...")
    combined_df = combined_df.sort_values(['CHROMOSOME', 'GENOME_START'])

    # Save to parent directory (mutation_data)
    combined_df.to_parquet('../cosmic_variants_combined.parquet', compression='snappy')

    total_size = os.path.getsize('../cosmic_variants_combined.parquet') / (1024**2)
    print(f"    ✓ Created ../cosmic_variants_combined.parquet")
    print(f"      → Total variants: {len(combined_df):,}")
    print(f"      → File size: {total_size:.1f}MB")
    print(f"      → Columns: {len(combined_df.columns)}")
    print(f"      → Sources: {len(parquet_files)} VCF files")
    print()

    # Clean up cosmic_data subfolder if successful
    print("  → Cleaning up cosmic_data subfolder...")
    os.chdir('..')
    if os.path.exists('cosmic_variants_combined.parquet'):
        shutil.rmtree('cosmic_data')
        print("    ✓ Removed cosmic_data subfolder")
    else:
        print("    ⚠ Combined file not created, keeping cosmic_data subfolder")

print("✓ Successfully converted COSMIC VCF files to Parquet with expanded columns")
PYTHON_SCRIPT
                fi
            else
                echo -e "${YELLOW}⚠${NC} Download incomplete, skipping Parquet conversion"
            fi
        else
            echo -e "${YELLOW}⚠${NC} Skipping COSMIC download"
            echo -e "${CYAN}  Pipeline will use ClinVar/gnomAD only${NC}"
        fi

        # Return to mutation_data directory (in case we were in cosmic_data)
        cd /lab/barcheese01/mdiberna/swissisoform/data/mutation_data 2>/dev/null || cd ../
    fi

    cd ../genome_data
else
    echo -e "${YELLOW}⚠${NC} Skipping COSMIC download (--skip-cosmic flag)"
fi

echo ""

# ============================================================================
# Verification
# ============================================================================

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Verification                                                ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

all_present=true

# Check genome
if [ -f "GRCh38.p7.genome.fa" ]; then
    SIZE=$(du -h "GRCh38.p7.genome.fa" | cut -f1)
    echo -e "${GREEN}✓${NC} GRCh38.p7.genome.fa (${SIZE})"
else
    echo -e "${RED}✗${NC} GRCh38.p7.genome.fa missing"
    all_present=false
fi

# Check source GENCODE versions
for version in "${unique_versions[@]}"; do
    version_num=$(echo "$version" | sed 's/v//')
    if [ -f "gencode.v${version_num}.annotation.gtf" ]; then
        SIZE=$(du -h "gencode.v${version_num}.annotation.gtf" | cut -f1)
        echo -e "${GREEN}✓${NC} gencode.v${version_num}.annotation.gtf (${SIZE})"
    else
        echo -e "${RED}✗${NC} gencode.v${version_num}.annotation.gtf missing"
        all_present=false
    fi
done

echo ""

# ============================================================================
# Summary
# ============================================================================

echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Setup Summary                                               ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""

if [ "$all_present" = true ]; then
    echo -e "${GREEN}✓ All required files are ready!${NC}"
    echo ""
    echo "Configuration:"
    echo -e "  ${GREEN}✓${NC} Dataset config: ../data/ribosome_profiling/dataset_config.yaml"
    echo -e "  ${GREEN}✓${NC} Datasets configured: $num_datasets"
    echo -e "  ${GREEN}✓${NC} All transcripts use ENST (Ensembl/GENCODE)"
    echo -e "  ${GREEN}✓${NC} GENCODE versions: v24 and/or v25"
    echo ""
    echo -e "${BLUE}Next steps:${NC}"
    echo "  1. Review config: cat ../data/ribosome_profiling/dataset_config.yaml"
    echo "  2. Run cleanup: bash scripts/1_cleanup_files.sh"
    echo ""
else
    echo -e "${RED}✗ Some required files are missing${NC}"
    echo -e "${YELLOW}  Check the verification section above${NC}"
    exit 1
fi

cd ../../scripts
