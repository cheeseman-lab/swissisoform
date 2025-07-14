# SwissIsoform

SwissIsoform is a Python package for analyzing alternative protein isoforms discovered through ribosome profiling data. It integrates genomic annotations, truncation sites, and mutation data to generate comprehensive protein sequence datasets for downstream analysis and machine learning applications.

## Features

- **Genome sequence and annotation handling** with GENCODE support
- **Alternative truncation site processing** from ribosome profiling BED files
- **Protein sequence generation** from canonical and truncated transcripts
- **Mutation integration** from ClinVar and gnomAD databases
- **Subcellular localization prediction** using DeepLoc
- **Batch processing** for large-scale dataset generation
- **Automated pipeline** with numbered workflow steps

## Quick Start

### Prerequisites
- Linux/Unix system with SLURM workload manager
- Conda package manager
- Access to GPU partition (for localization prediction)

### Installation

1. **Clone the repository:**
```bash
git clone https://github.com/cheeseman-lab/swissisoform.git
cd swissisoform
```

2. **Create and activate the conda environment:**
```bash
conda env create --file=environment.yml
conda activate swissisoform
```

3. **Install DeepLoc (for localization prediction):**
```bash
conda create -n deeploc python=3.8
conda activate deeploc
pip install DeepLoc-2.1.0.tar.gz
```

### Complete Pipeline

```bash
# Step 0: Download reference data
bash 0_download_genome.sh

# Step 1: Clean ribosome profiling data  
bash 1_cleanup_files.sh

# Step 3: Generate protein sequences (SLURM job)
sbatch 3_generate_proteins.sh

# Step 4: Predict subcellular localization (SLURM job)
sbatch 4_predict_localization.sh
```

## Pipeline Overview

### Step 0: Download Reference Data
**Script:** `0_download_genome.sh`  
**Runtime:** 10-30 minutes  
**Downloads:**
- Human reference genome (GRCh38.p7)
- GENCODE v25 annotations (primary)
- GENCODE v47 annotations (for gene name updates)

### Step 1: Data Cleanup
**Script:** `1_cleanup_files.sh`  
**Runtime:** 5-15 minutes  
**Requirements:** Place your ribosome profiling BED files in `data/ribosome_profiling/`
- `full_truncations_JL.bed` - All detected truncation sites
- `selected_truncations_JL.bed` - Curated truncation sites

**Outputs:**
- Cleaned GTF annotations with updated gene names
- Cleaned BED files with validated coordinates
- Gene lists for analysis (`gene_list.txt`, `gene_list_reduced.txt`)

### Step 3: Protein Sequence Generation
**Script:** `3_generate_proteins.sh` (SLURM)  
**Runtime:** 2-24 hours  
**Resources:** 4 parallel tasks, 24 CPUs, 32GB RAM

Generates 4 datasets simultaneously:
1. **Reduced pairs** - Canonical + truncated from curated sites
2. **Reduced mutations** - With missense/nonsense/frameshift variants
3. **Full pairs** - Canonical + truncated from all sites  
4. **Full mutations** - With mutation variants

### Step 4: Subcellular Localization Prediction
**Script:** `4_predict_localization.sh` (SLURM)  
**Runtime:** 1-6 hours  
**Resources:** 1 GPU, 8 CPUs, 48GB RAM

Runs DeepLoc predictions in both Fast and Accurate modes for all protein sequences.

## Data Structure

After running the pipeline, your data will be organized as:

```
swissisoform/
├── data/
│   ├── genome_data/
│   │   ├── GRCh38.p7.genome.fa
│   │   ├── gencode.v25.annotation.ensembl_cleaned.gtf
│   │   └── hela_top_transcript.txt
│   └── ribosome_profiling/
│       ├── full_truncations_JL_cleaned.bed
│       ├── selected_truncations_JL_cleaned.bed
│       ├── gene_list.txt (all genes)
│       └── gene_list_reduced.txt (curated genes)
├── results/
│   ├── reduced/                    # Curated truncation sites
│   │   ├── protein_sequences.fasta/csv
│   │   ├── protein_sequences_with_mutations.fasta/csv
│   │   └── *_results.csv          # DeepLoc predictions
│   └── full/                       # All truncation sites
│       ├── protein_sequences.fasta/csv
│       ├── protein_sequences_with_mutations.fasta/csv
│       └── *_results.csv          # DeepLoc predictions
└── scripts/
    ├── 0_download_genome.sh
    ├── 1_cleanup_files.sh
    ├── 3_generate_proteins.sh
    ├── 4_predict_localization.sh
    └── run_pipeline.sh
```

## Output Files

### Protein Sequence Datasets

Each dataset contains:

**Standard pairs** (`protein_sequences.*`):
- Canonical proteins (full-length)
- Truncated proteins (alternative start sites)

**Mutation variants** (`protein_sequences_with_mutations.*`):
- Canonical proteins
- Truncated proteins  
- Mutated variants (missense, nonsense, frameshift)

### Localization Predictions

**DeepLoc outputs** (`*_results.csv`):
- Fast mode: Quick predictions for all sequences
- Accurate mode: High-accuracy predictions for all sequences
- Includes confidence scores and subcellular compartment assignments

### File Formats

**FASTA**: Ready for sequence analysis tools
```
>GENE_TRANSCRIPT_VARIANT description
PROTEIN_SEQUENCE
```

**CSV**: Structured data with metadata
```csv
gene,transcript_id,variant_id,sequence,length,variant_type
NAXE,ENST00000123456,canonical,MAKTG...,245,canonical
NAXE,ENST00000123456,trunc_AUG_100_200,KTGFL...,180,truncated
```

## Advanced Usage

### Single Gene Analysis

For analyzing individual genes or small sets:

```python
import asyncio
from swissisoform.genome import GenomeHandler
from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.translation import TruncatedProteinGenerator
from swissisoform.mutations import MutationHandler
from swissisoform.utils import load_preferred_transcripts

# Initialize components
genome = GenomeHandler(
    'data/genome_data/GRCh38.p7.genome.fa',
    'data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf'
)

alt_isoforms = AlternativeIsoform()
alt_isoforms.load_bed('data/ribosome_profiling/selected_truncations_JL_cleaned.bed')

mutation_handler = MutationHandler()
preferred_transcripts = load_preferred_transcripts(
    'data/genome_data/hela_top_transcript.txt'
)

# Initialize protein generator
protein_generator = TruncatedProteinGenerator(
    genome_handler=genome,
    alt_isoform_handler=alt_isoforms,
    output_dir='output/',
    mutation_handler=mutation_handler
)

# Generate protein sequences for a gene
gene_name = "NAXE"
pairs = protein_generator.extract_gene_proteins(gene_name, preferred_transcripts)

if pairs:
    for pair in pairs:
        print(f"Transcript: {pair['transcript_id']}")
        print(f"Canonical protein length: {len(pair['canonical']['protein'])}")
        print(f"Truncated protein length: {len(pair['truncated']['protein'])}")
        print()
```

### Custom Analysis Pipeline

You can also run components individually:

```bash
# Generate only specific datasets
python3 translate.py gene_list.txt output_dir/ \
  --genome data/genome_data/GRCh38.p7.genome.fa \
  --annotation data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf \
  --bed data/ribosome_profiling/selected_truncations_JL_cleaned.bed \
  --include-canonical \
  --include-mutations \
  --impact-types "missense variant"

# Run visualization for specific genes
python3 analyze_truncations.py gene_list.txt output_dir/ \
  --visualize \
  --preferred-transcripts data/genome_data/hela_top_transcript.txt
```

## License

MIT License - see LICENSE file for details.

## Support

For questions, issues, or contributions:
- Open an issue on GitHub