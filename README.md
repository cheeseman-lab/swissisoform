# SwissIsoform

SwissIsoform is a Python package for analyzing alternative protein isoforms discovered through ribosome profiling data. It integrates genomic annotations, truncation sites, and mutation data to generate comprehensive protein sequence datasets for downstream analysis.

## Features

- **Genome sequence and annotation handling** with GENCODE support
- **Alternative truncation site processing** from ribosome profiling BED files
- **Mutation analysis** in truncation regions from ClinVar
- **Protein sequence generation** from canonical and truncated transcripts
- **Mutation integration** into protein sequences
- **Interactive visualizations** of transcript features and mutations
- **Subcellular localization prediction** using DeepLoc

## Installation

```bash
git clone https://github.com/cheeseman-lab/swissisoform.git
cd swissisoform

# Create conda environment
conda env create --file=environment.yml
conda activate swissisoform

# For localization prediction (optional)
# 1. Download DeepLoc 2.1 from: https://services.healthtech.dtu.dk/services/DeepLoc-2.0/
# 2. Place DeepLoc-2.1.0.tar.gz in the main swissisoform/ directory
# 3. Install DeepLoc environment
conda create -n deeploc python=3.8
conda activate deeploc
pip install DeepLoc-2.1.0.tar.gz
```

## Required Input Files

SwissIsoform requires only **two input files** to run the analysis:

1. **Ribosome profiling BED file**: `data/ribosome_profiling/truncations.bed`
   - Format: Standard BED with truncation sites
   - Example file provided: Contains ribosome profiling truncation sites
   
2. **Preferred transcripts file**: `data/genome_data/hela_top_transcript.txt`
   - Format: One transcript ID per line (e.g., `ENST00000000000.1`)
   - Example file provided: Contains HeLa cell preferred transcripts

**Note**: Example files for both requirements are included in the repository and can be used as-is or replaced with your own experimental data.

## Quick Start

### Full Pipeline (SLURM)

```bash
# Download reference data
bash scripts/0_download_genome.sh

# Clean ribosome profiling data  
bash scripts/1_cleanup_files.sh

# Analyze mutations (SLURM)
sbatch scripts/2_analyze_mutations.sh

# Generate protein sequences (SLURM)
sbatch scripts/3_generate_proteins.sh

# Predict localization (SLURM, GPU)
sbatch scripts/4_predict_localization.sh
```

### Interactive Analysis (Jupyter)

#### Mutation Analysis

```python
import asyncio
from swissisoform.genome import GenomeHandler
from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.mutations import MutationHandler
from swissisoform.utils import load_preferred_transcripts

# Initialize components
genome = GenomeHandler(
    'data/genome_data/GRCh38.p7.genome.fa',
    'data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf'
)

alt_isoforms = AlternativeIsoform()
alt_isoforms.load_bed('data/ribosome_profiling/truncations_cleaned.bed')

mutation_handler = MutationHandler()
preferred_transcripts = load_preferred_transcripts('data/genome_data/hela_top_transcript.txt')

# Analyze mutations for a gene
async def analyze_gene_mutations(gene_name):
    result = await mutation_handler.analyze_gene_mutations_comprehensive(
        gene_name=gene_name,
        genome_handler=genome,
        alt_isoform_handler=alt_isoforms,
        output_dir='output/',
        visualize=True,
        impact_types={"clinvar": ["missense variant"]},
        preferred_transcripts=preferred_transcripts
    )
    return result

# Run analysis
result = await analyze_gene_mutations("NAXE")
print(f"Found {result['transcript_truncation_pairs']} transcript-truncation pairs")
print(f"Mutations in truncation regions: {result['mutations_filtered']}")
```

#### Protein Generation

```python
from swissisoform.translation import TruncatedProteinGenerator

# Initialize protein generator
protein_generator = TruncatedProteinGenerator(
    genome_handler=genome,
    alt_isoform_handler=alt_isoforms,
    output_dir='output/',
    mutation_handler=mutation_handler
)

# Generate protein sequences with mutations
async def generate_proteins_with_mutations(gene_list):
    dataset = await protein_generator.create_protein_sequence_dataset_with_mutations(
        gene_list=gene_list,
        preferred_transcripts=preferred_transcripts,
        include_mutations=True,
        impact_types=["missense variant"],
        output_format="fasta,csv",
        min_length=50,
        max_length=2000
    )
    return dataset

# Generate dataset
gene_list = ["NAXE", "ERGIC3", "GARS1"]
dataset = await generate_proteins_with_mutations(gene_list)
print(f"Generated {len(dataset)} protein sequences")

# For pairs-only (no mutations)
pairs_dataset = protein_generator.create_protein_sequence_dataset_pairs(
    gene_list=gene_list,
    preferred_transcripts=preferred_transcripts,
    output_format="fasta,csv",
    min_length=50,
    max_length=2000
)
```

## Data Structure

```
swissisoform/
├── data/
│   ├── genome_data/
│   │   ├── GRCh38.p7.genome.fa
│   │   ├── gencode.v25.annotation.ensembl_cleaned.gtf
│   │   └── hela_top_transcript.txt
│   └── ribosome_profiling/
│       ├── truncations_cleaned.bed
│       ├── gene_list.txt
│       └── gene_list_reduced.txt
└── results/
    ├── reduced/                     # Curated gene set (22 genes)
    │   ├── mutations/
    │   │   ├── gene_level_results.csv
    │   │   ├── truncation_level_results.csv
    │   │   └── [gene_name]/         # Visualizations
    │   └── proteins/
    │       ├── protein_sequences_pairs.fasta/csv
    │       ├── protein_sequences_with_mutations.fasta/csv
    │       └── *_results.csv        # DeepLoc predictions
    └── full/                        # All genes (~1,100 genes)
        ├── mutations/
        └── proteins/
```

## Output Files

### Mutation Analysis
- **`gene_level_results.csv`**: Gene-level mutation summary
- **`truncation_level_results.csv`**: Detailed transcript-truncation pair analysis with mutation counts by impact type
- **Visualizations**: PDF files showing transcript structure, truncation sites, and mutations

### Protein Sequences
- **`protein_sequences_pairs.fasta/csv`**: Canonical + truncated protein pairs
- **`protein_sequences_with_mutations.fasta/csv`**: Sequences with integrated mutations

### Localization Predictions
- **`*_results.csv`**: DeepLoc predictions (Fast and Accurate modes)

## File Formats

**FASTA sequences:**
```
>GENE_TRANSCRIPT_VARIANT description
PROTEIN_SEQUENCE
```

**CSV data:**
```csv
gene,transcript_id,variant_id,sequence,length,variant_type
NAXE,ENST00000123456,canonical,MAKTG...,245,canonical
NAXE,ENST00000123456,trunc_AUG_100_200,KTGFL...,180,truncated
```

## Requirements

- Linux/Unix system with SLURM (for full pipeline)
- Conda package manager
- GPU access (for fast localization prediction)
- Internet connection (for ClinVar API access)

## License

MIT License - see LICENSE file for details.