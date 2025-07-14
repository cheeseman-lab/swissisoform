# SwissIsoform

SwissIsoform is a Python package for analyzing alternative protein isoforms discovered through ribosome profiling data. It integrates genomic annotations, truncation sites, and mutation data to generate comprehensive protein sequence datasets for downstream analysis and machine learning applications.

## Features

- **Genome sequence and annotation handling** with GENCODE support
- **Alternative truncation site processing** from ribosome profiling BED files
- **Comprehensive mutation analysis** in truncation regions from ClinVar and gnomAD
- **Protein sequence generation** from canonical and truncated transcripts
- **Mutation integration** into protein sequences with impact filtering
- **Interactive visualizations** of transcript features and mutations
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

# Step 2: Analyze mutations in truncation regions (SLURM job)
sbatch 2_analyze_mutations.sh

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

### Step 2: Mutation Analysis
**Script:** `2_analyze_mutations.sh` (SLURM)  
**Runtime:** 1-8 hours  
**Resources:** 2 parallel tasks, 16 CPUs, 24GB RAM

Analyzes mutations in alternative isoform truncation regions:
1. **Reduced dataset** - Curated truncation sites with selected genes
2. **Full dataset** - All truncation sites with complete gene list

**Features:**
- Fetches mutation data from ClinVar
- Filters for missense, nonsense, and frameshift variants
- Creates transcript-truncation pair analysis
- Generates visualizations for each gene
- Produces detailed mutation statistics

**Outputs:**
- Gene-level mutation summary (`gene_level_results.csv`)
- Transcript-truncation pair details (`truncation_level_results.csv`)
- Interactive visualizations (PDF files organized by gene/transcript)

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
│   ├── mutations/                  # Step 2 outputs
│   │   ├── reduced/
│   │   │   ├── gene_level_results.csv
│   │   │   ├── truncation_level_results.csv
│   │   │   └── [gene_name]/       # Visualizations by gene
│   │   └── full/
│   │       ├── gene_level_results.csv
│   │       ├── truncation_level_results.csv
│   │       └── [gene_name]/       # Visualizations by gene
│   ├── reduced/                    # Step 3 & 4 outputs (curated)
│   │   ├── protein_sequences.fasta/csv
│   │   ├── protein_sequences_with_mutations.fasta/csv
│   │   └── *_results.csv          # DeepLoc predictions
│   └── full/                       # Step 3 & 4 outputs (all sites)
│       ├── protein_sequences.fasta/csv
│       ├── protein_sequences_with_mutations.fasta/csv
│       └── *_results.csv          # DeepLoc predictions
└── scripts/
    ├── 0_download_genome.sh
    ├── 1_cleanup_files.sh
    ├── 2_analyze_mutations.sh
    ├── 3_generate_proteins.sh
    ├── 4_predict_localization.sh
    ├── analyze_mutations.py
    ├── analyze_truncations.py
    ├── translate.py
    └── run_pipeline.sh
```

## Output Files

### Mutation Analysis Results

**Gene-level summary** (`gene_level_results.csv`):
- Overall statistics per gene
- Total transcripts, truncation features, and transcript-truncation pairs
- Summary mutation counts

**Transcript-truncation details** (`truncation_level_results.csv`):
- Detailed analysis for each transcript-truncation pair
- Mutation counts by impact type (missense, nonsense, frameshift)
- ClinVar variant IDs for each mutation category
- Genomic coordinates of truncation regions

**Visualizations** (PDF files):
- Transcript structure with CDS, UTRs, start/stop codons
- Alternative truncation sites marked as brackets
- Mutation positions colored by impact type
- Separate files for full view and zoomed view of truncation regions

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

### Individual Step Analysis

You can run individual pipeline steps for custom analysis:

#### Mutation Analysis Only
```bash
# Analyze mutations for specific genes
python3 analyze_mutations.py gene_list.txt output_dir/ \
  --genome data/genome_data/GRCh38.p7.genome.fa \
  --annotation data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf \
  --bed data/ribosome_profiling/selected_truncations_JL_cleaned.bed \
  --visualize \
  --sources "clinvar" \
  --impact-types "missense variant" "nonsense variant" "frameshift variant"
```

#### Transcript-Mutation Visualization
```bash
# Generate visualizations for transcript features and mutations
python3 analyze_truncations.py gene_list.txt output_dir/ \
  --visualize \
  --preferred-transcripts data/genome_data/hela_top_transcript.txt
```

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

# Comprehensive mutation analysis for a gene
async def analyze_gene_mutations(gene_name):
    result = await mutation_handler.analyze_gene_mutations_comprehensive(
        gene_name=gene_name,
        genome_handler=genome,
        alt_isoform_handler=alt_isoforms,
        output_dir='output/',
        visualize=True,
        impact_types={"clinvar": ["missense variant", "nonsense variant"]},
        preferred_transcripts=preferred_transcripts
    )
    return result

# Run analysis
result = asyncio.run(analyze_gene_mutations("NAXE"))
print(f"Found {result['transcript_truncation_pairs']} transcript-truncation pairs")
print(f"Total mutations in truncation regions: {result['mutations_filtered']}")
```

### Custom Protein Generation

```python
# Initialize protein generator with mutation support
protein_generator = TruncatedProteinGenerator(
    genome_handler=genome,
    alt_isoform_handler=alt_isoforms,
    output_dir='output/',
    mutation_handler=mutation_handler
)

# Generate protein sequences with mutations
async def generate_proteins_with_mutations(gene_name):
    enhanced_pairs = await protein_generator.extract_gene_proteins_with_mutations(
        gene_name=gene_name,
        preferred_transcripts=preferred_transcripts,
        include_mutations=True,
        impact_types=["missense variant", "nonsense variant"]
    )
    
    if enhanced_pairs:
        for pair in enhanced_pairs:
            print(f"Transcript: {pair['transcript_id']}")
            print(f"Canonical protein length: {len(pair['canonical']['protein'])}")
            print(f"Truncated protein length: {len(pair['truncated_base']['protein'])}")
            print(f"Mutation variants: {len(pair['truncated_mutations'])}")
            print()

# Run generation
asyncio.run(generate_proteins_with_mutations("NAXE"))
```

### Custom Dataset Generation

```bash
# Generate only specific datasets
python3 translate.py gene_list.txt output_dir/ \
  --genome data/genome_data/GRCh38.p7.genome.fa \
  --annotation data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf \
  --bed data/ribosome_profiling/selected_truncations_JL_cleaned.bed \
  --include-canonical \
  --include-mutations \
  --impact-types "missense variant" "frameshift variant" \
  --min-length 50 \
  --max-length 2000
```

## Troubleshooting

### Common Issues

**Missing mutation data:**
- Check internet connection for ClinVar API access
- Verify gene names match HGNC symbols
- Consider rate limiting (pipeline includes automatic delays)

**Visualization errors:**
- Ensure matplotlib and required dependencies are installed
- Check that output directories have write permissions
- Verify transcript IDs exist in annotation file

**Memory issues:**
- Reduce number of parallel tasks in SLURM scripts
- Process smaller gene lists
- Increase memory allocation in SLURM headers

### Performance Tips

- Use preferred transcript lists to focus analysis
- Filter by specific mutation impact types to reduce processing time
- Run pipeline steps incrementally for large datasets
- Monitor SLURM job logs for progress and errors

## License

MIT License - see LICENSE file for details.

## Support

For questions, issues, or contributions:
- Open an issue on GitHub
- Check SLURM job logs in `out/` directory for detailed error messa