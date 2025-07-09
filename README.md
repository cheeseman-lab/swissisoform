# SwissIsoform

SwissIsoform is a Python package for analyzing and visualizing alternative protein isoforms discovered through ribosome profiling data. It provides tools for integrating genomic, mutation, and isoform data to generate comprehensive visualizations of protein variants.

## Features

- Genome sequence and annotation handling
- Alternative isoform visualization
- Integration of mutation data from multiple sources:
  - gnomAD
  - ClinVar
  - Custom aggregator files
- Configurable visualization options for different feature types
- Support for zoomed visualization of specific regions

## Installation

1. Clone the repository:
```bash
git clone https://github.com/cheeseman-lab/swissisoform.git
cd swissisoform
```

2. Create and activate the conda environment:
```bash
conda env create --file=environment.yml
conda activate swissisoform
```

## Data Requirements

The package requires several data files to function:

### 1. Reference Genome and Annotations

Download the required genome files:
```bash
cd data/genome_data
bash download_genome.sh
```

This will download:
- Human reference genome (hg38.fa)
- NCBI RefSeq annotations (hg38.ncbiRefSeq.gtf)

### 2. Alternative Start Sites

- Required format: BED file
- Expected fields: chromosome, start, end, name (format: ENSGXXXXX_GENE_STARTCODON_TYPE)

### 3. Mutation Data (Optional)

The package can integrate mutation data from:
- gnomAD (accessed via API)
- ClinVar (accessed via API)
- Custom aggregator files (CSV format)

## Usage

### Basic Example

Here's a comprehensive example demonstrating the new approach to analyzing and visualizing isoforms with mutation data:

```python
from swissisoform.genome import GenomeHandler
from swissisoform.visualize import GenomeVisualizer
from swissisoform.isoform import AlternativeIsoform
from swissisoform.mutations import MutationHandler
from swissisoform.utils import analyze_mutations
import os

# Initialize handlers
genome = GenomeHandler(
    'data/genome_data/hg38.fa',
    'data/genome_data/hg38.ncbiRefSeq.gtf'
)

# Load alternative isoform data
alt_isoforms = AlternativeIsoform()
alt_isoforms.load_bed('path/to/isoforms.bed')

# Initialize mutation handler
mutation_handler = MutationHandler()

# Process a specific gene
gene_name = "NAXE"
output_dir = f"./{gene_name}/"
os.makedirs(output_dir, exist_ok=True)

# Get alternative isoform features
alt_features = alt_isoforms.get_visualization_features(gene_name)

# Get transcript information
transcript_info = genome.get_transcript_ids(gene_name)

# Analyze mutations with different filtering options
mutations_unfiltered = await analyze_mutations(
    gene_name=gene_name, 
    mutation_handler=mutation_handler,
    alt_features=alt_features, 
    sources=["clinvar"],
)

# Filter mutations by impact type
impact_types = {
    "clinvar": ["missense variant", "nonsense variant", "frameshift variant"],
}
mutations_filtered = await analyze_mutations(
    gene_name=gene_name, 
    mutation_handler=mutation_handler,
    alt_features=alt_features, 
    sources=["clinvar"], 
    impact_types=impact_types,
)

# Create visualizations
visualizer = GenomeVisualizer(genome)

# Generate visualizations for each transcript
for _, transcript in transcript_info.iterrows():
    transcript_id = transcript['transcript_id']
    
    # Standard view with unfiltered mutations
    visualizer.visualize_transcript(
        gene_name=gene_name,
        transcript_id=transcript_id,
        alt_features=alt_features,
        mutations_df=mutations_unfiltered,
        output_file=f'{output_dir}{transcript_id}_unfiltered.png'
    )
    
    # Standard view with filtered mutations
    visualizer.visualize_transcript(
        gene_name=gene_name,
        transcript_id=transcript_id,
        alt_features=alt_features,
        mutations_df=mutations_filtered,
        output_file=f'{output_dir}{transcript_id}_filtered.png'
    )
    
    # Zoomed view for detailed analysis
    visualizer.visualize_transcript_zoomed(
        gene_name=gene_name,
        transcript_id=transcript_id,
        alt_features=alt_features,
        mutations_df=mutations_filtered,
        output_file=f'{output_dir}{transcript_id}_filtered_zoom.png',
        padding=100
    )
```

For more detailed examples and use cases, see the Jupyter notebooks in `notebooks/`.

## Deeploc installation

```bash
# Create and activate a conda environment (recommended)
conda create -n deeploc python=3.8
conda activate deeploc

# Then install DeepLoc 2.1
pip install DeepLoc-2.1.0.tar.gz

# Run results_reduced
deeploc2 -f scripts/results_reduced/protein_sequence_dataset_pairs.fasta -o scripts/results_reduced/

# Run results
deeploc2 -f scripts/results/protein_sequence_dataset_pairs.fasta -o scripts/results/
```

## License

MIT License - see LICENSE file for details.