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

Basic example of visualizing an isoform:

```python
from swissisoform.genome import GenomeHandler
from swissisoform.visualize import GenomeVisualizer
from swissisoform.isoform import AlternativeIsoform

# Initialize handlers
genome = GenomeHandler(
    'data/genome_data/hg38.fa',
    'data/genome_data/hg38.ncbiRefSeq.gtf'
)

# Load alternative isoform data
alt_isoforms = AlternativeIsoform()
alt_isoforms.load_bed('path/to/isoforms.bed')

# Get features for visualization
alt_features = alt_isoforms.get_visualization_features('GENE_NAME')

# Create visualization
visualizer = GenomeVisualizer(genome)
visualizer.visualize_transcript(
    'GENE_NAME', 
    'TRANSCRIPT_ID',
    alt_features=alt_features,
    output_file='output.png'
)
```

For more detailed examples and use cases, see the Jupyter notebook at `notebooks/visualize_isoforms_NAXE.ipynb`.

## License

MIT License - see LICENSE file for details.

