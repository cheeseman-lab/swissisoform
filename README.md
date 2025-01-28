# SwissIsoform

Tool for analyzing and visualizing alternative protein isoforms from ribosome profiling data.

## Installation

```bash
# Clone repository
git clone https://github.com/cheeseman-lab/swissisoform.git
cd swissisoform

# Create conda environment
conda env create --file=environment.yml

# Activate environment
conda activate swissisoform
```

## Data Requirements

1. Reference genome (FASTA)
2. Gene annotations (GTF)

```bash
# Move into genome_data folder
cd data/genome_data

# Download genome data
sh download_genome.sh
```

3. Alternative start sites (BED)



## Usage

See `notebooks/visualize_isoforms.ipynb` for example usage.