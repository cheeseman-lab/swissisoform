# CLAUDE.md

Instructions for Claude Code when working in this repository.

## Project

**SwissIsoform** v0.1.0 — Analysis pipeline for alternative protein isoforms from ribosome profiling data. Takes alternative truncation/start sites, generates protein sequences, annotates with mutations (ClinVar, gnomAD, COSMIC), predicts subcellular localization (DeepLoc), and produces prioritized summaries.

## Current Focus

**Branch**: `scripts_refactor` — Active cleanup and refactoring of pipeline scripts.
**TODO**: `TODO.md` — Prioritized task list (last updated 2026-01-16).
**Key issues**: Translation correctness for certain proteins, pipeline speed.

## Development

```bash
eval "$(conda shell.bash hook)" && conda activate swissisoform
uv pip install -e .  # if needed
```

Requires `.env` with API credentials (see `.env.example`): NCBI, gnomAD, COSMIC.

## Pipeline (6 steps, sequential)

```bash
scripts/0_download_genome.sh     # Download GENCODE genome + annotations
scripts/1_cleanup_files.sh       # Process input BED files, prepare isoforms
scripts/2_analyze_mutations.sh   # ClinVar, gnomAD, COSMIC mutation annotation
scripts/3_generate_proteins.sh   # Translate truncated + canonical sequences
scripts/4_predict_localization.sh  # DeepLoc 2.1 (requires GPU, separate conda env)
scripts/5_summarize_results.sh   # Aggregate results, generate reports
```

## Code Structure

```
swissisoform/
├── src/swissisoform/        # Main package
│   ├── alternative_isoforms.py  # BED file processing, isoform logic
│   ├── genome.py            # GENCODE genome/annotation handling
│   ├── mutations.py         # ClinVar/gnomAD/COSMIC integration
│   ├── translation.py       # Protein sequence generation
│   ├── summary.py           # Result aggregation
│   ├── visualize.py         # Plotting
│   ├── config.py            # Configuration management
│   └── utils.py
├── scripts/                 # Pipeline shell scripts (0-5)
├── data/
│   ├── genome_data/         # GENCODE reference files
│   ├── mutation_data/       # ClinVar, gnomAD, COSMIC downloads
│   └── ribosome_profiling/  # Input BED files
├── results/                 # Pipeline outputs
└── reports/                 # Generated reports
```

## Data Conventions

- **Input**: BED files with alternative truncation/start sites from ribosome profiling
- **Intermediate/output**: Parquet for large dataframes, TSV for small readable ones
- Custom parquet files auto-rename old column names to new prefixed format

## Code Style

- Linter: `ruff` (Google-style docstrings)
- Build: setuptools
- Tests: `pytest`

## Known Issues

- Some proteins not translating correctly (active investigation)
- Multi-dataset analysis with mixed GENCODE versions not yet supported (workaround: run separately)
- DeepLoc step requires separate conda environment with GPU access
