# SwissIsoform Test Suite

Comprehensive test suite for the SwissIsoform pipeline.

## Prerequisites

Before running the test suite, ensure you have:

1. **Processed HeLa BED file**: Run steps 0 and 1 to generate `data/ribosome_profiling/hela_isoforms_with_transcripts.bed`
   ```bash
   # Step 0: Download genome data
   bash scripts/0_setup_genome_data.sh

   # Step 1: Process ribosome profiling data
   bash scripts/1_cleanup_files.sh
   ```

2. **Test gene list**: Edit `tests/test_genes.txt` with genes of interest (one per line)

## Quick Start

```bash
# Run the full pipeline test
bash tests/run_pipeline_test.sh

# Or with verbose logging
bash tests/run_pipeline_test.sh -vv
```

**Review results** in `tests/test_run_TIMESTAMP/`

## Command Options

```bash
bash tests/run_pipeline_test.sh [-v|-vv] [--visualize] [--clean]

Options:
  -v           Info logging (shows progress and summaries)
  -vv          Debug logging (shows detailed validation)
  --visualize  Generate mutation visualization plots
  --clean      Remove previous test results before running
```

## What Gets Tested

The pipeline runs two stages:

### 1. Mutation Analysis
- Identifies alternative start sites (truncations/extensions)
- Queries mutation databases (ClinVar, gnomAD, COSMIC)
- Outputs: `mutations/gene_level_results.csv` and `mutations/isoform_level_results.csv`

### 2. Protein Generation
- Generates canonical and alternative protein sequences
- Applies mutations to create mutant variants
- Outputs: `proteins/protein_sequences_with_mutations.csv` and `.fasta`

## Output Structure

```
tests/test_run_YYYYMMDD_HHMMSS/
├── test_config.txt              # Test configuration
├── pipeline_test.log            # Complete log output
├── SUMMARY.md                   # Test summary
├── mutations/
│   ├── gene_level_results.csv           # Summary per gene
│   └── isoform_level_results.csv        # Transcript-isoform pairs
└── proteins/
    ├── protein_sequences_with_mutations.csv    # All sequences (CSV)
    └── protein_sequences_with_mutations.fasta  # All sequences (FASTA)
```

## Interpreting Results

### Gene-Level Results
Key columns: `gene_name`, `status`, `total_transcripts`, `alternative_features`, `mutations_*`

**Check for:**
- All genes have `status: success`
- `alternative_features > 0` for genes with isoforms
- Mutation counts by type (missense, nonsense, etc.)

### Protein Sequences
Key columns: `gene_name`, `transcript_id`, `feature_id`, `bed_name`, `feature_type`, `variant_type`, `sequence`, `length`

**Variant types:**
- `canonical` - Standard protein from annotated start
- `extension` - Protein from upstream alternative start (N-terminal extension)
- `truncation` - Protein from downstream alternative start (N-terminal truncation)
- `canonical_mutated` - Canonical + missense mutation
- `extension_mutated` - Extension + missense mutation
- `truncation_mutated` - Truncation + missense mutation

## Logging Levels

### Standard Output (default)
Shows warnings and errors only.

### Info Level (`-v`)
Adds:
- Progress indicators and summaries
- Validation disagreements (database vs. sequence analysis)
- Alternative start site mutations

### Debug Level (`-vv`)
Adds:
- Detailed explanations for disagreements
- Position mapping details
- Codon-level analysis

**Note**: Agreements between database and validation are silent (not logged) to reduce noise.

## Troubleshooting

### No results for a gene?
1. Check gene is in BED file:
   ```bash
   grep "GENENAME" data/ribosome_profiling/hela_isoforms_with_transcripts.bed
   ```
2. Verify exact gene symbol (case-sensitive)
3. Check error messages in `gene_level_results.csv`

### Missing alternative isoforms?
This is often expected! Alternatives may fail validation if:
- Out of frame (creates frameshifts)
- Contains premature stop codons
- Results in very short proteins (<20 AA)

Enable debug mode to see why:
```bash
bash tests/run_pipeline_test.sh -vv 2>&1 | grep "GENENAME" | grep -i "fail"
```

### No mutations found?
Normal! Many genes don't have mutations in databases, or mutations don't overlap alternative isoform regions.

## Current Test Genes

The test suite runs on genes listed in `tests/test_genes.txt`:
- **ADAR** (chr1, negative strand) - Multiple extensions and truncations, tests cache handling for complex features
- **DYNLRB1** (chr20, positive strand) - Multiple truncations
- **RBM10** (chrX, positive strand) - Multiple extensions, tests overlapping extension regions
- **TCP1** (chr6, negative strand) - 4 extensions with overlapping coordinates, critical for testing cache key fix

## Test Configuration

- **Dataset**: HeLa ribosome profiling data
- **Genome**: GRCh38.p7
- **Annotation**: GENCODE v25 (with v47 gene name updates if available)
- **BED file**: `data/ribosome_profiling/hela_isoforms_with_transcripts.bed`
- **Mutation sources**: ClinVar, gnomAD, COSMIC

## Best Practices

1. **Start small** - Test with 3-5 genes first
2. **Validate results** - Check outputs before expanding
3. **Use appropriate logging** - `-v` for quality control, `-vv` for debugging
4. **Document edge cases** - Note genes that fail and investigate why
5. **Keep successful lists** - Save working gene lists for regression testing
