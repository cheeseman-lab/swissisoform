# Debug Logging Guide for SwissIsoform Pipeline Tests

## Overview

The test pipeline now includes comprehensive debug logging to validate each mutation and track its effects through the pipeline.

## Running Tests with Different Verbosity Levels

### Minimal Output (Default)
```bash
bash tests/run_pipeline_test.sh
```
Shows only warnings and errors.

### Standard Output (Recommended)
```bash
bash tests/run_pipeline_test.sh --clean
```
Shows INFO level: disagreements, alt start sites, summaries.

### Full Debug Output
```bash
bash tests/run_pipeline_test.sh --debug --clean
```
Shows DEBUG level: includes agreements, detailed explanations.

### What the `--debug` flag does:

1. **Enables DEBUG-level logging** (`-vv`) in both:
   - `analyze_mutations.py`
   - `generate_proteins.py`

2. **Captures all output** to `pipeline_test.log` in the test run folder

3. **Shows detailed validation** for each mutation including:
   - Position mapping (genomic → coding)
   - Reference/alternate alleles
   - Original impact from database
   - Validated impact from sequence analysis
   - Agreement/disagreement tracking
   - Alternative start site detection

## Log Output Structure

All output is saved to: `tests/test_run_YYYYMMDD_HHMMSS/pipeline_test.log`

### Key Debug Sections:

#### 1. Mutation Fetching
```
INFO - Found 150 gnomAD variants for DMD
INFO - Found 89 ClinVar variants for DMD
INFO - Found 234 COSMIC variants for DMD
```

#### 2. Individual Mutation Validation (Disagreements Only)
```
INFO - Validating 473 mutations for ENST00000357033
INFO -   Disagreement: chr1:12345678:C>T - missense variant → synonymous variant
DEBUG -     Note: Database predicted protein change, codon analysis shows silent
INFO -   Disagreement: chr1:12346000:A>T - synonymous variant → missense variant [alt start]
DEBUG -     Note: Database predicted silent, codon analysis shows protein change
```

Note: Agreements are NOT logged to reduce verbosity. Only disagreements, alt start site mutations, and failures are shown.

#### 3. Alternative Start Site Detection
```
INFO -   Alt start site mutation: chr1:12345601:A>G at ATG codon (12345600)
```

#### 4. Validation Summary
```
INFO - Validation complete: 473 total, 450 agree (95.1%), 21 disagree, 2 failed, 15 in alt start sites
INFO - Disagreement breakdown:
INFO -   missense variant → synonymous variant: 12
INFO -   synonymous variant → missense variant: 6
INFO -   missense variant → nonsense variant: 3
```

## Mutation Sources

The test pipeline now queries all three variant sources:
- **ClinVar** - Clinical variants
- **gnomAD** - Population variants
- **COSMIC** - Cancer variants

Source tracking is included in the output CSV files with columns:
- `mutation_sources` - comma-separated list of sources
- `mutations_clinvar` - count from ClinVar
- `mutations_gnomad` - count from gnomAD
- `mutations_cosmic` - count from COSMIC

## Logging Levels

The pipeline uses three verbosity levels:

| Flag | Level | What's Logged |
|------|-------|--------------|
| None | WARNING | Errors and warnings only |
| `-v` | INFO | Summary statistics and progress |
| `-vv` | DEBUG | Full mutation-by-mutation validation |

The `--debug` flag in the test script automatically uses `-vv`.

## GQL Library Logging

The GraphQL library (used for gnomAD queries) has special handling:
- **Non-debug mode**: Suppressed to WARNING (no query/response spam)
- **Debug mode**: Set to INFO (shows request summary, not full JSON)

This prevents the log from being flooded with GraphQL internals while still allowing you to see what's happening.

## Logging Levels and Verbosity

### INFO Level (Default with -v)
Shows:
- Summary statistics (variant counts, validation results)
- **Disagreements only** (database vs validation mismatches)
- Alternative start site mutations
- Validation failures (warnings)
- One-line summary per transcript

### DEBUG Level (With -vv)
Adds:
- Detailed explanations for disagreements
- Cache building status
- Position mapping details
- Codon analysis (only for disagreements)

**Key Design Principle**: Agreements are silent. The focus is on quality control by highlighting **discrepancies** that need review.

## Troubleshooting

### No debug logs appearing

Check that you're using `--debug` flag:
```bash
bash tests/run_pipeline_test.sh --debug --visualize --clean
```

### Log file too large

Debug mode generates extensive output. For production runs without validation details, omit the `--debug` flag:
```bash
bash tests/run_pipeline_test.sh --visualize --clean
```

### Want to see only mutation validation

Use `grep` to filter the log:
```bash
grep "Validating\|AGREEMENT\|DISAGREEMENT" tests/test_run_*/pipeline_test.log
```

## Output Files

Each test run creates:
```
tests/test_run_YYYYMMDD_HHMMSS/
├── pipeline_test.log          # Complete output log with all debug info
├── test_config.txt            # Test configuration
├── SUMMARY.md                 # Test summary
├── mutations/
│   ├── gene_level_results.csv
│   └── isoform_level_results.csv
└── proteins/
    ├── protein_sequences_with_mutations.csv
    └── protein_sequences_with_mutations.fasta
```

## Next Steps

After running a test with debug logging:

1. **Review the log** for validation issues:
   ```bash
   less tests/test_run_*/pipeline_test.log
   ```

2. **Check for disagreements**:
   ```bash
   grep "DISAGREEMENT" tests/test_run_*/pipeline_test.log
   ```

3. **Verify alternative start sites**:
   ```bash
   grep "alternative start site" tests/test_run_*/pipeline_test.log
   ```

4. **Examine validation failures**:
   ```bash
   grep "VALIDATION FAILED" tests/test_run_*/pipeline_test.log
   ```
