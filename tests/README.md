# SwissIsoform Test Suite

This directory contains the test suite for the SwissIsoform pipeline.

- ABHD18 has a premature stop (positive extension)
- DMD has weird mutations that do not functionally work (negative strand truncation)
- PNPO is a positive strand truncation and extension gene
- SIX5 is a negative strand truncation and extension gene
- PSRC1 is a truncation with single bases at cutoffs between exons, represents an edge case

## Quick Start

1. **Add test genes** to `test_genes.txt`:
   ```bash
   echo "TP53" >> tests/test_genes.txt
   echo "PTEN" >> tests/test_genes.txt
   ```

2. **Run the test**:
   ```bash
   bash tests/run_pipeline_test.sh
   ```

3. **Review results** in `results/test_run_TIMESTAMP/`

## Running Tests

### Basic Test
```bash
bash tests/run_pipeline_test.sh
```

### With Debug Logging
```bash
bash tests/run_pipeline_test.sh --debug
```

### With Visualizations
```bash
bash tests/run_pipeline_test.sh --visualize
```

### Clean Previous Results
```bash
bash tests/run_pipeline_test.sh --clean
```

### Combine Options
```bash
bash tests/run_pipeline_test.sh --debug --visualize --clean
```

## Output Structure

After running, results are saved to `results/test_run_TIMESTAMP/`:

```
results/test_run_YYYYMMDD_HHMMSS/
├── test_config.txt           # Test configuration
├── SUMMARY.md                # Test summary
├── mutations/
│   ├── gene_level_results.csv         # Gene-level mutation summary
│   ├── truncation_level_results.csv   # Transcript-isoform pairs with mutations
│   └── visualizations/                # (if --visualize enabled)
└── proteins/
    ├── protein_sequences_with_mutations.csv    # Protein sequences (CSV format)
    └── protein_sequences_with_mutations.fasta  # Protein sequences (FASTA format)
```

## Interpreting Results

### Mutation Analysis (`mutations/gene_level_results.csv`)

Check for:
- Total mutations found per gene
- Mutations in alternative isoform regions
- Breakdown by mutation type (missense, nonsense, frameshift, etc.)

### Protein Generation (`proteins/protein_sequences_with_mutations.csv`)

Check for:
- **variant_type** column:
  - `canonical` - Standard protein from annotated start
  - `alternative` - Protein from alternative start (truncation/extension)
  - `canonical_mutated` - Canonical protein with mutation applied
  - `alternative_mutated` - Alternative protein with mutation applied

- **Sequence validation:**
  - All sequences should start with 'M' (methionine)
  - No premature stop codons (except at the end)
  - Reasonable length (>20 amino acids for valid proteins)

### Common Issues to Check

1. **No results for a gene**
   - Gene may not be in the BED file
   - Gene name spelling may differ
   - Check gene_level_results.csv for error messages

2. **Missing alternative isoforms**
   - May be out-of-frame (failing validation)
   - Check debug logs with `--debug` flag
   - Look for "validation failed" messages

3. **No mutations found**
   - Gene may not have mutations in the mutation database
   - Mutations may not overlap with alternative isoform regions
   - Normal for some genes

## Validating Test Results

### Step 1: Check Gene-Level Results

```bash
# View gene-level summary
column -t -s, results/test_run_*/mutations/gene_level_results.csv | less -S
```

Expected columns:
- `gene_name`: Gene symbol
- `status`: success/error
- `total_transcripts`: Number of transcripts analyzed
- `alternative_features`: Number of alternative start sites
- `transcript_feature_pairs`: Number of transcript-isoform pairs
- `mutations_*`: Mutation counts by type

### Step 2: Check Transcript-Isoform Pairs

```bash
# View detailed pairs
column -t -s, results/test_run_*/mutations/truncation_level_results.csv | less -S
```

Expected columns include:
- `transcript_id`, `feature_type`, `feature_position`
- `mutation_count_total`: Total mutations in this isoform region
- `mutations_*`: Breakdown by mutation type
- `variant_ids_*`: Specific variant identifiers

### Step 3: Check Protein Sequences

```bash
# Count sequences by type
tail -n +2 results/test_run_*/proteins/protein_sequences_with_mutations.csv | \
  cut -d',' -f6 | sort | uniq -c
```

Expected variant types:
- `canonical`: Standard proteins
- `alternative`: Truncated/extended proteins
- `canonical_mutated`: Mutated canonical
- `alternative_mutated`: Mutated alternative

### Step 4: Validate Protein Quality

```bash
# Check all proteins start with M
tail -n +2 results/test_run_*/proteins/protein_sequences_with_mutations.csv | \
  cut -d',' -f8 | grep -v '^M' || echo "All sequences start with M"

# Check for premature stops (should be empty)
tail -n +2 results/test_run_*/proteins/protein_sequences_with_mutations.csv | \
  cut -d',' -f8 | grep '\*' | grep -v '\*$' || echo "No premature stop codons"
```

## Debugging Failed Tests

### Enable Debug Logging

```bash
bash tests/run_pipeline_test.sh --debug 2>&1 | tee test_debug.log
```

### Check for Common Errors

1. **File Not Found**
   - Verify paths in test_config.txt
   - Check that genome data is downloaded

2. **Gene Not Found**
   - Verify gene name spelling
   - Check gene is in BED file: `grep "GENENAME" data/ribosome_profiling/isoforms_with_transcripts.bed`

3. **Validation Failures**
   - Look for "validation failed" in debug logs
   - Check if alternative isoforms are out-of-frame
   - This is expected for some edge cases

4. **No Mutations**
   - Normal for genes without mutations in database
   - Check mutation database is configured correctly