# SwissIsoform Test Suite Guide

Complete guide for testing the SwissIsoform pipeline with custom gene lists.

## Quick Start (3 Steps)

### 1. Find Good Test Genes

```bash
# Find suggested test genes
python tests/find_test_genes.py

# Save suggestions to a file
python tests/find_test_genes.py --output tests/suggested_genes.txt
```

This will show genes categorized by:
- Positive strand (both truncations and extensions)
- Negative strand (both truncations and extensions)
- Truncated only
- Extended only
- Complex genes (multiple transcripts/isoforms)

### 2. Add Genes to Test List

Edit `tests/test_genes.txt` and add your chosen genes:

```bash
# Example - add a few genes
cat > tests/test_genes.txt << 'EOF'
# Test genes for SwissIsoform pipeline

# Positive strand with both types
GENE1
GENE2

# Negative strand with both types
GENE3
GENE4

# Edge cases
GENE5
EOF
```

### 3. Run the Test

```bash
# Basic test
bash tests/run_pipeline_test.sh

# With debug output
bash tests/run_pipeline_test.sh --debug

# With visualizations
bash tests/run_pipeline_test.sh --visualize

# Clean old results and run with all options
bash tests/run_pipeline_test.sh --clean --debug --visualize
```

## What Gets Tested

The test suite runs both major pipeline components:

### 1. Mutation Analysis (`analyze_mutations.py`)
- Identifies transcript-isoform pairs for each gene
- Queries mutation databases (gnomAD, ClinVar, COSMIC)
- Finds mutations in alternative isoform regions
- Generates gene-level and pair-level summaries

**Output:**
- `mutations/gene_level_results.csv` - Summary per gene
- `mutations/truncation_level_results.csv` - Detailed pairs with mutations

### 2. Protein Generation (`generate_proteins.py`)
- Generates protein sequences for all isoforms
- Applies mutations to create mutant proteins
- Validates protein sequences (start with M, no premature stops)
- Exports in CSV and FASTA formats

**Output:**
- `proteins/protein_sequences_with_mutations.csv` - All sequences (CSV)
- `proteins/protein_sequences_with_mutations.fasta` - All sequences (FASTA)

## Test Scenarios to Cover

### Recommended Test Cases

| Scenario | Purpose | Expected Result |
|----------|---------|-----------------|
| **Positive strand + both types** | Test forward strand processing | Both truncated and extended isoforms generated |
| **Negative strand + both types** | Test reverse complement handling | Both truncated and extended isoforms generated |
| **Truncated only** | Test truncation-specific logic | Only truncated isoforms generated |
| **Extended only** | Test extension-specific logic | Only extended isoforms generated |
| **Out-of-frame isoform** | Test validation logic | Should fail validation, no protein generated |
| **Gene with many mutations** | Test mutation integration | Many mutated protein variants |
| **Gene with no mutations** | Test without mutations | Only canonical/alternative pairs, no mutants |

## Files Created

```
tests/
├── README.md                 # Detailed documentation
├── TEST_SUITE_GUIDE.md      # This file (quick reference)
├── test_genes.txt           # Your test gene list (edit this!)
├── find_test_genes.py       # Helper to find good genes
└── run_pipeline_test.sh     # Main test script

results/test_run_TIMESTAMP/
├── test_config.txt          # Test configuration
├── SUMMARY.md               # Test summary
├── mutations/
│   ├── gene_level_results.csv
│   └── truncation_level_results.csv
└── proteins/
    ├── protein_sequences_with_mutations.csv
    └── protein_sequences_with_mutations.fasta
```

## Command Options

### `run_pipeline_test.sh`

```bash
--debug      # Enable detailed logging (shows all processing steps)
--visualize  # Generate visualization plots for mutations
--clean      # Remove previous test results before running
--help       # Show help message
```

### `find_test_genes.py`

```bash
--bed PATH           # Path to BED file (default: data/ribosome_profiling/isoforms_with_transcripts.bed)
--output PATH        # Save suggested genes to file
--max-per-category N # Maximum genes to show per category (default: 10)
```

## Interpreting Results

### Gene-Level Results (`mutations/gene_level_results.csv`)

Key columns:
- `gene_name` - Gene symbol
- `status` - success or error
- `total_transcripts` - Number of transcripts analyzed
- `alternative_features` - Number of alternative start sites found
- `transcript_feature_pairs` - Transcript-isoform pairs generated
- `mutations_*` - Mutation counts by type

**What to check:**
- ✓ All genes have `status: success`
- ✓ `transcript_feature_pairs > 0` for genes with isoforms
- ✓ Mutation counts match expectations

### Transcript-Isoform Pairs (`mutations/truncation_level_results.csv`)

Key columns:
- `transcript_id` - Ensembl transcript ID
- `feature_type` - Truncated or Extended
- `feature_position` - Genomic position of alternative start
- `mutation_count_total` - Total mutations in this region
- `variant_ids_*` - Specific variant IDs by type

**What to check:**
- ✓ Both canonical and alternative pairs present
- ✓ Mutation counts reasonable (0 is normal for some genes)
- ✓ Variant IDs present when mutations found

### Protein Sequences (`proteins/protein_sequences_with_mutations.csv`)

Key columns:
- `gene` - Gene symbol
- `transcript_id` - Ensembl transcript ID
- `variant_id` - Isoform identifier
- `variant_type` - Type of sequence:
  - `canonical` - Standard protein
  - `alternative` - Truncated/extended protein
  - `canonical_mutated` - Canonical + mutation
  - `alternative_mutated` - Alternative + mutation
- `sequence` - Amino acid sequence
- `length` - Sequence length

**What to check:**
- ✓ All sequences start with 'M' (methionine)
- ✓ No premature stop codons (internal '*')
- ✓ Reasonable lengths (>20 AA for valid proteins)
- ✓ Both canonical and alternative variants present

## Validation Checks

### Quick Validation Commands

```bash
# Count sequences by type
tail -n +2 results/test_run_*/proteins/protein_sequences_with_mutations.csv | \
  cut -d',' -f6 | sort | uniq -c

# Check all sequences start with M
tail -n +2 results/test_run_*/proteins/protein_sequences_with_mutations.csv | \
  cut -d',' -f8 | grep -v '^M' | wc -l
# Should output: 0

# Check for premature stops
tail -n +2 results/test_run_*/proteins/protein_sequences_with_mutations.csv | \
  cut -d',' -f8 | grep '\*[A-Z]' | wc -l
# Should output: 0

# Count genes processed
tail -n +2 results/test_run_*/mutations/gene_level_results.csv | \
  cut -d',' -f1 | wc -l

# Count successful genes
tail -n +2 results/test_run_*/mutations/gene_level_results.csv | \
  grep ',success,' | wc -l
```

## Troubleshooting

### No results for a gene?

1. **Check gene is in BED file:**
   ```bash
   grep "^GENENAME" data/ribosome_profiling/isoforms_with_transcripts.bed
   ```

2. **Check gene name spelling:**
   - Use exact gene symbol from BED file
   - Case-sensitive!

3. **Check error messages:**
   ```bash
   grep "GENENAME" results/test_run_*/mutations/gene_level_results.csv
   ```

### Missing alternative isoforms?

This is often expected! Alternative isoforms may fail validation if:
- Out of frame (creates frameshifts)
- Generates premature stop codons
- Results in very short proteins (<20 AA)

Enable debug mode to see why:
```bash
bash tests/run_pipeline_test.sh --debug 2>&1 | grep "GENENAME" | grep -i "fail"
```

### No mutations found?

This is normal! Many genes don't have mutations in the databases, or mutations don't overlap with alternative isoform regions.

Check mutation databases are configured:
```bash
ls -lh data/mutation_data/
```

## Example Workflow

```bash
# 1. Find good test genes
python tests/find_test_genes.py --output tests/suggested_genes.txt

# 2. Pick a few genes and add to test list
head -20 tests/suggested_genes.txt | grep -v '^#' > tests/test_genes.txt

# 3. Run test with debug to see details
bash tests/run_pipeline_test.sh --debug --clean

# 4. Check results
cd results/test_run_*/
less mutations/gene_level_results.csv
less proteins/protein_sequences_with_mutations.csv

# 5. Validate protein sequences
grep -c '^M' proteins/protein_sequences_with_mutations.csv
```

## Best Practices

1. **Start small** - Test with 3-5 genes first
2. **Check results** - Validate before adding more genes
3. **Use debug mode** - Understand what's happening
4. **Keep test genes** - Save successful gene lists for regression testing
5. **Document edge cases** - Note genes that fail and why

## Next Steps

After testing:

1. **Review results** - Check all expected outputs generated
2. **Validate sequences** - Ensure protein quality
3. **Check edge cases** - Note any failures and investigate
4. **Expand tests** - Add more genes to cover more scenarios
5. **Automate** - Set up regular test runs

## Getting Help

- **Test documentation:** `tests/README.md`
- **Logging guide:** `docs/LOGGING.md`
- **Pipeline docs:** Main project README
- **Debug logs:** Run with `--debug` flag

---

**Remember:** The test suite is designed to be flexible. Start with a few genes, validate the results, then expand your testing as needed!
