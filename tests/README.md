# SwissIsoform Test Suite

Comprehensive test suite for the SwissIsoform pipeline.

## Quick Start

1. **Edit test gene list**: Add genes to `test_genes.txt` (one per line)
   ```bash
   echo "ABHD18" >> tests/test_genes.txt
   echo "ADAR" >> tests/test_genes.txt
   ```

2. **Run the test**:
   ```bash
   bash tests/run_pipeline_test.sh
   ```

3. **Review results** in `tests/test_run_TIMESTAMP/`

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
Key columns: `gene`, `transcript_id`, `variant_type`, `sequence`, `length`

**Variant types:**
- `canonical` - Standard protein from annotated start
- `alternative` - Protein from alternative start (truncation/extension)
- `canonical_mutated` - Canonical + mutation
- `alternative_mutated` - Alternative + mutation

**Validate:**
- All sequences start with 'M' (methionine)
- No premature stop codons (internal '*')
- Reasonable lengths (>20 amino acids)

## Validation Commands

```bash
# Count sequences by variant type
tail -n +2 tests/test_run_*/proteins/protein_sequences_with_mutations.csv | \
  cut -d',' -f6 | sort | uniq -c

# Check all sequences start with M
tail -n +2 tests/test_run_*/proteins/protein_sequences_with_mutations.csv | \
  cut -d',' -f8 | grep -v '^M' | wc -l
# Should output: 0

# Check for premature stops
tail -n +2 tests/test_run_*/proteins/protein_sequences_with_mutations.csv | \
  cut -d',' -f8 | grep '\*[A-Z]' | wc -l
# Should output: 0
```

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

The test suite runs on genes listed in `test_genes.txt`:
- **ABHD18** - Positive strand with premature stop (extension)
- **ACSF3** - Positive strand
- **ADAR** - Positive strand
- **CLDND1** - Positive strand
- **CIB1** - Positive strand
- **FGF2** - Positive strand
- **MAPK14** - Negative strand

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
