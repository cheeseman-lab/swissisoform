# Validation Output Format Guide

## Overview

Validation output is now optimized for **database lookup** and **quality control**. Each variant is logged with information you can use to search in ClinVar, gnomAD, or COSMIC.

## Output Format

### Disagreements (INFO Level)
```
{source}: {variant_id} | {chr}:{pos}:{ref}>{alt} | {database_impact} → {validated_impact} [tags]
```

**Example:**
```
INFO - ClinVar: VCV000012345 | chr17:43044295:C>T | missense variant → synonymous variant
INFO - gnomAD: 17-43044300-G-A | chr17:43044300:G>A | synonymous variant → missense variant [alt start]
INFO - COSMIC: COSV12345678 | chr13:32900001:A>G | missense variant → nonsense variant
```

### Agreements (DEBUG Level Only)
```
DEBUG - {source}: {variant_id} | {chr}:{pos}:{ref}>{alt} | Validated: {impact}
```

**Example:**
```
DEBUG - ClinVar: VCV000067890 | chr17:43045802:G>A | Validated: missense variant
```

### Alt Start Site Mutations (INFO Level)
```
{source}: {variant_id} | {chr}:{pos}:{ref}>{alt} | Alt start site at {codon} codon
```

**Example:**
```
INFO - gnomAD: 17-43044295-C-T | chr17:43044295:C>T | Alt start site at ATG codon
```

### Validation Failures (WARNING Level)
```
{source}: {variant_id} | {chr}:{pos}:{ref}>{alt} | Validation failed: {error}
```

**Example:**
```
WARNING - COSMIC: COSV99999999 | chr17:43044000:T>C | Validation failed: Position not in transcript
```

## Fields Explained

| Field | Description | Purpose |
|-------|-------------|---------|
| `source` | ClinVar, gnomAD, COSMIC, or Custom | Know which database to search |
| `variant_id` | Database-specific identifier | Direct lookup in the source database |
| `chr:pos:ref>alt` | Genomic coordinates | Universal variant identifier (VCF format) |
| `database_impact` | Original impact from source | What the database says |
| `validated_impact` | Our validation result | What sequence analysis shows |
| `[tags]` | Special annotations | `[alt start]` for alternative start sites |

## How to Use This Information

### 1. Look Up a Variant in ClinVar
```bash
# From log: "ClinVar: VCV000012345 | chr17:43044295:C>T | ..."

# Search by accession:
https://www.ncbi.nlm.nih.gov/clinvar/VCV000012345

# Or by position:
https://www.ncbi.nlm.nih.gov/clinvar/?term=17[chr]+AND+43044295[chrpos37]
```

### 2. Look Up a Variant in gnomAD
```bash
# From log: "gnomAD: 17-43044300-G-A | chr17:43044300:G>A | ..."

# Search by position (gnomAD browser):
https://gnomad.broadinstitute.org/variant/17-43044300-G-A
```

### 3. Look Up a Variant in COSMIC
```bash
# From log: "COSMIC: COSV12345678 | chr13:32900001:A>G | ..."

# Search by COSMIC ID:
https://cancer.sanger.ac.uk/cosmic/search?q=COSV12345678

# Or by position in COSMIC browser
```

### 4. Extract All Disagreements
```bash
# From your test log:
grep "→" tests/test_run_*/pipeline_test.log

# Example output:
# INFO - ClinVar: VCV000012345 | chr17:43044295:C>T | missense variant → synonymous variant
# INFO - gnomAD: 17-43044300-G-A | chr17:43044300:G>A | synonymous variant → missense variant
```

### 5. Extract Variants by Source
```bash
# All ClinVar disagreements:
grep "ClinVar:.*→" tests/test_run_*/pipeline_test.log

# All gnomAD disagreements:
grep "gnomAD:.*→" tests/test_run_*/pipeline_test.log

# All COSMIC disagreements:
grep "COSMIC:.*→" tests/test_run_*/pipeline_test.log
```

### 6. Extract Variants by Gene
```bash
# From header: "Validating 50 mutations for BRCA1 (ENST00000357654)"
grep -A 100 "Validating.*BRCA1" tests/test_run_*/pipeline_test.log | grep "→"
```

## Complete Example

```
INFO - Validating 150 mutations for BRCA1 (ENST00000357654)
INFO -   Region: upstream_extension at 43044295-43045802
INFO - ClinVar: VCV000012345 | chr17:43044295:C>T | missense variant → synonymous variant
DEBUG -   Database predicted protein change, validation shows silent
INFO - gnomAD: 17-43044300-G-A | chr17:43044300:G>A | Alt start site at ATG codon
INFO - gnomAD: 17-43044350-T-C | chr17:43044350:T>C | synonymous variant → missense variant [alt start]
DEBUG -   Database predicted silent, validation shows protein change
DEBUG - ClinVar: VCV000012350 | chr17:43044400:A>G | Validated: missense variant
DEBUG - ClinVar: VCV000012351 | chr17:43044450:G>C | Validated: missense variant
... [143 agreements not shown at INFO level]
INFO - Validation complete: 150 total, 147 agree (98.0%), 3 disagree, 0 failed, 2 in alt start sites
INFO - Disagreement breakdown:
INFO -   missense variant → synonymous variant: 1
INFO -   synonymous variant → missense variant: 1
```

## Quality Control Workflow

### Quick Review
```bash
# 1. Check overall agreement rate
grep "Validation complete" tests/test_run_*/pipeline_test.log

# 2. Count disagreements
grep -c "→" tests/test_run_*/pipeline_test.log

# 3. See disagreement patterns
grep "Disagreement breakdown:" -A 10 tests/test_run_*/pipeline_test.log
```

### Deep Dive
```bash
# 1. Extract all disagreements to a file
grep "→" tests/test_run_*/pipeline_test.log > disagreements.txt

# 2. Review each variant
# Copy variant_id or chr:pos:ref>alt to search in database

# 3. Determine if disagreement is:
#    - Database annotation error
#    - Validation error
#    - Biological complexity (e.g., alternative splicing)
```

## Benefits of New Format

1. **Database Searchable**: Each line contains IDs you can paste directly into databases
2. **Parseable**: Pipe-delimited format easy to grep/awk/parse
3. **Concise**: Agreements are silent (DEBUG only), focus on issues
4. **Informative**: Contains gene, source, ID, position, and impact in one line
5. **Professional**: No emojis, no all-caps, clean format

## Comparison to Old Format

### Old (Verbose, Hard to Search)
```
DEBUG - [1/150] Validating VCV000012345
DEBUG - Position: 43044295
DEBUG - Mutation: C>T
DEBUG - Original impact: 'missense variant'
DEBUG - Validated impact: 'synonymous variant'
DEBUG - 🔄 DISAGREEMENT: missense variant → synonymous variant
DEBUG - 📝 Explanation: Database predicted protein change, but codon analysis shows silent mutation
```

**Issues**:
- 7 lines per variant
- Hard to grep
- No chromosome info
- Can't copy/paste to search

### New (Concise, Searchable)
```
INFO - ClinVar: VCV000012345 | chr17:43044295:C>T | missense variant → synonymous variant
```

**Benefits**:
- 1 line per variant
- Easy to grep
- All info for database lookup
- Copy chr:pos:ref>alt or ID to search
