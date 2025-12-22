# 5' UTR Length Change Analysis

Analysis of ClinVar variants that affect 5' UTR length in MANE Select transcripts, identifying variants where indels/duplications cause the UTR to cross biologically relevant length thresholds.

## Approach

### Genomic Coordinate Validation

This analysis uses **genomic coordinate validation** rather than relying solely on HGVS notation. This provides more accurate and defensible results.

#### Why Genomic Validation?

Previous approaches relied on HGVS notation (e.g., `c.-50del`) to identify 5' UTR variants. This is problematic because:

1. **HGVS positions may not match GENCODE annotations** - Different transcript versions can have different UTR boundaries
2. **No validation that genomic coordinates fall within UTR** - A variant annotated as `c.-50` might not actually be in the annotated 5' UTR
3. **Size calculations from HGVS can be inaccurate** - Especially for insertions and duplications where ClinVar's genomic coordinates represent insertion points, not sequence lengths

#### Our Approach

1. **Extract 5' UTR genomic intervals from GENCODE** (v47, MANE Select transcripts)
2. **Validate each ClinVar variant** by checking if BOTH start AND stop genomic coordinates fall within the annotated 5' UTR regions
3. **Calculate size changes using a hybrid approach**:
   - **Deletions**: Use genomic coordinates (accurate)
   - **Insertions**: Parse HGVS for inserted sequence length
   - **Duplications**: Parse HGVS for duplicated range

### Filtering Criteria

1. **UTR containment**: Both variant start AND stop positions must be within GENCODE-annotated 5' UTR regions
2. **Minimum size**: Configurable (default: 10bp)
3. **Threshold crossing**: Check if variant causes UTR to cross specified thresholds (default: 10, 20, 30, 40, 50bp)

### Threshold Crossing Logic

- **Deletions**: WT UTR > threshold AND Mutant UTR <= threshold (UTR shrinks below threshold)
- **Insertions/Duplications**: WT UTR < threshold AND Mutant UTR >= threshold (UTR grows above threshold)

## Data Sources

| Source | Description | URL |
|--------|-------------|-----|
| ClinVar | Variant annotations | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz |
| GENCODE | v47 GTF annotations | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz |
| RefSeq | GRCh38 GTF (for Ensembl-RefSeq mapping) | https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz |

## Results Summary

### Filtering Funnel

| Step | Filter | Deletions | Insertions | Duplications | Total |
|------|--------|-----------|------------|--------------|-------|
| 0 | Total GRCh38 indels | - | - | - | 236,935 |
| 1 | **START in 5' UTR** | 872 | 86 | 504 | **1,462** |
| 2 | **BOTH start AND stop in UTR** | 465 | 82 | 338 | **885** |
| 3 | Size >= 10bp | 79 | 8 | 56 | 143 |

The filtering funnel shows how many variants are retained at each step:
- **Step 1**: Variants where the START position falls within a GENCODE 5' UTR region (more permissive)
- **Step 2**: Variants where BOTH start AND stop positions are within 5' UTR (stricter, excludes UTR→CDS spanning variants)
- **Step 3**: Size filter applied

### Threshold Crossings (size >= 10bp)

| Threshold | Variants | Deletions | Insertions | Duplications |
|-----------|----------|-----------|------------|--------------|
| 10bp | 1 | 1 | 0 | 0 |
| 20bp | 2 | 1 | 1 | 0 |
| 30bp | 8 | 5 | 1 | 2 |
| 40bp | 7 | 2 | 1 | 4 |
| 50bp | 7 | 6 | 1 | 0 |

### Notable Variants

16 variants cross at least one threshold:

| Gene | HGVS | Type | WT UTR | Mut UTR | Thresholds Crossed | Clinical Significance |
|------|------|------|--------|---------|-------------------|----------------------|
| BCKDHA | c.12_13ins... (73bp) | insertion | 10 | 83 | 20,30,40,50 | Pathogenic |
| CRYGD | c.-44_-23del | deletion | 51 | 29 | 30,40,50 | Benign |
| VHL | c.-39_-10del | deletion | 70 | 40 | 40,50 | Uncertain significance |
| GALNT12 | c.2_13dup | duplication | 29 | 41 | 30,40 | Uncertain significance |
| SGSH | c.-39_-16dup24 | duplication | 20 | 44 | 30,40 | Uncertain significance |

## Output Files

### Intermediate Files (Filtering Funnel)

| File | Description | Count |
|------|-------------|-------|
| `5utr_variants_start_in_utr.csv` | Variants where START is in 5' UTR | 1,462 |
| `5utr_variants_both_in_utr.csv` | Variants where BOTH start AND stop in UTR | 885 |
| `5utr_variants_full.csv` | Full dataset with UTR calculations | 885 |
| `5utr_variants_min10bp.csv` | Size >= 10bp (includes threshold columns) | 143 |
| `5utr_variants_min10bp_cross{N}bp.csv` | Variants crossing specific threshold | varies |
| `5utr_length_changes_min10bp_cross{N}bp.png` | Slope plot for each threshold | - |

### Output Columns

| Column | Description |
|--------|-------------|
| `gene_name` | Gene symbol |
| `hgvsc` | HGVS cDNA notation |
| `variant_type` | deletion, insertion, or duplication |
| `start_in_utr` | Boolean: is START position in 5' UTR? |
| `stop_in_utr` | Boolean: is STOP position in 5' UTR? |
| `both_in_utr` | Boolean: are BOTH start AND stop in UTR? |
| `wt_utr_length` | Wild-type 5' UTR length (bp) |
| `mutant_utr_length` | Mutant 5' UTR length (bp) |
| `size_change` | Net change in bp (negative for deletions) |
| `crosses_Xbp` | Boolean: does variant cross X bp threshold? |
| `thresholds_crossed` | Comma-separated list of crossed thresholds |
| `genomic_coords` | Genomic coordinates (chrN:start-stop) |
| `clinical_significance` | ClinVar classification |
| `clinvar_url` | Link to ClinVar record |

## Usage

```bash
# Run with default parameters (min-size=10, thresholds=10,20,30,40,50)
python run_analysis.py

# Custom parameters
python run_analysis.py --min-size 5 --thresholds "20,40,60,80"

# Force re-download of data files
python run_analysis.py --force-download
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--min-size` | 10 | Minimum variant size in bp |
| `--thresholds` | "10,20,30,40,50" | Comma-separated threshold values |
| `--force-download` | False | Force re-download of data files |

## Methods Details

### 5' UTR Identification

For each MANE Select transcript in GENCODE v47:
1. Parse UTR features and start_codon positions
2. Determine 5' vs 3' UTR based on strand:
   - **Plus strand (+)**: 5' UTR regions are BEFORE start codon (end < start_codon_pos)
   - **Minus strand (-)**: 5' UTR regions are AFTER start codon (start > start_codon_pos)
3. Store genomic intervals: (chromosome, start, end, strand)

### Variant Validation

For each ClinVar indel:
1. Get genomic coordinates (Chromosome, Start, Stop)
2. Look up gene's 5' UTR intervals
3. Check if BOTH Start AND Stop positions fall within any UTR interval
4. Only include variants that pass this check

### Size Calculation

| Variant Type | Method | Rationale |
|--------------|--------|-----------|
| Deletion | Genomic coordinates | ClinVar Start/Stop accurately represents deleted region |
| Insertion | Parse HGVS sequence | ClinVar coords show insertion point, not inserted length |
| Duplication | Parse HGVS range | ClinVar coords show insertion point, not duplicated length |

## Biological Context

The 5' UTR plays a critical role in translation regulation. Short 5' UTRs (< 40-50bp) may:
- Reduce translation efficiency
- Affect ribosome scanning
- Impact upstream open reading frame (uORF) usage

Variants that cause the 5' UTR to cross these length thresholds may have functional consequences even if they don't directly affect the coding sequence.

## Requirements

- Python 3.8+
- pandas
- matplotlib
- wget (for downloading data files)
