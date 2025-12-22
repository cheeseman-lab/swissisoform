# 5' UTR Length Change Analysis

Analysis of ClinVar variants that affect 5' UTR length in MANE Select transcripts, identifying variants where indels cause the UTR to cross biologically relevant length thresholds.

## Data Sources

- **ClinVar** (downloaded December 15, 2024)
  - Total variants: 8,375,541
  - GRCh38 variants: 4,155,543
  - URL: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

- **GENCODE v47**
  - Used for MANE Select transcript annotations and 5' UTR genomic intervals
  - URL: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz

- **RefSeq GRCh38**
  - Used for Ensembl-to-RefSeq transcript ID mapping
  - URL: https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz

- **HeLa Ribosome Profiling Data**
  - Used to identify genes with truncated isoforms (N-terminally truncated proteins from downstream start sites)
  - 1,276 genes with truncated isoforms
  - Source: `data/hela_truncated_genes.txt` (derived from `hela_isoforms_with_transcripts.bed`)

## Approach

### Genomic Coordinate Validation

This analysis uses **genomic coordinate validation** rather than relying solely on HGVS notation.

Previous approaches relied on HGVS notation (e.g., `c.-50del`) to identify 5' UTR variants. This is problematic because:
- HGVS positions may not match GENCODE annotations (different transcript versions can have different UTR boundaries)
- No validation that genomic coordinates fall within the annotated UTR
- Size calculations from HGVS can be inaccurate for insertions and duplications

Our approach:
- Extract 5' UTR genomic intervals from GENCODE v47 MANE Select transcripts
- Validate each ClinVar variant by checking if BOTH start AND stop genomic coordinates fall within the annotated 5' UTR regions
- Calculate size changes using a hybrid approach:
  - **Deletions**: Use genomic coordinates (accurate)
  - **Insertions**: Parse HGVS for inserted sequence length
  - **Duplications**: Parse HGVS for duplicated range
  - **Indels (delins)**: Parse HGVS for net size change (inserted - deleted)

### Threshold Crossing Logic

- **Deletions/shrinking indels**: WT UTR > threshold AND Mutant UTR <= threshold (UTR shrinks below threshold)
- **Insertions/duplications/growing indels**: WT UTR < threshold AND Mutant UTR >= threshold (UTR grows above threshold)

## Results Summary

### Filtering Funnel

**Step 0: Total GRCh38 indels in ClinVar**
- 254,668 variants (deletions, insertions, duplications, indels)

**Step 1: START position in 5' UTR** (1,579 variants)
- Variant start coordinate falls within a GENCODE 5' UTR region
- Deletions: 872
- Insertions: 86
- Duplications: 504
- Indels: 117

**Step 2: BOTH start AND stop in 5' UTR** (988 variants)
- Both coordinates must be within GENCODE 5' UTR regions
- This step excludes 591 variants where start is in UTR but stop is not:
  - 84 deletions spanning from UTR into CDS (e.g., `c.-16_10del`)
  - ~320 large structural variants (whole-gene deletions/duplications, often 10-300kb)
  - These are correctly excluded because the 5' UTR length change is irrelevant when the entire gene is affected
- Deletions: 465
- Insertions: 82
- Duplications: 338
- Indels: 103

**Step 3: Size change >= 10bp** (166 variants)
- Deletions: 79
- Insertions: 8
- Duplications: 56
- Indels: 23

**Step 4: Threshold crossings** (17 unique variants cross at least one threshold)
- 10bp threshold: 2 variants (1 deletion, 1 indel)
- 20bp threshold: 3 variants (1 deletion, 1 insertion, 1 indel)
- 30bp threshold: 8 variants (5 deletions, 1 insertion, 2 duplications)
- 40bp threshold: 7 variants (2 deletions, 1 insertion, 4 duplications)
- 50bp threshold: 7 variants (6 deletions, 1 insertion)

### Threshold-Crossing Variants

**17 variants cross at least one biologically relevant threshold:**

By clinical significance:
- Uncertain significance: 8
- Likely benign: 4
- Benign: 2
- Conflicting classifications: 2
- Pathogenic: 1

By variant type:
- Deletions: 11
- Duplications: 4
- Indels: 1
- Insertions: 1

Notable variants:
- **BCKDHA c.12_13ins...** (73bp insertion): 10bp → 83bp, crosses 20/30/40/50bp thresholds, **Pathogenic**
- **CRYGD c.-44_-23del**: 51bp → 29bp, crosses 30/40/50bp thresholds, Benign
- **ATR c.-29_-9delinsT**: 29bp → 9bp, crosses 10/20bp thresholds, Conflicting
- **VHL c.-39_-10del**: 70bp → 40bp, crosses 40/50bp thresholds, Uncertain significance
- **SGSH c.-39_-16dup24**: 20bp → 44bp, crosses 30/40bp thresholds, Uncertain significance

### Truncated Isoform Overlap

Variants are annotated with whether the gene has a truncated isoform in HeLa ribosome profiling data. A truncated isoform indicates that the gene uses alternative downstream start sites, producing N-terminally truncated proteins.

- 30 of 166 variants (>= 10bp) are in genes with truncated isoforms
- **5 threshold-crossing variants** are in genes with truncated isoforms:
  - VHL: 4 deletions crossing 40/50bp thresholds
  - CSTB: 1 deletion crossing 30bp threshold

This is biologically relevant because 5' UTR length changes may affect which start site is used, potentially shifting the balance between full-length and truncated protein isoforms.

## Output Files

Intermediate files (filtering funnel):
- `5utr_variants_start_in_utr.csv` - Variants where START is in 5' UTR (1,579 variants)
- `5utr_variants_both_in_utr.csv` - Variants where BOTH start AND stop in UTR (988 variants)
- `5utr_variants_full.csv` - Full dataset with UTR calculations (988 variants)
- `5utr_variants_min10bp.csv` - Size >= 10bp with threshold columns (166 variants)
- `5utr_variants_min10bp_cross{N}bp.csv` - Variants crossing specific threshold
- `5utr_length_changes_min10bp_cross{N}bp.png` - Slope plot for each threshold

Output columns:
- `gene_name` - Gene symbol
- `hgvsc` - HGVS cDNA notation
- `variant_type` - deletion, insertion, duplication, or indel
- `start_in_utr` / `stop_in_utr` / `both_in_utr` - UTR containment flags
- `wt_utr_length` - Wild-type 5' UTR length (bp)
- `mutant_utr_length` - Mutant 5' UTR length (bp)
- `size_change` - Net change in bp (negative for deletions)
- `crosses_Xbp` - Boolean flags for each threshold
- `thresholds_crossed` - Comma-separated list of crossed thresholds
- `has_truncated_isoform` - Gene has a truncated isoform in HeLa ribosome profiling data
- `genomic_coords` - Genomic coordinates (chrN:start-stop)
- `clinical_significance` - ClinVar classification
- `clinvar_url` - Link to ClinVar record

## Usage

```bash
# Run with default parameters (min-size=10, thresholds=10,20,30,40,50)
python run_analysis.py

# Custom parameters
python run_analysis.py --min-size 5 --thresholds "20,40,60,80"

# Force re-download of data files
python run_analysis.py --force-download
```

Arguments:
- `--min-size` (default: 10) - Minimum variant size in bp
- `--thresholds` (default: "10,20,30,40,50") - Comma-separated threshold values
- `--force-download` - Force re-download of data files

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
