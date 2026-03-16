# SwissIsoform TODO

## High Priority

### Mutation calling improvements & validation runs
**Goal**: Enhance mutation data quality, then validate with full test runs

**Mutation enhancements**:
- [x] Add homozygous counts from gnomAD (`allele_count_hom`, `allele_count_hemi`)
- [x] Filter out 0 allele count variants from gnomAD
- [x] Add additional mutation-level data for ClinVar, gnomAD, COSMIC
  - gnomAD: `allele_count_hemi` (hemizygous for X-linked)
  - ClinVar: `review_status`, `submission_count`
- [ ] Integrate pediatric cancer data (MSK - full pipeline support)

**Validation runs** (after cleanup):
- [ ] Full HeLa test run
- [ ] MSK novel data: reformat as parquet, run as `hela_msk` (similar to `hela_bch` workflow)

---

## Medium Priority

### Multi-dataset analysis with mixed GENCODE versions
**Goal**: Support combined BED files from multiple datasets using different GENCODE GTF versions

**Issue**: Combined BED file contains features from multiple cell lines (hela, rpe1, ipsc, hff), each processed with a different GENCODE GTF version. Current pipeline assumes single GTF per analysis run.

**Tasks**:
- [ ] Investigate GENCODE version differences across datasets (v25, v44, etc.)
- [ ] Determine approach: multi-GTF support vs. unified preprocessing
- [ ] Implement transcript ID mapping if needed
- [ ] Update step 2 to handle combined datasets
- [ ] Test with `combined_isoforms_with_transcripts.bed`

**Workaround**: Run individual datasets separately

---

## Low Priority

### Documentation cleanup
- [ ] Clean up repository (remove old test files, unused scripts)
- [ ] Merge `data/README.md` into main `README.md`
- [ ] Add clear quickstart guide and pipeline workflow docs

### uORF protein sequence support
**Goal**: Generate protein sequences and mutational outputs for uORF isoforms
**Status**: uORFs detected but excluded from analysis pipeline; add support when needed

### Add phyloP conservation scores to protein generation
- [ ] Integrate phyloP scores into step 3 (generate proteins)
- [ ] Score variants/regions by evolutionary conservation

### Improve result summarization (step 5)
**Goal**: Make summary output more interpretable and actionable for prioritizing isoforms
- [ ] Generate plots/visualizations
- [ ] Rank isoforms by biological relevance metrics
- [ ] Surface high-priority candidates more clearly

### Interactive viewer
- [ ] Explore options for interactive visualization of results

---

*Last updated: 2026-01-16*

---

## Recent Changes (2026-01-16)

**Mutation data improvements**:
- Renamed per-mutation columns with source prefixes for clarity:
  - `gnomad_allele_frequency`, `gnomad_allele_count`, `gnomad_allele_count_hom`, `gnomad_allele_count_hemi`
  - `cosmic_sample_count`
  - `clinvar_clinical_significance`, `clinvar_review_status`, `clinvar_submission_count`
- Added `gnomad_allele_count_hemi` (hemizygous count for X-linked variants)
- gnomAD variants with 0 allele count are now filtered out
- Added aggregated `gnomad_allele_count_hom` to isoform-level results
- Custom parquet files auto-rename old column names to new prefixed names

**New ClinVar fields for clinical interpretation**:
- `clinvar_star_rating` (0-4) - confidence level derived from review_status
- `clinvar_condition` - associated disease/condition name
- `clinvar_last_evaluated` - date of last evaluation
