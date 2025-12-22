# 5' UTR Variants in Genes with Truncated Isoforms: Biological Analysis

This document analyzes the 30 variants (≥10bp) in genes that also have truncated isoforms detected in HeLa ribosome profiling data. The hypothesis is that 5' UTR length changes may affect translation start site selection, shifting the balance between full-length and N-terminally truncated protein isoforms.

## Tiered Analysis of Variants

### Tier 1: Strong Evidence for Alternative Translation (VHL)

VHL dominates with **10 variants**, all in a gene with well-characterized alternative translation:

**VHL Biology:**
- Produces two major isoforms: **VHL30** (full-length, 213 aa) and **VHL19** (160 aa, from internal Met54)
- Also has a **CUG-initiated VHLα isoform** with distinct tumor suppressor functions
- Clinical evidence: Homozygous mutation abolishing VHL19 (p.Met54Ile) causes erythrocytosis with high EPO but NO tumors - demonstrating isoform balance is functionally critical

| Gene | Variant | Size | Clinical | Mechanism |
|------|---------|------|----------|-----------|
| VHL | c.-39_-10del | -30bp | VUS | Disrupts CUG start region |
| VHL | c.-50_-28del | -23bp | VUS | Disrupts CUG start region |
| VHL | c.-54_-35del | -20bp | VUS | Disrupts CUG start region |
| VHL | c.-62_-43del | -20bp | VUS | Disrupts CUG start region |
| VHL | c.-61_-51dup | +11bp | Conflicting | May alter spacing to CUG |
| VHL | c.-54_-44dup | +11bp | Likely Path | May alter spacing to CUG |
| VHL | c.-11_10dup | +11bp | VUS | Near canonical AUG |
| VHL | c.-12_2dup | +12bp | Conflicting | Spans AUG start |
| VHL | c.-50_-28dup | +23bp | VUS | Duplicates CUG region |
| VHL | c.-54_-35dup | +20bp | VUS | Duplicates CUG region |

**Threshold Crossings:** 4 VHL deletions cross the 40bp and/or 50bp thresholds, directly in the eIF4F blind spot range.

---

### Tier 2: Known Regulatory Elements / Strong Candidates

| Gene | Variant | Size | Clinical | Rationale |
|------|---------|------|----------|-----------|
| CSTB | c.-36_-24del | -13bp | Benign | Near dodecamer repeat region; EPM1 caused by reduced expression. Complex 5' UTR with repeat expansions. |
| GNAS | c.-30_-4del | -27bp | Likely benign | Complex imprinted locus with multiple alternative promoters producing distinct isoforms (Gsα, XLαs, NESP55). Extensive alternative initiation documented. |
| GNAS | c.-23_-14del | -10bp | Likely benign | Same mechanism - GNAS is a paradigm for alternative translation initiation. |
| DMD | c.-164_-132del | -33bp | VUS | Multiple tissue-specific isoforms (Dp427m/b/p, Dp260, Dp140, Dp116, Dp71) from internal promoters. IRES-mediated translation from exon 5 produces functional N-truncated dystrophin. 5' UTR critical for regulation. |

**CSTB Note:** Crosses the 30bp threshold. CSTB mutations cause progressive myoclonic epilepsy (EPM1), and the gene has documented translational regulation through its 5' UTR repeat region.

**GNAS Note:** One of the most complex loci in the human genome for alternative isoform production. Changes to 5' UTR could affect which of multiple overlapping transcripts is translated.

**DMD Note:** The 33bp deletion is substantial. DMD isoform balance is tissue-specific and functionally critical - muscle vs brain vs retinal isoforms have distinct roles.

---

### Tier 3: Plausible Candidates (Less Literature Support)

| Gene | Variant | Size | Clinical | Notes |
|------|---------|------|----------|-------|
| ANKH | c.-128_-80delinsGCC | -46bp | VUS | Pyrophosphate transporter; large deletion near start. Mutations cause craniometaphyseal dysplasia. |
| ANKH | c.-17_-2dup | +16bp | Likely benign | Near start codon |
| NQO1 | c.-36_-9del | -28bp | Likely benign | Antioxidant enzyme with documented translational regulation; polymorphisms affect cancer risk. |
| DSG2 | c.-41_-21del | -21bp | VUS | Desmosomal protein; arrhythmogenic cardiomyopathy gene. |
| DSG2 | c.-11_11dup | +11bp | VUS | Spans start codon - could affect initiation efficiency. |
| DHCR7 | c.-182_-171del | -12bp | VUS | Smith-Lemli-Opitz syndrome gene; cholesterol biosynthesis. |
| PRKAR1A | c.-32_-16dup | +17bp | VUS | PKA regulatory subunit; Carney complex. |
| ACTB | c.-25_-9dup | +17bp | Likely benign | Beta-actin; housekeeping gene with tight translational control. |
| TMCO1 | c.-54_-35dup | +20bp | Likely Path | ER calcium channel; associated with glaucoma and cerebro-facio-thoracic dysplasia. |
| RECQL | c.-421_-407del | -15bp | Benign | RecQ helicase family; genome stability. |
| CHMP1A | c.-94_-81del | -14bp | Likely benign | ESCRT component; chromatin remodeling. |
| SH3PXD2B | c.-142_-131dup | +12bp | Benign | Frank-ter Haar syndrome. |
| LAPTM4B | c.-67_-52dup | +16bp | Benign | Lysosomal protein; overexpressed in cancers. |
| AGTPBP1 | c.-185_-184ins | +15bp | Benign | Neurodegeneration gene (Purkinje cell degeneration). |
| AGTPBP1 | c.11_26dup | +16bp | VUS | Same gene, different variant near start. |
| NUCKS1 | c.-42_-41ins | +26bp | VUS | Nuclear casein kinase substrate; cell cycle regulation. |

---

## Summary Statistics

| Tier | Variants | Genes | Key Findings |
|------|----------|-------|--------------|
| Tier 1 | 10 | 1 (VHL) | CUG-initiated isoform with distinct functions; clinical evidence for isoform balance importance |
| Tier 2 | 4 | 3 (CSTB, GNAS, DMD) | Well-documented alternative translation or tissue-specific isoforms |
| Tier 3 | 16 | 13 | Plausible candidates based on gene function and variant location |

---

## Variants Crossing Biologically Relevant Thresholds

These variants cause the 5' UTR to cross length thresholds associated with altered translation efficiency:

| Gene | Variant | WT UTR | Mut UTR | Thresholds Crossed | Has Truncated Isoform |
|------|---------|--------|---------|-------------------|----------------------|
| VHL | c.-39_-10del | 70bp | 40bp | 40, 50 | Yes |
| VHL | c.-50_-28del | 70bp | 47bp | 50 | Yes |
| VHL | c.-54_-35del | 70bp | 50bp | 50 | Yes |
| VHL | c.-62_-43del | 70bp | 50bp | 50 | Yes |
| CSTB | c.-36_-24del | 39bp | 26bp | 30 | Yes |

All 5 threshold-crossing variants with truncated isoforms are in **Tier 1 or Tier 2** genes with documented alternative translation mechanisms.

---

## Literature Sources

- Von Hippel-Lindau tumor suppressor isoforms: [Wikipedia](https://en.wikipedia.org/wiki/Von_Hippel%E2%80%93Lindau_tumor_suppressor)
- VHL nonstop mutations and isoform balance: [Science Advances](https://www.science.org/doi/10.1126/sciadv.adr6375)
- VHL internal translation initiation: [PubMed 26224408](https://pubmed.ncbi.nlm.nih.gov/26224408/)
- Ribosome scanning and start codon selection: [Nature Communications](https://www.nature.com/articles/s41467-021-26923-3)
- 5' UTR translational control: [PMC7422601](https://pmc.ncbi.nlm.nih.gov/articles/PMC7422601/)
- Leaky scanning mechanism: [PMC390293](https://pmc.ncbi.nlm.nih.gov/articles/PMC390293/)
- Pathogenic 5' UTR variants census: [Frontiers in Molecular Biosciences](https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2023.1257550/full)
- DMD IRES-mediated translation: [Nature Medicine](https://www.nature.com/articles/nm.3628)
- uORF annotation in OMIM genes: [Nucleic Acids Research](https://academic.oup.com/nar/article/51/3/1229/6991039)
