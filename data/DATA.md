# SwissIsoform Data Requirements

This document describes the required data formats for the SwissIsoform pipeline.

## Directory Structure

```
data/
├── genome_data/              # Reference genome and annotations (GRCh38)
├── mutation_data/            # Mutation databases (ClinVar, gnomAD, COSMIC, or custom)
└── ribosome_profiling/       # Alternative start site (TIS) data from ribosome profiling
```

---

## Alternative Start Site Data (BED Format)

### Overview

Alternative translation initiation sites (TIS) from ribosome profiling experiments should be provided as BED files in `data/ribosome_profiling/`. The pipeline supports two formats.

**Important**: The transcript IDs in your BED file must match the GENCODE GTF annotation version used for pipeline processing. Mismatched versions will cause errors during sequence extraction and validation.

### Transcript and GTF Requirements

1. **Transcript ID matching**: Each feature's transcript ID (e.g., `ENST00000370214.8`) must exist in your GENCODE GTF annotation file
2. **GENCODE version consistency**: Use the same GENCODE version for:
   - Generating TIS data (ribosome profiling analysis)
   - Processing data with step 1 (`1_cleanup_files.sh`)
   - Running the mutation analysis pipeline

3. **Feature type alignment**: The transcript and feature type must be biologically consistent:
   - **Annotated**: Uses the annotated start codon from the GTF
   - **Extended**: Uses an upstream start codon (N-terminal extension)
   - **Truncated**: Uses a downstream start codon (N-terminal truncation)

   The transcript ID you specify should correspond to the isoform that contains the alternative start site.

**Example GTF versions**:
- GENCODE v25 (used for HeLa data processing)
- GENCODE v44 (used for RPE1 data processing)

### Generating TIS Data

TIS data can be generated from ribosome profiling experiments using [RiboTISH](https://github.com/zhpn1024/ribotish), which identifies translation initiation sites and outputs them in BED format compatible with this pipeline.

### Format 1: BED5 (Recommended)

Minimal format - sequences will be extracted from the genome during step 1 processing.

**Tab-separated, no header, 5 columns:**

| Column | Description | Example |
|--------|-------------|---------|
| `chr` | Chromosome | `chr1` |
| `start` | Start position (0-based, inclusive) | `94418479` |
| `end` | End position (0-based, exclusive) | `94418481` |
| `name` | Feature identifier (see format below) | `ATG_Annotated_ABCD3_ENST00000370214.8` |
| `strand` | Strand | `+` or `-` |

**Feature name format:**
```
{StartCodon}_{Type}_{GeneName}_{TranscriptID}
```

Where:
- **StartCodon**: `ATG`, `GTG`, `CTG`, `TTG`, `ACG`, `AGG`, etc.
- **Type**: `Annotated`, `Extended`, `Truncated`, or `uORF`
- **GeneName**: Gene symbol (e.g., `ABCD3`)
- **TranscriptID**: Ensembl transcript ID with version from your GTF (e.g., `ENST00000370214.8`)

**Example:**
```
chr1	229558650	229558652	ATG_Annotated_ABCB10_ENST00000344517.4	-
chr1	94418479	94418481	ATG_Annotated_ABCD3_ENST00000370214.8	+
chr1	94418395	94418397	GTG_Extended_ABCD3_ENST00000370214.8	+
```

### Format 2: BED7 (With Pre-computed Sequences)

Tab-separated, no header, 7 columns: chr, start, end, name, strand, dna_sequence, protein_sequence

Use this format if you've already extracted and translated sequences. The DNA sequence must be the full coding sequence from start to stop codon (on positive strand, reverse complemented for negative strand genes). This is not used in the pipeline, but can be used for validation that translation is working accurately.

### File Naming

Name files as: `{DATASET}_{Study}.bed`

Examples:
- `HELA_Ly2024.bed`
- `IPSC_Chen2020.bed`
- `RPE1_Ly2024.bed`

The dataset name (before underscore) will be used in output file naming.

### Processing

After placing BED files in `data/ribosome_profiling/`, run:

```bash
bash scripts/1_cleanup_files.sh
```

This processes raw BED files and:
- Extracts sequences from the genome (for BED5 format)
- Adds transcript and RefSeq annotations from GTF
- Validates transcript IDs exist in GTF
- Generates analysis-ready files (e.g., `hela_isoforms_with_transcripts.bed`)

### Requirements

- **Coordinate system**: BED format (0-based start, 0-based exclusive end)
- **Genome assembly**: GRCh38/hg38
- **GTF compatibility**: Transcript IDs must match your GENCODE GTF version
- **Unique features**: Each feature name must be unique
- **Valid coordinates**: start < end, both non-negative
- **Chromosome format**: `chr1`-`chr22`, `chrX`, `chrY`, `chrM`

---

## Custom Mutation Data (Parquet Format)

### Overview

Custom mutation databases can be provided as Parquet files in `data/mutation_data/`. These are used alongside or instead of ClinVar, gnomAD, and COSMIC.

### Usage

```bash
# Via SLURM
cd scripts/

sbatch --export=DATASET=your_dataset,OUTPUT_NAME="analysis_your_dataset",SOURCES="custom",CUSTOM_PARQUET="data/mutation_data/your_variants.parquet" \
  2_analyze_mutations.sh
```

### Required Columns

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `chromosome` | str | Chromosome (with or without 'chr' prefix) | `chr11` or `11` |
| `position` | int | Genomic position (1-based, GRCh38) | `154627850` |
| `reference` | str | Reference allele | `G` |
| `alternate` | str | Alternate allele | `A` |
| `variant_id` | str | Unique variant identifier | `ADAR_11_154627850_G_A` |
| `gene_name` | str | Gene symbol (must match your gene list) | `ADAR` |
| `source` | str | **Must be "custom"** | `custom` |

### Optional Columns

| Column | Type | Description |
|--------|------|-------------|
| `occurrence_count` | int | Number of times variant was observed (will be summed) |
| `source_dataset` | str | Original dataset name for tracking (e.g., `Dataset_A`) |

**Important**:
- All variants must have `source = "custom"` for proper aggregation
- Use `source_dataset` to preserve original source information if combining multiple datasets
- File format: Apache Parquet (`.parquet`)
- Coordinates: GRCh38/hg38 (1-based)

### Output Integration

Custom variants will appear in mutation results with these columns:
- `count_custom` - Total number of custom variants in the feature region
- `ids_custom_missense_variant` - List of missense variant IDs
- `ids_custom_nonsense_variant` - List of nonsense variant IDs
- (Similar `ids_custom_*` columns for other variant types)

If `occurrence_count` is provided, the pipeline will sum it and report in results.

---

## Reference Genome

The pipeline uses **GRCh38/hg38** (Genome Reference Consortium Human Build 38).

All custom mutation and TIS data must use GRCh38 coordinates. To convert from GRCh37/hg19, use [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver).

---

Last updated: 2025-11-10
