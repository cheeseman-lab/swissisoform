# SwissIsoform v2 — Orchestration Design for Parallel Rewrite

**Date:** 2026-04-14
**Author:** Matteo Di Bernardo + Claude Code
**Status:** Design approved, pending implementation plan
**Repo:** `swissisoform-v2` (new, clean build at `/lab/barcheese01/mdiberna/swissisoform-v2/`)

---

## 1. Goal

Rewrite SwissIsoform from scratch as a modular 13-module pipeline, consolidating code from three repositories (`swissisoform`, `tiap`, `coTISja`) into a unified codebase with rich domain objects, comprehensive test fixtures, and a pipeline orchestrator. The rewrite is executed entirely by Claude Code sessions — a senior engineer session orchestrating parallel junior engineer sessions.

### What This Plan Covers

- Domain model design (Gene, TIS, Isoform hierarchy)
- Module interface contracts and column schemas
- Test architecture (two-tier: synthetic fast + biological correctness)
- Wave-based execution strategy for parallel Claude Code sessions
- Junior session prompt architecture

### What This Plan Does NOT Cover

- TIS caller replacement (Ribo-TISH wrapping for now; novel caller is a separate research track)
- Foundation model integration (Evo 2, Alpha Genome — reach goal)
- Data portal development
- coTISja (Paper 2) modeling work

### Reference Documents

| Document | Path | Purpose |
|----------|------|---------|
| Scientific overview | `swissisoform/docs/tis_projects_overview.md` | Biology, 2-paper structure, evidence scoring |
| Technical spec v2 | `swissisoform/docs/tis_technical_spec.md` | Module specs, code provenance, implementation phases |
| Execution plan | `swissisoform/docs/CLAUDE.md` | Step-by-step for sequential Claude Code execution |
| Gap analysis | `swissisoform/docs/swissisoform_review.md` | Missing modules, scoring expansion, tool recommendations |

### Source Repositories (Read-Only Reference)

| Repo | Path | What It Contributes |
|------|------|---------------------|
| swissisoform (current) | `/lab/barcheese01/mdiberna/swissisoform/` | BED parsing, translation, mutations, genome handling |
| TIAP | `/lab/barcheese01/mdiberna/tiap/` | Modular annotation pipeline (14 modules), pipeline.py pattern |
| coTISja | `/lab/barcheese01/smaffa/coTISja/` | Ribo-TISH filtering, Kozak annotation, expression normalization |

---

## 2. Architecture Decisions

### 2.1 Fresh Repository

The rewrite lives in a new repo at `/lab/barcheese01/mdiberna/swissisoform-v2/`. Rationale:

- **No legacy code in context.** Junior sessions won't waste tokens reading 16.5K lines of monolithic code (mutations.py alone is 4385 lines)
- **No import confusion.** Clean `src/swissisoform/` with only new modular code
- **Clean git history.** Starts with the domain model, not 20+ commits of legacy iteration
- **Smaller worktrees.** Each junior's worktree clone is fast and light

The old repos stay on disk as read-only reference. Juniors read source files via absolute paths.

### 2.2 Caller Strategy

The pipeline accepts clean Ribo-TISH output as input (Level 1: wrap existing caller). The domain model and module interfaces are caller-agnostic — a `TranslationInitiationSite` doesn't know or care whether it came from Ribo-TISH, a periodicity-aware caller, or a future multi-evidence model. Improving the caller (CAGE integration, novel statistical models) is a parallel research track that does not block this rewrite.

### 2.3 Domain Model over Flat DataFrames

The data is inherently hierarchical: a gene has many TIS, a TIS has expression across cell lines, each TIS maps to mutations. A flat DataFrame keyed on `variant_id` loses these relationships. The design uses rich domain objects that serialize to/from relational Parquet tables.

---

## 3. Domain Model

### 3.1 Core Classes

```python
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any


class ORFType(Enum):
    """Classification of the ORF produced by a TIS."""
    ANNOTATED = "Annotated"
    EXTENDED = "Extended"
    TRUNCATED = "Truncated"
    UORF = "uORF"
    UOORF = "uoORF"                    # upstream ORF overlapping annotated CDS
    INTERNAL_OUT_OF_FRAME = "internal_oof"
    THREE_UTR_ORF = "3utr_orf"
    ALT_ORF = "altORF"


@dataclass
class CellLineExpression:
    """Expression of a TIS in one cell line."""
    raw_count: int
    cpm: float
    p_value: float
    initiation_efficiency: float | None = None  # ratio to canonical TIS


@dataclass
class TranslationInitiationSite:
    """A single called TIS — the atomic unit of the pipeline.

    Each TIS represents a candidate translation start site detected by the
    caller (Ribo-TISH or successor). Modules annotate TIS objects by writing
    to the `annotations` dict under their MODULE_NAME key.
    """
    # Identity
    tis_id: str                          # "chr17:7676594:-" (genomic coordinate)
    gene_name: str
    transcript_id: str                   # ENST ID (MANE select or CAGE-resolved)
    chrom: str
    position: int
    strand: str                          # "+" or "-"
    start_codon: str                     # "AUG", "CUG", "GUG", "UUG", etc.
    orf_type: ORFType

    # Expression (per cell line)
    expression: dict[str, CellLineExpression] = field(default_factory=dict)

    # Protein products
    canonical_protein: str = ""          # full canonical AA sequence
    isoform_protein: str = ""            # full isoform AA sequence
    differential_sequence: str = ""      # the unique part
    shared_sequence: str = ""            # the overlapping part

    # Initiation context
    kozak_context: str | None = None     # -6 to +5 around start

    # Module annotations — each module writes to annotations[MODULE_NAME]
    annotations: dict[str, dict[str, Any]] = field(default_factory=dict)


@dataclass
class Gene:
    """A gene with its canonical reference and all detected TIS."""
    gene_name: str
    canonical_transcript_id: str
    canonical_protein: str
    tis_sites: list[TranslationInitiationSite] = field(default_factory=list)

    # Gene-level annotations (from gene_ref module)
    gene_annotations: dict[str, Any] = field(default_factory=dict)


@dataclass
class VariantAnnotation:
    """A single variant mapped to an isoform-specific region.

    Stored in a separate table with tis_id as foreign key.
    """
    tis_id: str                          # FK to TranslationInitiationSite
    source: str                          # "gnomad", "clinvar", "cosmic"
    variant_id: str                      # source-specific ID
    chrom: str
    position: int
    ref: str
    alt: str
    in_differential_region: bool         # does it fall in the unique part?
    metadata: dict[str, Any] = field(default_factory=dict)
```

### 3.2 Serialization Contract

Domain objects serialize to three Parquet tables:

| Table | Key | Content |
|-------|-----|---------|
| `tis_table.parquet` | `tis_id` | One row per TIS. Flattened annotations become prefixed columns (`biophysics_pI`, `conservation_diamond_score`, etc.). Expression becomes wide columns (`expr_HeLa_cpm`, `expr_K562_raw_count`, etc.) |
| `gene_table.parquet` | `gene_name` | One row per gene. Gene-level annotations from gene_ref module. |
| `variant_table.parquet` | `tis_id` + `variant_id` | One row per variant per TIS. Foreign key to tis_table. |

Round-trip functions in `src/swissisoform/io/parquet.py`:
- `tis_to_dataframe(sites: list[TranslationInitiationSite]) -> pd.DataFrame`
- `dataframe_to_tis(df: pd.DataFrame) -> list[TranslationInitiationSite]`
- `genes_to_dataframe(genes: list[Gene]) -> pd.DataFrame`

### 3.3 Design Rationale

- **`annotations` is a flat dict, not typed fields per module.** Adding a module never changes the domain class — it writes to a new key. Junior sessions can't break the model.
- **`ORFType` is an enum, not a string.** Classification errors caught at construction time.
- **Expression is a nested dict, not wide columns.** `expression["HeLa"].cpm` scales to 8+ cell lines without schema changes.
- **Variants live outside the TIS object.** The one-to-many relationship (many variants per TIS) doesn't belong inside the TIS dataclass — it's a separate table with FK.

---

## 4. Module Interface Contract

### 4.1 Protocol

```python
from typing import Protocol


class ModuleProtocol(Protocol):
    """Interface every annotation module must implement."""

    MODULE_NAME: str                     # e.g., "biophysics"
    OUTPUT_COLUMNS: list[str]            # columns produced when serialized
    SCOPE: str                           # "A", "B", "C", or "D"

    def __init__(self, config: PipelineConfig): ...

    def run(
        self, tis_sites: list[TranslationInitiationSite]
    ) -> list[TranslationInitiationSite]:
        """Annotate TIS sites, writing to tis.annotations[MODULE_NAME].

        Invariants:
        - len(output) == len(input)  — never drop a site
        - Every site gets annotations[MODULE_NAME] set (NaN/None for missing)
        - No mutation of fields outside annotations[MODULE_NAME]

        Note for Scope D modules (gene_ref, crossval): these operate at the
        gene level but still receive/return list[TIS]. They should group
        sites by gene_name internally, compute gene-level annotations once,
        and write the same gene-level result to all TIS for that gene.
        """
        ...
```

### 4.2 Scope Definitions

| Scope | Data Accessed | Examples |
|-------|---------------|----------|
| A — Differential region | `differential_sequence` vs. `shared_sequence` | Biophysics, motifs, VEP enrichment |
| B — Whole-protein paired | `canonical_protein` vs. `isoform_protein` | Localization, functional annotation, structure |
| C — Genomic position | `tis_id`, `chrom`, `position`, `strand` | TIS discovery, initiation context, PhyloP |
| D — Gene-level | `gene_name` only | Gene reference, cross-validation |

### 4.3 Module Registry

Each module's `OUTPUT_COLUMNS` are pre-defined. Junior sessions implement against this contract — they cannot invent or rename columns.

| Module | `MODULE_NAME` | Source | Scope | `OUTPUT_COLUMNS` |
|--------|---------------|--------|-------|-------------------|
| Paired comparison engine | `paired` | New | — | (utility, not a pipeline module) |
| Core Identity & Expression | `core_identity` | swissisoform | — | `variant_id`, `orf_type`, `canonical_protein_length`, `isoform_protein_length`, `differential_length_aa`, `in_frame`, `large_truncation_warning` |
| Initiation Context | `initiation_context` | coTISja | C | `kozak_context`, `kozak_hamming_major`, `kozak_hamming_partial`, `kozak_hamming_full`, `utr5_gc_content`, `upstream_aug_count`, `upstream_non_aug_count`, `gc_window_50bp`, `gc_window_250bp` |
| Biophysics | `biophysics` | TIAP | A | `biophysics_pI`, `biophysics_gravy`, `biophysics_instability`, `biophysics_disorder_fraction`, `biophysics_llps_score`, `biophysics_mol_weight`, `biophysics_charge_ph7`, `biophysics_aromaticity`, `biophysics_helix_fraction`, `biophysics_sheet_fraction`, `biophysics_coil_fraction`, `biophysics_extinction_coeff`, `biophysics_half_life`, `biophysics_aliphatic_index`, `biophysics_boman_index`, `biophysics_hydrophobic_moment`, `biophysics_flexibility`, `biophysics_isoelectric_point`, `biophysics_aa_composition_chi2`, `biophysics_aa_composition_pvalue` |
| Motifs | `motifs` | TIAP | A | `motifs_cdk_count`, `motifs_atm_count`, `motifs_14_3_3_count`, `motifs_eb1_count`, `motifs_rgg_count`, `motifs_sh3_count`, `motifs_heme_count`, `motifs_pip_count`, `motifs_apim_count`, `motifs_znf_count`, `motifs_ring_count`, `motifs_nls_count`, `motifs_nes_count`, `motifs_sumo_count`, `motifs_phospho_count`, `motifs_total_slim_count`, `motifs_total_slim_density` |
| Functional Annotation | `functional` | New (stub) | B | `functional_domains_canonical`, `functional_domains_isoform`, `functional_domains_gained`, `functional_domains_lost`, `functional_domains_disrupted`, `functional_signal_peptide_canonical`, `functional_signal_peptide_isoform`, `functional_signal_peptide_changed`, `functional_tm_topology_canonical`, `functional_tm_topology_isoform`, `functional_tm_topology_changed`, `functional_targeting_canonical`, `functional_targeting_isoform`, `functional_targeting_changed` |
| Structure | `structure` | New (stub) | A+B | `structure_plddt_diffseq`, `structure_plddt_canonical`, `structure_plddt_isoform`, `structure_plddt_delta`, `structure_tm_score`, `structure_rmsd_shared`, `structure_extension_contacts`, `structure_method`, `structure_ss_helix_diff`, `structure_ss_sheet_diff`, `structure_ss_coil_diff` |
| Localization | `localization` | TIAP | B | `localization_canonical`, `localization_isoform`, `localization_changed`, `localization_canonical_confidence`, `localization_isoform_confidence` |
| Variant Effect Prediction | `vep` | New (stub) | A+B | `vep_am_mean_pathogenicity_unique`, `vep_am_mean_pathogenicity_shared`, `vep_am_enrichment_ratio`, `vep_am_pathogenic_count_unique`, `vep_esm1b_damaging_isoform_specific`, `vep_esm1b_mean_llr_unique`, `vep_esm1b_mean_llr_shared` |
| Conservation | `conservation` | TIAP | A+C | `conservation_diamond_score`, `conservation_diamond_label`, `conservation_tblastn_primate_intact`, `conservation_tblastn_mammal_intact`, `conservation_phylop_tis`, `conservation_phylop_kozak_mean`, `conservation_phylop_unique_mean`, `conservation_phylocsf_score` |
| Cross-Validation | `crossval` | TIAP | D | `crossval_chen_weissman`, `crossval_ingolia_mouse`, `crossval_fedorova_extension`, `crossval_qti_stress`, `crossval_kagan_immunity`, `crossval_fedorova_proteomics_count`, `crossval_pepquery2_validated` |
| Gene Reference | `gene_ref` | TIAP | D | `generef_uniprot`, `generef_hpa`, `generef_depmap`, `generef_gygi`, `generef_leonetti`, `generef_ly`, `generef_chx`, `generef_omim` |
| Clinical Variation | `clinical` | swissisoform | A+C | `clinical_gnomad_gene_missense_count`, `clinical_gnomad_gene_lof_count`, `clinical_gnomad_gene_constraint_oe`, `clinical_clinvar_gene_pathogenic_count`, `clinical_cosmic_gene_mutation_count`, `clinical_gnomad_unique_variant_density`, `clinical_gnomad_shared_variant_density`, `clinical_clinvar_unique_pathogenic`, `clinical_cosmic_unique_somatic` |
| Evidence Scoring | `scoring` | TIAP | — | `scoring_existence_score`, `scoring_functional_impact_score`, `scoring_combined_score`, `scoring_confidence_tier`, `scoring_evidence_flags` |

---

## 5. Test Architecture

### 5.1 Two Tiers

**Tier 1 — Synthetic unit tests (fast, < 5s total, no I/O):**
- Pre-computed protein sequences for 8 TIS across 3 genes
- Hardcoded constants — not computed at test time
- Cover: extension (+ strand AUG, + strand CUG, - strand), truncation (+ strand, - strand), uORF, uoORF, annotated
- Every module test: columns present, no rows lost, known-value spot checks
- This is what overnight junior sessions run

**Tier 2 — Biological correctness tests (slow, requires genome data):**
- Run against real GENCODE genome + GTF
- Validate translation against known proteins (UniProt cross-check)
- Test negative strand handling, non-AUG codons, stop codon inclusion
- Marked `@pytest.mark.slow`
- Run during Wave 1 and integration only

### 5.2 Fixture Design

The synthetic fixtures must exercise the biological edge cases that cause real bugs:

```python
# conftest.py — synthetic TIS fixtures

SYNTHETIC_GENES = {
    "TESTGENE_POS": {
        # Positive strand gene with known canonical protein
        "canonical_transcript": "ENST_TEST_POS",
        "chrom": "chr1",
        "strand": "+",
        "canonical_protein": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS*",
    },
    "TESTGENE_NEG": {
        # Negative strand gene — catches reverse complement bugs
        "canonical_transcript": "ENST_TEST_NEG",
        "chrom": "chr2",
        "strand": "-",
        "canonical_protein": "MDLSALREVELIQNMHRKEVTHPVFD*",
    },
    "TESTGENE_MULTI": {
        # Gene with multiple TIS — tests per-gene grouping
        "canonical_transcript": "ENST_TEST_MULTI",
        "chrom": "chr3",
        "strand": "+",
        "canonical_protein": "MSSGNAKIGHPAPNFKATAVMHNEFIASK*",
    },
}

SYNTHETIC_TIS = [
    # 1. Annotated — canonical AUG, should match canonical protein exactly
    TranslationInitiationSite(
        tis_id="chr1:1000:+", gene_name="TESTGENE_POS",
        transcript_id="ENST_TEST_POS", chrom="chr1", position=1000,
        strand="+", start_codon="AUG", orf_type=ORFType.ANNOTATED,
        canonical_protein="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS*",
        isoform_protein="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS*",
        differential_sequence="", shared_sequence="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS*",
        expression={"HeLa": CellLineExpression(1000, 50.0, 0.001, 1.0),
                    "K562": CellLineExpression(800, 40.0, 0.002, 1.0)},
    ),
    # 2. Extension (+ strand, AUG) — adds N-terminal prefix
    TranslationInitiationSite(
        tis_id="chr1:950:+", gene_name="TESTGENE_POS",
        transcript_id="ENST_TEST_POS", chrom="chr1", position=950,
        strand="+", start_codon="AUG", orf_type=ORFType.EXTENDED,
        canonical_protein="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS*",
        isoform_protein="MRGSHHHHHGSMEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS*",
        differential_sequence="MRGSHHHHHGS",
        shared_sequence="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS*",
        expression={"HeLa": CellLineExpression(250, 12.5, 0.005, 0.25),
                    "K562": CellLineExpression(200, 10.0, 0.008, 0.25)},
    ),
    # 3. Extension (+ strand, CUG) — non-AUG start, leucine instead of methionine
    TranslationInitiationSite(
        tis_id="chr1:940:+", gene_name="TESTGENE_POS",
        transcript_id="ENST_TEST_POS", chrom="chr1", position=940,
        strand="+", start_codon="CUG", orf_type=ORFType.EXTENDED,
        canonical_protein="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS*",
        isoform_protein="LRRPPAGAMEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS*",
        differential_sequence="LRRPPAGA",
        shared_sequence="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS*",
        expression={"HeLa": CellLineExpression(100, 5.0, 0.01, 0.1),
                    "K562": CellLineExpression(50, 2.5, 0.05, 0.06)},
    ),
    # 4. Extension (- strand) — catches reverse complement bugs
    TranslationInitiationSite(
        tis_id="chr2:5000:-", gene_name="TESTGENE_NEG",
        transcript_id="ENST_TEST_NEG", chrom="chr2", position=5000,
        strand="-", start_codon="AUG", orf_type=ORFType.EXTENDED,
        canonical_protein="MDLSALREVELIQNMHRKEVTHPVFD*",
        isoform_protein="MAGTKLMDLSALREVELIQNMHRKEVTHPVFD*",
        differential_sequence="MAGTKLM",
        shared_sequence="DLSALREVELIQNMHRKEVTHPVFD*",
        expression={"HeLa": CellLineExpression(500, 25.0, 0.001, 0.5),
                    "K562": CellLineExpression(400, 20.0, 0.003, 0.5)},
    ),
    # 5. Truncation (+ strand) — shorter than canonical
    TranslationInitiationSite(
        tis_id="chr1:1030:+", gene_name="TESTGENE_POS",
        transcript_id="ENST_TEST_POS", chrom="chr1", position=1030,
        strand="+", start_codon="AUG", orf_type=ORFType.TRUNCATED,
        canonical_protein="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS*",
        isoform_protein="MEPPLSQETFSDLWKLLPENNVLSPLPS*",
        differential_sequence="EEPQSDPSV",  # the lost N-terminal region
        shared_sequence="MEPPLSQETFSDLWKLLPENNVLSPLPS*",
        expression={"HeLa": CellLineExpression(150, 7.5, 0.01, 0.15),
                    "K562": CellLineExpression(100, 5.0, 0.02, 0.13)},
    ),
    # 6. Truncation (- strand)
    TranslationInitiationSite(
        tis_id="chr2:4970:-", gene_name="TESTGENE_NEG",
        transcript_id="ENST_TEST_NEG", chrom="chr2", position=4970,
        strand="-", start_codon="AUG", orf_type=ORFType.TRUNCATED,
        canonical_protein="MDLSALREVELIQNMHRKEVTHPVFD*",
        isoform_protein="MEVELIQNMHRKEVTHPVFD*",
        differential_sequence="DLSALR",  # lost region
        shared_sequence="MEVELIQNMHRKEVTHPVFD*",
        expression={"HeLa": CellLineExpression(300, 15.0, 0.005, 0.3),
                    "K562": CellLineExpression(200, 10.0, 0.01, 0.25)},
    ),
    # 7. uORF — upstream ORF, entirely in 5' UTR
    TranslationInitiationSite(
        tis_id="chr3:800:+", gene_name="TESTGENE_MULTI",
        transcript_id="ENST_TEST_MULTI", chrom="chr3", position=800,
        strand="+", start_codon="AUG", orf_type=ORFType.UORF,
        canonical_protein="MSSGNAKIGHPAPNFKATAVMHNEFIASK*",
        isoform_protein="MPKLQRST*",  # short uORF product
        differential_sequence="MPKLQRST*",  # entire product is unique
        shared_sequence="",  # no overlap with canonical CDS
        expression={"HeLa": CellLineExpression(30, 1.5, 0.05, 0.03),
                    "K562": CellLineExpression(10, 0.5, 0.1, 0.01)},
    ),
    # 8. uoORF — upstream ORF overlapping canonical CDS
    TranslationInitiationSite(
        tis_id="chr3:850:+", gene_name="TESTGENE_MULTI",
        transcript_id="ENST_TEST_MULTI", chrom="chr3", position=850,
        strand="+", start_codon="AUG", orf_type=ORFType.UOORF,
        canonical_protein="MSSGNAKIGHPAPNFKATAVMHNEFIASK*",
        isoform_protein="MFLGRTSSSGNAKIGHPAPNFKATAVMHNEFIASK*",
        differential_sequence="MFLGRT",
        shared_sequence="SSSGNAKIGHPAPNFKATAVMHNEFIASK*",
        expression={"HeLa": CellLineExpression(20, 1.0, 0.08, 0.02),
                    "K562": CellLineExpression(5, 0.25, 0.2, 0.005)},
    ),
]
```

### 5.3 Per-Module Test Contract

Every module test file follows this pattern:

```python
class TestModuleX:
    def test_output_columns_present(self, config, synthetic_tis):
        result = ModuleX(config).run(synthetic_tis)
        for site in result:
            assert ModuleX.MODULE_NAME in site.annotations
            for col in ModuleX.OUTPUT_COLUMNS:
                key = col.removeprefix(f"{ModuleX.MODULE_NAME}_")
                assert key in site.annotations[ModuleX.MODULE_NAME]

    def test_no_sites_lost(self, config, synthetic_tis):
        result = ModuleX(config).run(synthetic_tis)
        assert len(result) == len(synthetic_tis)

    def test_handles_empty_differential(self, config, synthetic_tis):
        """Annotated isoforms (empty differential_sequence) don't crash."""
        result = ModuleX(config).run(synthetic_tis)
        annotated = [s for s in result if s.orf_type == ORFType.ANNOTATED]
        assert len(annotated) > 0  # fixture includes annotated

    def test_known_value(self, config):
        """Spot-check a known computation."""
        ...  # module-specific
```

### 5.4 pytest Configuration

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
markers = [
    "slow: integration tests (full pipeline, real genome data)",
    "gpu: requires GPU (structure prediction, embeddings)",
    "network: requires network access (API calls, downloads)",
]
addopts = "-m 'not slow and not gpu and not network' --tb=short -q"
```

---

## 6. Paired Comparison Engine

Shared utility used by modules with Scope A or B. Lives in `src/swissisoform/compare/paired.py`.

```python
class PairedComparison:
    """Canonical-vs-isoform delta logic shared across modules."""

    @staticmethod
    def compare_categorical(canonical: str, isoform: str) -> dict:
        """For localization, signal peptide, targeting: did it change?"""
        return {"changed": canonical != isoform,
                "canonical": canonical, "isoform": isoform}

    @staticmethod
    def compare_sets(canonical_set: set, isoform_set: set) -> dict:
        """For domains: what was gained, lost, disrupted?"""
        return {"gained": list(isoform_set - canonical_set),
                "lost": list(canonical_set - isoform_set),
                "shared": list(canonical_set & isoform_set)}

    @staticmethod
    def compare_scalar_regions(
        unique_value: float, shared_value: float
    ) -> dict:
        """For Scope A: enrichment of metric in unique vs shared region."""
        ratio = unique_value / shared_value if shared_value != 0 else float("inf")
        return {"unique": unique_value, "shared": shared_value,
                "ratio": ratio, "enriched": ratio > 1.0}

    @staticmethod
    def compare_structures(
        canonical_metrics: dict, isoform_metrics: dict
    ) -> dict:
        """Structure comparison metrics."""
        return {
            "plddt_delta": (isoform_metrics.get("plddt", 0)
                           - canonical_metrics.get("plddt", 0)),
            "tm_score": isoform_metrics.get("tm_score"),
            "rmsd_shared": isoform_metrics.get("rmsd_shared"),
            "extension_contacts": isoform_metrics.get("extension_contacts", 0),
        }
```

---

## 7. Pipeline Configuration

```python
@dataclass
class PipelineConfig:
    """Configuration for the full pipeline."""
    cell_lines: list[str] = field(
        default_factory=lambda: [
            "HeLa", "K562", "U2OS", "RPE1_Async",
            "RPE1_Que", "RPE1_Sen", "iPSC", "HFF",
        ]
    )
    genome_fasta: Path | None = None
    gtf_path: Path | None = None
    output_dir: Path = Path("results")

    # Module-specific config (modules check for their key)
    conservation: ConservationConfig | None = None
    structure: StructureConfig | None = None
    scoring: ScoringConfig | None = None
    clinical: ClinicalConfig | None = None

@dataclass
class ConservationConfig:
    diamond_db: Path | None = None
    tblastn_db: Path | None = None
    phylop_bigwig: Path | None = None

@dataclass
class StructureConfig:
    method: str = "chai1"               # chai1, boltz2, af3
    device: str = "cuda"
    batch_size: int = 4

@dataclass
class ScoringConfig:
    min_cell_lines: int = 3
    existence_high_threshold: int = 5
    functional_high_threshold: int = 3
    truncation_max_aa: int = 200

@dataclass
class ClinicalConfig:
    gnomad_api_url: str = "https://gnomad.broadinstitute.org/api"
    clinvar_email: str = ""
    cosmic_db: Path | None = None
```

---

## 8. Wave Execution Strategy

### 8.1 Overview

```
Wave 1 — Foundation (interactive + 2 juniors)
    ├── Senior: scaffold, domain model, fixtures, paired engine
    ├── Junior W1-A: core_identity.py (hardest port)
    └── Junior W1-B: initiation_context.py (from coTISja)
         ↓
    REVIEW GATE: domain model validated with real porting
         ↓
Wave 2 — Parallel module ports (overnight, 8 juniors)
    ├── W2-A: biophysics.py      (TIAP, low complexity)
    ├── W2-B: motifs.py           (TIAP, low complexity)
    ├── W2-C: conservation.py     (TIAP, medium — stubs)
    ├── W2-D: crossval.py         (TIAP, medium — 5 matchers)
    ├── W2-E: gene_ref.py         (TIAP, low complexity)
    ├── W2-F: localization.py     (TIAP, low complexity)
    ├── W2-G: clinical.py         (swissisoform, HIGH — biggest refactor)
    └── W2-H: scoring.py          (TIAP, medium — 13→17 criteria)
         ↓
    REVIEW GATE: senior merges, full Tier 1 test suite
         ↓
Wave 3 — Integration (interactive + 2-3 juniors)
    ├── Senior: pipeline orchestrator (pipeline.py)
    ├── Junior W3-A: stub modules (functional, VEP, structure)
    ├── Junior W3-B: integration tests
    └── Junior W3-C: Tier 2 biological correctness tests
         ↓
    FINAL GATE: all tests green, old code archived
```

### 8.2 Wave 1 — Foundation

**Senior engineer (interactive, with user review):**

1. `init-codebase` scaffold at `/lab/barcheese01/mdiberna/swissisoform-v2/`
2. `src/swissisoform/models.py` — domain model (Section 3)
3. `src/swissisoform/modules/base.py` — module protocol + `validate_module_output()` helper
4. `src/swissisoform/compare/paired.py` — paired comparison engine (Section 6)
5. `tests/conftest.py` — Tier 1 synthetic fixtures (Section 5.2)
6. `src/swissisoform/io/parquet.py` — serialization round-trip
7. `src/swissisoform/config.py` — `PipelineConfig` and sub-configs (Section 7)
8. Commit: `chore: scaffold codebase with domain model and test fixtures`

**Junior W1-A: `core_identity.py`** (can run overnight after senior finishes)
- Read from: `/lab/barcheese01/mdiberna/swissisoform/src/swissisoform/alternative_isoforms.py` and `translation.py`
- Port: BED parsing, ORF classification, protein extraction, differential sequence extraction
- Write: `src/swissisoform/modules/core_identity.py` + `tests/test_core_identity.py`
- Run both Tier 1 (synthetic) and Tier 2 (real genome) tests
- Complexity: HIGH — translation logic is inherently difficult

**Junior W1-B: `initiation_context.py`**
- Read from: `/lab/barcheese01/smaffa/coTISja/src/scripts/analysis_pipeline_helpers.py`
- Port: Kozak extraction, Hamming distance, GC content, upstream counts
- Write: `src/swissisoform/modules/initiation_context.py` + `tests/test_initiation_context.py`
- Tier 1 tests only
- Complexity: LOW-MEDIUM

**Review gate:** User verifies domain model works with real ported code. Core identity module produces correct translations for known genes on both strands.

### 8.3 Wave 2 — Parallel Module Ports

All 8 juniors run in parallel git worktrees. Each gets the same preamble (Section 9) plus module-specific instructions.

| Junior | Module | Source Path | Complexity | Notes |
|--------|--------|-------------|------------|-------|
| W2-A | `biophysics.py` | `/lab/barcheese01/mdiberna/tiap/src/tiap/modules/biophysical.py` | Low | Direct port + AA composition contrast |
| W2-B | `motifs.py` | `/lab/barcheese01/mdiberna/tiap/src/tiap/modules/motifs.py` | Low | Direct port of 15 SLiM families |
| W2-C | `conservation.py` | `/lab/barcheese01/mdiberna/tiap/src/tiap/modules/conservation.py` | Medium | Port DIAMOND + tBLASTn/PhyloP stubs |
| W2-D | `crossval.py` | `/lab/barcheese01/mdiberna/tiap/src/tiap/modules/crossval.py` + `massspec.py` | Medium | 5 dataset matchers + digestion logic |
| W2-E | `gene_ref.py` | `/lab/barcheese01/mdiberna/tiap/src/tiap/modules/gene_ref.py` | Low | Direct port of 8 annotators |
| W2-F | `localization.py` | `/lab/barcheese01/mdiberna/tiap/src/tiap/modules/localization.py` | Low | Already does paired comparison |
| W2-G | `clinical.py` | `/lab/barcheese01/mdiberna/swissisoform/src/swissisoform/mutations.py` | **High** | 4385-line refactor into clean interface |
| W2-H | `scoring.py` | `/lab/barcheese01/mdiberna/tiap/src/tiap/modules/scoring.py` | Medium | Expand 13→17, dual-axis framework |

**Review gate:** Senior merges all branches, runs `pytest` (full Tier 1 suite), resolves conflicts.

### 8.4 Wave 3 — Integration

**Senior: Pipeline orchestrator** (`src/swissisoform/pipeline.py`)
- Wire all modules into `run_cpu()`, `run_gpu()`, `merge_and_score()`
- Master table assembly

**Junior W3-A: Stub modules**
- `modules/functional.py` (InterPro, SignalP, DeepTMHMM, TargetP)
- `modules/vep.py` (AlphaMissense, ESM1b)
- `modules/structure.py` (Chai-1/Boltz-2)
- Each: correct interface, `OUTPUT_COLUMNS`, NaN fills, passing tests

**Junior W3-B: Integration tests**
- `tests/test_integration.py` marked `@pytest.mark.slow`
- Full pipeline on synthetic data: BED input → all modules → scored output
- Round-trip serialization test

**Junior W3-C: Tier 2 biological correctness**
- Expanded tests against real HeLa data
- UniProt cross-check for known proteins
- Regression test against existing pipeline results

**Final gate:** All tests green. Old monolithic code archived.

---

## 9. Junior Session Prompt Template

### 9.1 Preamble (Identical for All Wave 2 Juniors)

```
You are implementing a module for the SwissIsoform v2 pipeline.

PROJECT: SwissIsoform annotates protein isoforms from alternate translation
initiation sites (TIS). You are building one module of a 13-module pipeline.

REPO: /lab/barcheese01/mdiberna/swissisoform-v2/

DOMAIN MODEL: Read src/swissisoform/models.py for the TranslationInitiationSite,
Gene, and VariantAnnotation dataclasses. Your module receives and returns
list[TranslationInitiationSite].

MODULE PROTOCOL: Read src/swissisoform/modules/base.py. Your module must:
- Define MODULE_NAME, OUTPUT_COLUMNS, SCOPE as class attributes
- Implement run(tis_sites: list[TranslationInitiationSite]) -> list[TranslationInitiationSite]
- Write ONLY to site.annotations[MODULE_NAME] for each site
- Never drop sites (len(output) == len(input))
- Use NaN/None for values that can't be computed

TEST FIXTURES: Read tests/conftest.py for the synthetic_tis fixture.
Your tests must use this fixture.

PAIRED COMPARISON: If your module does canonical-vs-isoform comparison,
use src/swissisoform/compare/paired.py — do not implement your own delta logic.

CODE STYLE: ruff (line length 100), Google-style docstrings, type annotations.

COMMIT: Single commit with message "feat({module-name}): {description}"
```

### 9.2 Module-Specific Section (Example: W2-A Biophysics)

```
MODULE: biophysics (Module 3a)
SCOPE: A (operates on differential_sequence vs. shared_sequence)

SOURCE: Read the following file and port its logic:
  /lab/barcheese01/mdiberna/tiap/src/tiap/modules/biophysical.py

Also read for context:
  /lab/barcheese01/mdiberna/tiap/src/tiap/modules/differential.py

WHAT TO PORT:
- 18 biophysical properties (pI, GRAVY, disorder fraction, instability,
  LLPS score, molecular weight, charge at pH 7, aromaticity, helix/sheet/coil
  fractions, extinction coefficient, half-life, aliphatic index, Boman index,
  hydrophobic moment, flexibility, isoelectric point)
- Compute on differential_sequence (the unique region)

WHAT'S NEW (not in TIAP):
- AA composition contrast: per-residue frequencies in differential_sequence
  vs. shared_sequence, with chi-squared test
- Output: biophysics_aa_composition_chi2, biophysics_aa_composition_pvalue

OUTPUT_COLUMNS (exact names, all prefixed with "biophysics_"):
  biophysics_pI, biophysics_gravy, biophysics_instability,
  biophysics_disorder_fraction, biophysics_llps_score, biophysics_mol_weight,
  biophysics_charge_ph7, biophysics_aromaticity, biophysics_helix_fraction,
  biophysics_sheet_fraction, biophysics_coil_fraction,
  biophysics_extinction_coeff, biophysics_half_life,
  biophysics_aliphatic_index, biophysics_boman_index,
  biophysics_hydrophobic_moment, biophysics_flexibility,
  biophysics_isoelectric_point, biophysics_aa_composition_chi2,
  biophysics_aa_composition_pvalue

WRITE:
  src/swissisoform/modules/biophysics.py
  tests/test_biophysics.py

TESTS MUST INCLUDE:
- test_output_columns_present
- test_no_sites_lost
- test_handles_empty_differential (annotated isoforms)
- test_known_pI (compute pI for "MRGSHHHHHGS", verify against known value)
- test_known_gravy (compute GRAVY for a short hydrophobic sequence)
- test_aa_composition_contrast (sequences with different compositions)

DONE WHEN:
- pytest tests/test_biophysics.py -v passes (< 5s)
- All OUTPUT_COLUMNS present in annotations["biophysics"]
- Single commit made

CONSTRAINTS:
- Do NOT modify models.py, conftest.py, base.py, or any other module
- Do NOT add dependencies not already in pyproject.toml
- Do NOT import from the old swissisoform code — port, don't wrap
```

---

## 10. Target Codebase Structure

```
swissisoform-v2/
├── pyproject.toml
├── CLAUDE.md
├── config.yaml
├── src/swissisoform/
│   ├── __init__.py
│   ├── models.py                    # Domain model (Gene, TIS, VariantAnnotation)
│   ├── config.py                    # PipelineConfig + sub-configs
│   ├── pipeline.py                  # Pipeline orchestrator
│   ├── modules/
│   │   ├── __init__.py
│   │   ├── base.py                  # ModuleProtocol + validate_module_output()
│   │   ├── core_identity.py         # Module 1: BED parsing, classification, proteins
│   │   ├── initiation_context.py    # Module 2: Kozak, Hamming, GC, upstream counts
│   │   ├── biophysics.py            # Module 3a: 18 biophysical properties
│   │   ├── motifs.py                # Module 3b: 15 SLiM families
│   │   ├── functional.py            # Module 4: InterPro, SignalP, etc. (stub)
│   │   ├── structure.py             # Module 5: Chai-1/Boltz-2 (stub)
│   │   ├── localization.py          # Module 6: DeepLoc paired comparison
│   │   ├── vep.py                   # Module 7: AlphaMissense, ESM1b (stub)
│   │   ├── conservation.py          # Module 8: DIAMOND + stubs
│   │   ├── crossval.py              # Module 9: 5 dataset matchers
│   │   ├── gene_ref.py              # Module 10: 8 database annotators
│   │   ├── clinical.py              # Module 11: gnomAD/ClinVar/COSMIC
│   │   └── scoring.py               # Module 12: dual-axis 17-criteria
│   ├── compare/
│   │   ├── __init__.py
│   │   └── paired.py                # Canonical-vs-isoform delta logic
│   └── io/
│       ├── __init__.py
│       └── parquet.py               # Serialization round-trip
├── tests/
│   ├── conftest.py                  # Shared synthetic fixtures
│   ├── test_paired.py
│   ├── test_core_identity.py
│   ├── test_initiation_context.py
│   ├── test_biophysics.py
│   ├── test_motifs.py
│   ├── test_functional.py
│   ├── test_structure.py
│   ├── test_localization.py
│   ├── test_vep.py
│   ├── test_conservation.py
│   ├── test_crossval.py
│   ├── test_gene_ref.py
│   ├── test_clinical.py
│   ├── test_scoring.py
│   └── test_integration.py          # @pytest.mark.slow
├── data/                             # Symlinks or paths to shared data
├── scripts/
│   └── run_pipeline.py              # CLI entry point
└── structures/                       # PDB/mmCIF outputs (future)
```

---

## 11. Risk Register

| Risk | Impact | Mitigation |
|------|--------|------------|
| Domain model doesn't capture a relationship needed by a module | All Wave 2 juniors build on wrong foundation | Wave 1 validates with core_identity (hardest module) before fan-out |
| Junior sessions diverge on interpretation of `OUTPUT_COLUMNS` | Merge conflicts, inconsistent serialization | Exact column names pre-defined in this spec; module protocol enforces naming |
| Translation bugs on negative strand / non-AUG codons | Incorrect protein sequences propagate to all annotation modules | Tier 2 tests run in Wave 1 with known genes on both strands |
| `clinical.py` refactor too complex for single junior session | Incomplete or broken mutation handling | Give W2-G extra context in prompt; accept stubs for isoform-level features |
| Worktree merge conflicts | Integration delay | Each junior writes to different files; only conftest.py and `__init__.py` are shared (read-only) |
| conftest.py synthetic sequences don't exercise a real edge case | A module passes tests but fails on real data | Tier 2 tests in Wave 3 catch this; fixtures designed for known failure modes |
