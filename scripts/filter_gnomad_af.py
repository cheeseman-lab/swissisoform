#!/usr/bin/env python3
"""Filter ClinVar/COSMIC/clinical mutations by removing those that overlap with common gnomAD variants.

This script filters mutations by removing those that match gnomAD variants at or above a specified
allele frequency threshold. It operates on an existing isoform_level_results file (typically with
MANE annotations) and updates only the count columns while preserving all metadata.

Two Operating Modes:

1. Self-Contained Mode:
   Dataset contains gnomAD, ClinVar, and COSMIC variants.
   - Filter ClinVar against gnomAD with --clinvar-threshold (optional)
   - Filter COSMIC against gnomAD with --cosmic-threshold (optional)
   - Can specify one or both thresholds
   - gnomAD variants are RETAINED (used as filter reference only)

2. Cross-Dataset Mode:
   Test dataset contains clinical mutations (BCH/MSK).
   Reference dataset provides gnomAD variants.
   - Filter test mutations against gnomAD with --threshold

Usage:
    # Self-contained mode - both thresholds (using mane file)
    python scripts/filter_gnomad_af.py results/hela \\
        --mane-file results/hela/mutations/isoform_level_results_mane.csv \\
        --clinvar-threshold 0.01 --cosmic-threshold 0.001

    # Cross-dataset mode (using mane file)
    python scripts/filter_gnomad_af.py results/hela_bch \\
        --mane-file results/hela_bch/mutations/isoform_level_results_mane.csv \\
        --gnomad-source results/hela --threshold 0.01

    # Dry run for statistics only
    python scripts/filter_gnomad_af.py results/hela \\
        --mane-file results/hela/mutations/isoform_level_results_mane.csv \\
        --clinvar-threshold 0.01 --dry-run
"""

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Impact types in standard order
IMPACT_TYPES = [
    "missense_variant",
    "nonsense_variant",
    "frameshift_variant",
    "inframe_deletion",
    "inframe_insertion",
    "synonymous_variant",
]

# Sources in standard order
SOURCES = ["clinvar", "gnomad", "cosmic"]


def find_mutation_csvs(dataset_dir: Path) -> List[Path]:
    """Find all mutation CSV files in the dataset directory.

    Returns only original mutation files (ending in _mutations.csv),
    excluding any previously filtered files.
    """
    mutations_dir = dataset_dir / "mutations"
    if not mutations_dir.exists():
        raise FileNotFoundError(f"Mutations directory not found: {mutations_dir}")

    csv_files = list(mutations_dir.glob("*/*/*.csv"))
    # Filter to only original mutation files (not filtered versions)
    csv_files = [
        f
        for f in csv_files
        if f.name.endswith("_mutations.csv") and "_filtered_" not in f.name
    ]
    return sorted(csv_files)


def parse_feature_info_from_path(csv_path: Path) -> Optional[Dict]:
    """Extract feature info from mutation CSV path.

    Path format: mutations/{gene}/{transcript}/{transcript}_{feature_type}_{codon}_{start}_{end}_mutations.csv
    """
    filename = csv_path.stem  # Remove .csv
    filename = filename.replace("_mutations", "")  # Remove _mutations suffix

    parts = filename.split("_")
    # Expected format: ENST00000209873.8_extension_CTG_53321463_53321480
    # transcript_id_feature_type_codon_start_end

    if len(parts) < 5:
        logger.warning(f"Cannot parse feature info from: {csv_path}")
        return None

    transcript_id = parts[0]
    feature_type = parts[1]
    codon = parts[2]
    feature_start = int(parts[3])
    feature_end = int(parts[4])

    gene_name = csv_path.parent.parent.name

    return {
        "gene_name": gene_name,
        "transcript_id": transcript_id,
        "feature_type": feature_type,
        "feature_start": feature_start,
        "feature_end": feature_end,
        "codon": codon,
        "feature_id": f"{feature_type}_{codon}_{feature_start}_{feature_end}",
    }


def build_gnomad_lookup(dataset_dir: Path) -> Dict[Tuple[str, str, str, str], float]:
    """Build lookup dictionary of ALL gnomAD variants from dataset.

    Args:
        dataset_dir: Path to dataset containing gnomAD variants

    Returns:
        Dict mapping (chr, pos, ref, alt) to gnomAD allele frequency
    """
    logger.info(f"Building gnomAD lookup from {dataset_dir}...")
    csv_files = find_mutation_csvs(dataset_dir)
    logger.info(f"Found {len(csv_files)} mutation files")

    gnomad_lookup = {}
    total_gnomad = 0

    for csv_path in csv_files:
        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            logger.warning(f"Error reading {csv_path}: {e}")
            continue

        if df.empty or "source" not in df.columns:
            continue

        # Get gnomAD variants only
        gnomad_df = df[df["source"].str.lower() == "gnomad"].copy()
        total_gnomad += len(gnomad_df)

        if gnomad_df.empty:
            continue

        af_col = "gnomad_allele_frequency"

        # Add to lookup (include ALL gnomAD variants, no pre-filtering)
        for _, row in gnomad_df.iterrows():
            key = (
                str(row.get("chromosome", "")),
                str(row.get("position", "")),
                str(row.get("reference", "")),
                str(row.get("alternate", "")),
            )
            af_value = row.get(af_col, None)
            if pd.notna(af_value):
                gnomad_lookup[key] = float(af_value)
            else:
                # Store None for variants without AF (won't match threshold filter)
                gnomad_lookup[key] = None

    logger.info(
        f"gnomAD lookup: {total_gnomad} total variants, {len(gnomad_lookup)} unique positions"
    )
    return gnomad_lookup


def get_filter_suffix(
    clinvar_threshold: Optional[float] = None,
    cosmic_threshold: Optional[float] = None,
    threshold: Optional[float] = None,
) -> str:
    """Generate filename suffix from thresholds.

    Examples:
        clinvar_threshold=0.01, cosmic_threshold=0.001 -> "clinvar0.01_cosmic0.001"
        clinvar_threshold=0.01 -> "clinvar0.01"
        threshold=0.01 -> "gnomad0.01"
    """
    parts = []

    if threshold is not None:
        # Cross-dataset mode
        parts.append(f"gnomad{threshold}")
    else:
        # Self-contained mode
        if clinvar_threshold is not None:
            parts.append(f"clinvar{clinvar_threshold}")
        if cosmic_threshold is not None:
            parts.append(f"cosmic{cosmic_threshold}")

    return "_".join(parts)


def filter_by_gnomad_overlap(
    df: pd.DataFrame,
    gnomad_lookup: Dict[Tuple[str, str, str, str], float],
    clinvar_threshold: Optional[float] = None,
    cosmic_threshold: Optional[float] = None,
    cross_dataset_threshold: Optional[float] = None,
) -> Tuple[pd.DataFrame, Dict]:
    """Filter mutations by removing those matching common gnomAD variants.

    Self-contained mode (clinvar_threshold and/or cosmic_threshold set):
        - gnomAD variants: ALWAYS KEEP
        - ClinVar variants: REMOVE if matches gnomAD with AF >= clinvar_threshold
        - COSMIC variants: REMOVE if matches gnomAD with AF >= cosmic_threshold
        - Other sources: KEEP

    Cross-dataset mode (cross_dataset_threshold set):
        - ALL mutations: REMOVE if matches gnomAD with AF >= threshold

    Args:
        df: DataFrame with mutations
        gnomad_lookup: Dict mapping (chr, pos, ref, alt) to gnomAD AF
        clinvar_threshold: AF threshold for filtering ClinVar (self-contained mode)
        cosmic_threshold: AF threshold for filtering COSMIC (self-contained mode)
        cross_dataset_threshold: AF threshold for filtering all mutations (cross-dataset mode)

    Returns:
        Tuple of (filtered DataFrame, statistics dict)
    """
    stats = {
        "total": len(df),
        "gnomad_total": 0,
        "gnomad_retained": 0,
        "clinvar_total": 0,
        "clinvar_removed": 0,
        "clinvar_retained": 0,
        "cosmic_total": 0,
        "cosmic_removed": 0,
        "cosmic_retained": 0,
        "other_total": 0,
        "other_removed": 0,
        "other_retained": 0,
    }

    if df.empty:
        return df, stats

    df = df.copy()
    keep_mask = pd.Series([True] * len(df), index=df.index)

    # Helper to check gnomAD match
    def get_gnomad_af(row):
        key = (
            str(row.get("chromosome", "")),
            str(row.get("position", "")),
            str(row.get("reference", "")),
            str(row.get("alternate", "")),
        )
        return gnomad_lookup.get(key, None)

    # Cross-dataset mode: filter all mutations
    if cross_dataset_threshold is not None:
        for idx, row in df.iterrows():
            source = str(row.get("source", "")).lower()
            gnomad_af = get_gnomad_af(row)

            # Track source counts
            if source == "gnomad":
                stats["gnomad_total"] += 1
            elif source == "clinvar":
                stats["clinvar_total"] += 1
            elif source == "cosmic":
                stats["cosmic_total"] += 1
            else:
                stats["other_total"] += 1

            # Check if should be removed
            if gnomad_af is not None and gnomad_af >= cross_dataset_threshold:
                keep_mask[idx] = False
                if source == "gnomad":
                    pass  # gnomad_retained calculated at end
                elif source == "clinvar":
                    stats["clinvar_removed"] += 1
                elif source == "cosmic":
                    stats["cosmic_removed"] += 1
                else:
                    stats["other_removed"] += 1

    # Self-contained mode: filter ClinVar/COSMIC against gnomAD
    else:
        for idx, row in df.iterrows():
            source = str(row.get("source", "")).lower()
            gnomad_af = get_gnomad_af(row)

            if source == "gnomad":
                stats["gnomad_total"] += 1
                # ALWAYS keep gnomAD variants
                keep_mask[idx] = True

            elif source == "clinvar":
                stats["clinvar_total"] += 1
                if clinvar_threshold is not None:
                    if gnomad_af is not None and gnomad_af >= clinvar_threshold:
                        keep_mask[idx] = False
                        stats["clinvar_removed"] += 1
                # If no clinvar_threshold, keep all ClinVar

            elif source == "cosmic":
                stats["cosmic_total"] += 1
                if cosmic_threshold is not None:
                    if gnomad_af is not None and gnomad_af >= cosmic_threshold:
                        keep_mask[idx] = False
                        stats["cosmic_removed"] += 1
                # If no cosmic_threshold, keep all COSMIC

            else:
                stats["other_total"] += 1
                # Keep unknown sources

    # Calculate retained counts
    filtered_df = df[keep_mask]

    if "source" in filtered_df.columns:
        source_counts = filtered_df["source"].str.lower().value_counts()
        stats["gnomad_retained"] = source_counts.get("gnomad", 0)
        stats["clinvar_retained"] = source_counts.get("clinvar", 0)
        stats["cosmic_retained"] = source_counts.get("cosmic", 0)
        stats["other_retained"] = (
            len(filtered_df)
            - stats["gnomad_retained"]
            - stats["clinvar_retained"]
            - stats["cosmic_retained"]
        )

    return filtered_df, stats


def recalculate_counts(df: pd.DataFrame) -> Dict:
    """Recalculate count columns from filtered mutation DataFrame.

    Args:
        df: DataFrame with filtered mutations for a feature

    Returns:
        Dict with recalculated count columns only (not metadata)
    """
    result = {}

    # Determine impact column
    impact_col = "impact_validated" if "impact_validated" in df.columns else "impact"

    # Detect available sources from data
    if df.empty or "source" not in df.columns:
        detected_sources = []
    else:
        detected_sources = [s.lower() for s in df["source"].dropna().unique()]

    # Total mutations
    result["total_mutations"] = len(df)

    # Source databases
    if df.empty or "source" not in df.columns:
        result["source_databases"] = ""
    else:
        mutation_sources = set(s for s in df["source"].dropna().unique())
        result["source_databases"] = (
            ",".join(sorted(mutation_sources)) if mutation_sources else ""
        )

    # Count by impact type
    for impact in IMPACT_TYPES:
        if df.empty:
            result[f"count_{impact}"] = 0
        else:
            impact_df = df[df[impact_col] == impact.replace("_", " ")]
            result[f"count_{impact}"] = len(impact_df)

    # Start site counts
    if df.empty or "in_alt_start_site" not in df.columns:
        result["count_in_alt_start_site"] = 0
        result["alt_start_variant_ids"] = ""
    else:
        alt_start_df = df[df["in_alt_start_site"] == True]
        result["count_in_alt_start_site"] = len(alt_start_df)
        if "variant_id" in alt_start_df.columns:
            ids = alt_start_df["variant_id"].dropna().unique().tolist()
            result["alt_start_variant_ids"] = ",".join(
                str(id) for id in ids if str(id).strip()
            )
        else:
            result["alt_start_variant_ids"] = ""

    if df.empty or "in_canonical_start_site" not in df.columns:
        result["count_in_canonical_start_site"] = 0
        result["canonical_start_variant_ids"] = ""
    else:
        can_start_df = df[df["in_canonical_start_site"] == True]
        result["count_in_canonical_start_site"] = len(can_start_df)
        if "variant_id" in can_start_df.columns:
            ids = can_start_df["variant_id"].dropna().unique().tolist()
            result["canonical_start_variant_ids"] = ",".join(
                str(id) for id in ids if str(id).strip()
            )
        else:
            result["canonical_start_variant_ids"] = ""

    if df.empty or "is_alt_start_loss" not in df.columns:
        result["alternative_start_loss_count"] = 0
        result["alternative_start_loss_variant_ids"] = ""
    else:
        alt_loss_df = df[df["is_alt_start_loss"] == True]
        result["alternative_start_loss_count"] = len(alt_loss_df)
        if "variant_id" in alt_loss_df.columns:
            ids = alt_loss_df["variant_id"].dropna().unique().tolist()
            result["alternative_start_loss_variant_ids"] = ",".join(
                str(id) for id in ids if str(id).strip()
            )
        else:
            result["alternative_start_loss_variant_ids"] = ""

    if df.empty or "is_canonical_start_loss" not in df.columns:
        result["canonical_start_loss_count"] = 0
        result["canonical_start_loss_variant_ids"] = ""
    else:
        can_loss_df = df[df["is_canonical_start_loss"] == True]
        result["canonical_start_loss_count"] = len(can_loss_df)
        if "variant_id" in can_loss_df.columns:
            ids = can_loss_df["variant_id"].dropna().unique().tolist()
            result["canonical_start_loss_variant_ids"] = ",".join(
                str(id) for id in ids if str(id).strip()
            )
        else:
            result["canonical_start_loss_variant_ids"] = ""

    # Per-source counts and aggregates
    gnomad_allele_count_sum = 0
    gnomad_allele_count_hom_sum = 0
    cosmic_sample_count_sum = 0
    gnomad_allele_count_by_impact = {
        f"gnomad_allele_count_{imp}": 0 for imp in IMPACT_TYPES
    }
    cosmic_sample_count_by_impact = {
        f"cosmic_sample_count_{imp}": 0 for imp in IMPACT_TYPES
    }

    for source in SOURCES:
        if df.empty or "source" not in df.columns:
            source_df = pd.DataFrame()
        else:
            source_df = df[df["source"].str.lower() == source]
        result[f"count_{source}"] = len(source_df)

        if source == "gnomad" and len(source_df) > 0:
            if "gnomad_allele_count" in source_df.columns:
                gnomad_allele_count_sum = int(
                    source_df["gnomad_allele_count"].fillna(0).sum()
                )
            if "gnomad_allele_count_hom" in source_df.columns:
                gnomad_allele_count_hom_sum = int(
                    source_df["gnomad_allele_count_hom"].fillna(0).sum()
                )

        if source == "cosmic" and len(source_df) > 0:
            if "cosmic_sample_count" in source_df.columns:
                cosmic_sample_count_sum = int(
                    source_df["cosmic_sample_count"].fillna(0).sum()
                )

        # Source × Impact matrix
        for impact in IMPACT_TYPES:
            impact_normalized = impact.replace(" ", "_").lower()
            if len(source_df) > 0:
                source_impact_df = source_df[
                    source_df[impact_col] == impact.replace("_", " ")
                ]
            else:
                source_impact_df = pd.DataFrame()
            result[f"count_{source}_{impact_normalized}"] = len(source_impact_df)

            # Collect IDs
            if len(source_impact_df) > 0 and "variant_id" in source_impact_df.columns:
                ids = source_impact_df["variant_id"].dropna().unique().tolist()
                result[f"ids_{source}_{impact_normalized}"] = ",".join(
                    str(id) for id in ids if str(id).strip()
                )
            else:
                result[f"ids_{source}_{impact_normalized}"] = ""

            # Per-impact allele/sample counts
            if source == "gnomad" and len(source_impact_df) > 0:
                if "gnomad_allele_count" in source_impact_df.columns:
                    gnomad_allele_count_by_impact[
                        f"gnomad_allele_count_{impact_normalized}"
                    ] = int(source_impact_df["gnomad_allele_count"].fillna(0).sum())

            if source == "cosmic" and len(source_impact_df) > 0:
                if "cosmic_sample_count" in source_impact_df.columns:
                    cosmic_sample_count_by_impact[
                        f"cosmic_sample_count_{impact_normalized}"
                    ] = int(source_impact_df["cosmic_sample_count"].fillna(0).sum())

    # Add gnomAD aggregate columns
    result["gnomad_allele_count"] = gnomad_allele_count_sum
    result["gnomad_allele_count_hom"] = gnomad_allele_count_hom_sum
    result.update(gnomad_allele_count_by_impact)

    # Add COSMIC aggregate columns
    result["cosmic_sample_count"] = cosmic_sample_count_sum
    result.update(cosmic_sample_count_by_impact)

    # All variant IDs
    if df.empty or "variant_id" not in df.columns:
        result["all_variant_ids"] = ""
    else:
        all_ids = df["variant_id"].dropna().unique().tolist()
        result["all_variant_ids"] = ",".join(
            str(id) for id in all_ids if str(id).strip()
        )

    return result


def get_mutation_csv_path(
    dataset_dir: Path, gene_name: str, transcript_id: str, feature_id: str
) -> Optional[Path]:
    """Get path to mutation CSV file for a feature.

    Args:
        dataset_dir: Path to dataset directory
        gene_name: Gene name
        transcript_id: Transcript ID (with version)
        feature_id: Feature ID (format: feature_type_codon_start_end)

    Returns:
        Path to mutation CSV file, or None if not found
    """
    # Parse feature_id to get components
    # Format: extension_CTG_53321463_53321480
    parts = feature_id.split("_")
    if len(parts) < 4:
        return None

    feature_type = parts[0]
    codon = parts[1]
    start = parts[2]
    end = parts[3]

    # Construct expected filename
    # Format: ENST00000209873.8_extension_CTG_53321463_53321480_mutations.csv
    filename = f"{transcript_id}_{feature_type}_{codon}_{start}_{end}_mutations.csv"

    # Look for the file
    mutations_dir = dataset_dir / "mutations" / gene_name / transcript_id
    csv_path = mutations_dir / filename

    if csv_path.exists():
        return csv_path

    # Try without version in transcript path
    transcript_base = transcript_id.rsplit(".", 1)[0]
    mutations_dir_alt = dataset_dir / "mutations" / gene_name / transcript_base
    csv_path_alt = mutations_dir_alt / filename

    if csv_path_alt.exists():
        return csv_path_alt

    return None


def get_column_order() -> List[str]:
    """Get the standard column order for isoform_level_results."""
    columns = []

    # Section 1: Feature identification
    columns.extend(
        [
            "gene_name",
            "transcript_id",
            "feature_id",
            "bed_name",
            "feature_type",
            "feature_start",
            "feature_end",
            "feature_length_aa",
            "aa_difference_from_canonical",
            "total_mutations",
            "source_databases",
        ]
    )

    # Section 2: Summary counts by impact
    columns.extend(
        [
            "count_missense_variant",
            "count_nonsense_variant",
            "count_frameshift_variant",
            "count_inframe_deletion",
            "count_inframe_insertion",
            "count_synonymous_variant",
            "count_in_alt_start_site",
            "count_in_canonical_start_site",
            "alternative_start_loss_count",
            "canonical_start_loss_count",
        ]
    )

    # Section 3: Per-source totals and per-impact allele/sample counts
    columns.extend(
        [
            "count_clinvar",
            "count_gnomad",
            "gnomad_allele_count",
            "gnomad_allele_count_hom",
            "gnomad_allele_count_missense_variant",
            "gnomad_allele_count_nonsense_variant",
            "gnomad_allele_count_frameshift_variant",
            "gnomad_allele_count_inframe_deletion",
            "gnomad_allele_count_inframe_insertion",
            "gnomad_allele_count_synonymous_variant",
            "count_cosmic",
            "cosmic_sample_count",
            "cosmic_sample_count_missense_variant",
            "cosmic_sample_count_nonsense_variant",
            "cosmic_sample_count_frameshift_variant",
            "cosmic_sample_count_inframe_deletion",
            "cosmic_sample_count_inframe_insertion",
            "cosmic_sample_count_synonymous_variant",
        ]
    )

    # Section 4: Source × Impact matrix (count + IDs pairs)
    for source in SOURCES:
        for impact in IMPACT_TYPES:
            columns.append(f"count_{source}_{impact}")
            columns.append(f"ids_{source}_{impact}")

    # Section 5: Start codons and variant IDs
    columns.extend(
        [
            "alternative_start_codon",
            "canonical_start_codon",
            "alt_start_variant_ids",
            "alternative_start_loss_variant_ids",
            "canonical_start_variant_ids",
            "canonical_start_loss_variant_ids",
            "all_variant_ids",
        ]
    )

    return columns


def process_self_contained(
    dataset_dir: Path,
    mane_file: Path,
    clinvar_threshold: Optional[float] = None,
    cosmic_threshold: Optional[float] = None,
    dry_run: bool = False,
) -> Dict:
    """Process dataset in self-contained mode.

    Filters ClinVar and/or COSMIC against gnomAD variants in the same dataset.
    gnomAD variants are always retained. Uses existing MANE file as base.

    Returns:
        Statistics dict
    """
    # Load MANE file as base
    logger.info(f"Loading MANE file: {mane_file}")
    mane_df = pd.read_csv(mane_file)
    logger.info(f"  Loaded {len(mane_df)} features")

    # Build gnomAD lookup from same dataset
    gnomad_lookup = build_gnomad_lookup(dataset_dir)

    total_stats = {
        "gnomad_lookup_size": len(gnomad_lookup),
        "gnomad_total": 0,
        "gnomad_retained": 0,
        "clinvar_total": 0,
        "clinvar_removed": 0,
        "clinvar_retained": 0,
        "cosmic_total": 0,
        "cosmic_removed": 0,
        "cosmic_retained": 0,
        "other_total": 0,
        "other_retained": 0,
        "features_processed": 0,
        "features_not_found": 0,
    }

    filter_suffix = get_filter_suffix(
        clinvar_threshold=clinvar_threshold, cosmic_threshold=cosmic_threshold
    )

    # Process each feature in the MANE file
    for idx, row in mane_df.iterrows():
        gene_name = row["gene_name"]
        transcript_id = row["transcript_id"]
        feature_id = row["feature_id"]

        # Find corresponding mutation CSV
        csv_path = get_mutation_csv_path(
            dataset_dir, gene_name, transcript_id, feature_id
        )
        if csv_path is None:
            total_stats["features_not_found"] += 1
            continue

        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            logger.warning(f"Error reading {csv_path}: {e}")
            total_stats["features_not_found"] += 1
            continue

        # Apply filtering
        filtered_df, stats = filter_by_gnomad_overlap(
            df,
            gnomad_lookup,
            clinvar_threshold=clinvar_threshold,
            cosmic_threshold=cosmic_threshold,
        )

        # Accumulate stats
        for key in [
            "gnomad_total",
            "gnomad_retained",
            "clinvar_total",
            "clinvar_removed",
            "clinvar_retained",
            "cosmic_total",
            "cosmic_removed",
            "cosmic_retained",
            "other_total",
            "other_retained",
        ]:
            total_stats[key] += stats.get(key, 0)
        total_stats["features_processed"] += 1

        if not dry_run:
            # Write filtered per-gene CSV
            filtered_csv_path = csv_path.parent / csv_path.name.replace(
                "_mutations.csv", f"_mutations_filtered_{filter_suffix}.csv"
            )
            filtered_df.to_csv(filtered_csv_path, index=False)

        # Recalculate counts and update mane_df
        new_counts = recalculate_counts(filtered_df)
        for col, val in new_counts.items():
            if col in mane_df.columns:
                mane_df.at[idx, col] = val

    if dry_run:
        return total_stats

    if total_stats["features_not_found"] > 0:
        logger.warning(
            f"Could not find mutation CSVs for {total_stats['features_not_found']} features"
        )

    # Add traceability columns
    if clinvar_threshold is not None:
        mane_df["clinvar_gnomad_af_threshold"] = clinvar_threshold
    if cosmic_threshold is not None:
        mane_df["cosmic_gnomad_af_threshold"] = cosmic_threshold

    # Write filtered output
    output_path = (
        dataset_dir
        / "mutations"
        / f"isoform_level_results_filtered_{filter_suffix}.csv"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mane_df.to_csv(output_path, index=False)
    logger.info(f"Filtered results written to: {output_path}")

    return total_stats


def process_cross_dataset(
    dataset_dir: Path,
    mane_file: Path,
    gnomad_source_dir: Path,
    threshold: float,
    dry_run: bool = False,
) -> Dict:
    """Process dataset in cross-dataset mode.

    Filters test mutations against gnomAD variants from another dataset.
    Removes mutations that match gnomAD variants with AF >= threshold.
    Uses existing MANE file as base.

    Returns:
        Statistics dict
    """
    # Load MANE file as base
    logger.info(f"Loading MANE file: {mane_file}")
    mane_df = pd.read_csv(mane_file)
    logger.info(f"  Loaded {len(mane_df)} features")

    # Build gnomAD lookup from source dataset
    gnomad_lookup = build_gnomad_lookup(gnomad_source_dir)

    total_stats = {
        "gnomad_lookup_size": len(gnomad_lookup),
        "gnomad_source": str(gnomad_source_dir),
        "total_mutations": 0,
        "removed": 0,
        "retained": 0,
        "clinvar_total": 0,
        "clinvar_removed": 0,
        "clinvar_retained": 0,
        "cosmic_total": 0,
        "cosmic_removed": 0,
        "cosmic_retained": 0,
        "other_total": 0,
        "other_removed": 0,
        "other_retained": 0,
        "features_processed": 0,
        "features_not_found": 0,
    }

    filter_suffix = get_filter_suffix(threshold=threshold)

    # Process each feature in the MANE file
    for idx, row in mane_df.iterrows():
        gene_name = row["gene_name"]
        transcript_id = row["transcript_id"]
        feature_id = row["feature_id"]

        # Find corresponding mutation CSV
        csv_path = get_mutation_csv_path(
            dataset_dir, gene_name, transcript_id, feature_id
        )
        if csv_path is None:
            total_stats["features_not_found"] += 1
            continue

        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            logger.warning(f"Error reading {csv_path}: {e}")
            total_stats["features_not_found"] += 1
            continue

        original_count = len(df)

        # Apply filtering
        filtered_df, stats = filter_by_gnomad_overlap(
            df,
            gnomad_lookup,
            cross_dataset_threshold=threshold,
        )

        # Accumulate stats
        total_stats["total_mutations"] += original_count
        total_stats["retained"] += len(filtered_df)
        total_stats["removed"] += original_count - len(filtered_df)

        for key in [
            "clinvar_total",
            "clinvar_removed",
            "clinvar_retained",
            "cosmic_total",
            "cosmic_removed",
            "cosmic_retained",
            "other_total",
            "other_removed",
            "other_retained",
        ]:
            total_stats[key] += stats.get(key, 0)
        total_stats["features_processed"] += 1

        if not dry_run:
            # Write filtered per-gene CSV
            filtered_csv_path = csv_path.parent / csv_path.name.replace(
                "_mutations.csv", f"_mutations_filtered_{filter_suffix}.csv"
            )
            filtered_df.to_csv(filtered_csv_path, index=False)

        # Recalculate counts and update mane_df
        new_counts = recalculate_counts(filtered_df)
        for col, val in new_counts.items():
            if col in mane_df.columns:
                mane_df.at[idx, col] = val

    if dry_run:
        return total_stats

    if total_stats["features_not_found"] > 0:
        logger.warning(
            f"Could not find mutation CSVs for {total_stats['features_not_found']} features"
        )

    # Add traceability columns
    mane_df["gnomad_af_threshold"] = threshold
    mane_df["gnomad_source"] = str(gnomad_source_dir)

    # Write filtered output
    output_path = (
        dataset_dir
        / "mutations"
        / f"isoform_level_results_filtered_{filter_suffix}.csv"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mane_df.to_csv(output_path, index=False)
    logger.info(f"Filtered results written to: {output_path}")

    return total_stats


def print_summary(
    stats: Dict,
    dataset_dir: Path,
    clinvar_threshold: Optional[float] = None,
    cosmic_threshold: Optional[float] = None,
    cross_dataset_threshold: Optional[float] = None,
):
    """Print filtering summary statistics."""
    print("\n" + "=" * 60)
    print("gnomAD Overlap Filtering Summary")
    print("=" * 60)
    print(f"Dataset: {dataset_dir}")

    if cross_dataset_threshold is not None:
        print(
            f"Mode: Cross-dataset (gnomAD source: {stats.get('gnomad_source', 'N/A')})"
        )
        print(f"Threshold: AF >= {cross_dataset_threshold} (remove common variants)")
    else:
        print("Mode: Self-contained")
        if clinvar_threshold is not None:
            print(f"ClinVar threshold: AF >= {clinvar_threshold}")
        if cosmic_threshold is not None:
            print(f"COSMIC threshold: AF >= {cosmic_threshold}")

    print(f"\nFeatures processed: {stats.get('features_processed', 0)}")
    print(f"gnomAD lookup size: {stats.get('gnomad_lookup_size', 0):,}")
    print()

    if cross_dataset_threshold is not None:
        # Cross-dataset summary
        total = stats.get("total_mutations", 0)
        removed = stats.get("removed", 0)
        retained = stats.get("retained", 0)
        pct_removed = (removed / total * 100) if total > 0 else 0

        print("Overall:")
        print(f"  Total mutations: {total:,}")
        print(f"  Removed (matched common gnomAD): {removed:,} ({pct_removed:.1f}%)")
        print(f"  Retained: {retained:,}")
    else:
        # Self-contained summary - gnomAD is always retained
        print("gnomAD variants (reference - always retained):")
        print(f"  Total: {stats.get('gnomad_total', 0):,}")
        print(f"  Retained: {stats.get('gnomad_retained', 0):,}")

    # ClinVar stats
    if clinvar_threshold is not None or cross_dataset_threshold is not None:
        clinvar_total = stats.get("clinvar_total", 0)
        clinvar_removed = stats.get("clinvar_removed", 0)
        clinvar_retained = stats.get("clinvar_retained", 0)
        pct_removed = (
            (clinvar_removed / clinvar_total * 100) if clinvar_total > 0 else 0
        )

        print(f"\nClinVar variants:")
        print(f"  Total: {clinvar_total:,}")
        print(f"  Removed: {clinvar_removed:,} ({pct_removed:.1f}%)")
        print(f"  Retained: {clinvar_retained:,}")

    # COSMIC stats
    if cosmic_threshold is not None or cross_dataset_threshold is not None:
        cosmic_total = stats.get("cosmic_total", 0)
        cosmic_removed = stats.get("cosmic_removed", 0)
        cosmic_retained = stats.get("cosmic_retained", 0)
        pct_removed = (cosmic_removed / cosmic_total * 100) if cosmic_total > 0 else 0

        print(f"\nCOSMIC variants:")
        print(f"  Total: {cosmic_total:,}")
        print(f"  Removed: {cosmic_removed:,} ({pct_removed:.1f}%)")
        print(f"  Retained: {cosmic_retained:,}")

    # Other sources
    other_total = stats.get("other_total", 0)
    if other_total > 0:
        other_removed = stats.get("other_removed", 0)
        other_retained = stats.get("other_retained", 0)
        pct_removed = (other_removed / other_total * 100) if other_total > 0 else 0

        print(f"\nOther sources:")
        print(f"  Total: {other_total:,}")
        print(f"  Removed: {other_removed:,} ({pct_removed:.1f}%)")
        print(f"  Retained: {other_retained:,}")

    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Filter mutations by removing those matching common gnomAD variants",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Self-contained mode - filter ClinVar and COSMIC (using mane file)
  python scripts/filter_gnomad_af.py results/hela \\
      --mane-file results/hela/mutations/isoform_level_results_mane.csv \\
      --clinvar-threshold 0.01 --cosmic-threshold 0.001

  # Cross-dataset mode (filter test dataset against gnomAD from reference)
  python scripts/filter_gnomad_af.py results/hela_bch \\
      --mane-file results/hela_bch/mutations/isoform_level_results_mane.csv \\
      --gnomad-source results/hela --threshold 0.01

  # Dry run (statistics only, no output files)
  python scripts/filter_gnomad_af.py results/hela \\
      --mane-file results/hela/mutations/isoform_level_results_mane.csv \\
      --clinvar-threshold 0.01 --dry-run
        """,
    )

    parser.add_argument(
        "dataset_dir",
        type=Path,
        help="Path to dataset directory to process",
    )

    # MANE file argument (required)
    parser.add_argument(
        "--mane-file",
        type=Path,
        required=True,
        help="Path to isoform_level_results_mane.csv file to use as base",
    )

    # Self-contained mode arguments
    parser.add_argument(
        "--clinvar-threshold",
        type=float,
        default=None,
        help="AF threshold for filtering ClinVar variants (self-contained mode)",
    )
    parser.add_argument(
        "--cosmic-threshold",
        type=float,
        default=None,
        help="AF threshold for filtering COSMIC variants (self-contained mode)",
    )

    # Cross-dataset mode arguments
    parser.add_argument(
        "--gnomad-source",
        type=Path,
        default=None,
        help="Path to dataset with gnomAD data (enables cross-dataset mode)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="AF threshold for cross-dataset filtering",
    )

    # Common arguments
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report statistics without writing output files",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose output",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate arguments
    if not args.dataset_dir.exists():
        logger.error(f"Dataset directory not found: {args.dataset_dir}")
        sys.exit(1)

    if not args.mane_file.exists():
        logger.error(f"MANE file not found: {args.mane_file}")
        sys.exit(1)

    # Determine mode and validate
    cross_dataset_mode = args.gnomad_source is not None
    self_contained_mode = (
        args.clinvar_threshold is not None or args.cosmic_threshold is not None
    )

    if cross_dataset_mode and self_contained_mode:
        logger.error(
            "Cannot use --gnomad-source with --clinvar-threshold or --cosmic-threshold"
        )
        logger.error(
            "Choose either cross-dataset mode (--gnomad-source + --threshold) or self-contained mode (--clinvar-threshold and/or --cosmic-threshold)"
        )
        sys.exit(1)

    if cross_dataset_mode:
        if args.threshold is None:
            logger.error("Cross-dataset mode requires --threshold")
            sys.exit(1)
        if not args.gnomad_source.exists():
            logger.error(f"gnomAD source directory not found: {args.gnomad_source}")
            sys.exit(1)
    elif not self_contained_mode:
        logger.error(
            "Must specify either --clinvar-threshold and/or --cosmic-threshold (self-contained mode)"
        )
        logger.error("or --gnomad-source with --threshold (cross-dataset mode)")
        sys.exit(1)

    # Process
    if cross_dataset_mode:
        logger.info("Running in cross-dataset mode")
        stats = process_cross_dataset(
            dataset_dir=args.dataset_dir,
            mane_file=args.mane_file,
            gnomad_source_dir=args.gnomad_source,
            threshold=args.threshold,
            dry_run=args.dry_run,
        )
        print_summary(
            stats,
            args.dataset_dir,
            cross_dataset_threshold=args.threshold,
        )
    else:
        logger.info("Running in self-contained mode")
        stats = process_self_contained(
            dataset_dir=args.dataset_dir,
            mane_file=args.mane_file,
            clinvar_threshold=args.clinvar_threshold,
            cosmic_threshold=args.cosmic_threshold,
            dry_run=args.dry_run,
        )
        print_summary(
            stats,
            args.dataset_dir,
            clinvar_threshold=args.clinvar_threshold,
            cosmic_threshold=args.cosmic_threshold,
        )

    if args.dry_run:
        print("\n[DRY RUN - No output files written]")


if __name__ == "__main__":
    main()
