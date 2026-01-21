#!/usr/bin/env python3
"""Summary analysis of SwissIsoform pipeline results.

This script analyzes mutation analysis results and localization predictions
to provide comprehensive summaries, identify genes with interesting
localization changes, and generate ranked CSV outputs for analysis.

Generates the following outputs:
  - Text summaries (mutation_summary.txt, localization_summary.txt)
  - Basic CSVs (genes_with_localization_changes.csv, gene_level_summary.csv)
  - Ranked CSVs (summary_localization_changes.csv, summary_deleterious_burden.csv)
  - Detailed analysis CSVs (differential_localization_pairs_{mode}.csv,
    differential_localization_mutations_{mode}.csv)

Usage:
    python summarize_results.py                  # Uses DATASET env var or defaults to hela
    python summarize_results.py hela_bch         # Analyze specific dataset
    DATASET=hela_msk python summarize_results.py # Via environment variable
"""

import os
import sys
import warnings
from pathlib import Path
import pandas as pd
import re

# Suppress pandas warnings
warnings.filterwarnings("ignore", category=FutureWarning)

from swissisoform.summary import SummaryAnalyzer

# Localization compartments
COMPARTMENTS = [
    "Cytoplasm",
    "Nucleus",
    "Extracellular",
    "Cell membrane",
    "Mitochondrion",
    "Plastid",
    "Endoplasmic reticulum",
    "Lysosome/Vacuole",
    "Golgi apparatus",
    "Peroxisome",
]


# ============================================================================
# Ranked Summary Generation (from summarize_dataset.py)
# ============================================================================


def get_primary_loc(row):
    """Get the compartment with highest probability."""
    probs = {c: row[c] for c in COMPARTMENTS if c in row.index}
    return max(probs, key=probs.get) if probs else None


def get_max_shift(row1, row2):
    """Calculate maximum localization probability shift between two predictions."""
    max_shift = 0
    for c in COMPARTMENTS:
        if c in row1.index and c in row2.index:
            shift = abs(row2[c] - row1[c])
            if shift > max_shift:
                max_shift = shift
    return max_shift


def calculate_localization_shift(row1, row2):
    """Calculate detailed localization probability shifts between two sequences.

    Returns dict with shift statistics including total shift, max shift value,
    and which compartment had the maximum shift.
    """
    shifts = {}
    total_shift = 0
    for c in COMPARTMENTS:
        if c in row1.index and c in row2.index:
            shift = row2[c] - row1[c]
            shifts[c] = shift
            total_shift += abs(shift)

    if not shifts:
        return {
            "shifts": {},
            "total_shift": 0,
            "max_shift_comp": None,
            "max_shift_value": 0,
        }

    max_shift_comp = max(shifts, key=lambda x: abs(shifts[x]))
    max_shift_value = max(abs(v) for v in shifts.values())

    return {
        "shifts": shifts,
        "total_shift": total_shift,
        "max_shift_comp": max_shift_comp,
        "max_shift_value": max_shift_value,
    }


def parse_protein_id(protein_id):
    """Parse Protein_ID to extract gene, transcript, feature_type, and mutation info."""
    parts = protein_id.split("_")
    gene = parts[0]
    transcript = parts[1] if len(parts) > 1 else None

    is_mutant = "_mutated_" in protein_id
    is_canonical = protein_id.endswith("_canonical")

    if is_canonical:
        feature_type = "canonical"
    elif "_extension_" in protein_id:
        feature_type = "extension"
    elif "_truncation_" in protein_id:
        feature_type = "truncation"
    else:
        feature_type = "unknown"

    mutation_info = None
    if is_mutant:
        match = re.search(r"_mutated_(\d+)_([^_]+)_([A-Z]\d+[A-Z])", protein_id)
        if match:
            mutation_info = {
                "position": match.group(1),
                "change": match.group(2),
                "aa_change": match.group(3),
            }

    return {
        "gene": gene,
        "transcript": transcript,
        "feature_type": feature_type,
        "is_mutant": is_mutant,
        "is_canonical": is_canonical,
        "mutation_info": mutation_info,
    }


def detect_source_prefixes(df):
    """Detect all available mutation source prefixes from column names.

    Returns list of source prefixes (e.g., ['custom'], ['clinvar', 'cosmic']).
    """
    if df is None or df.empty:
        return []

    cols = df.columns.tolist()
    found = []

    # Check for common source prefixes
    for prefix in ["custom", "clinvar", "cosmic", "gnomad"]:
        if f"count_{prefix}_missense_variant" in cols:
            found.append(prefix)

    return found


def find_mutation_file(mut_dir):
    """Find the appropriate filtered mutation file.

    Returns (file_path, expected_sources) where expected_sources is a list of
    source prefixes that should be used based on the filename.
    """
    # Try different file patterns in order of preference
    # Format: (filename, [expected_sources])
    patterns = [
        ("isoform_level_results_mane.csv", None),  # MANE-annotated, auto-detect sources
        ("isoform_level_results_filtered_gnomad5e-05.csv", ["custom"]),
        (
            "isoform_level_results_filtered_clinvar5e-05_cosmic5e-05.csv",
            ["clinvar", "cosmic"],
        ),
        ("isoform_level_results_filtered_*.csv", None),  # Auto-detect
        ("isoform_level_results.csv", None),  # Unfiltered, auto-detect
    ]

    for pattern, expected_sources in patterns:
        if "*" in pattern:
            # Glob pattern
            matches = list(mut_dir.glob(pattern))
            if matches:
                # Use first match, auto-detect source
                return matches[0], None
        else:
            file_path = mut_dir / pattern
            if file_path.exists():
                return file_path, expected_sources

    return None, None


def generate_localization_changes(pairs_df, output_path):
    """Generate summary_localization_changes.csv."""
    pairs_df = pairs_df.copy()

    parsed = pairs_df["Protein_ID"].apply(parse_protein_id)
    pairs_df["gene"] = parsed.apply(lambda x: x["gene"])
    pairs_df["transcript"] = parsed.apply(lambda x: x["transcript"])
    pairs_df["feature_type"] = parsed.apply(lambda x: x["feature_type"])
    pairs_df["is_canonical"] = parsed.apply(lambda x: x["is_canonical"])
    pairs_df["is_mutant"] = parsed.apply(lambda x: x["is_mutant"])

    results = []
    for (gene, transcript), group in pairs_df.groupby(["gene", "transcript"]):
        canonical = group[group["is_canonical"]]
        alternatives = group[(~group["is_canonical"]) & (~group["is_mutant"])]

        if len(canonical) == 0 or len(alternatives) == 0:
            continue

        can_row = canonical.iloc[0]
        can_primary = get_primary_loc(can_row)

        for _, alt_row in alternatives.iterrows():
            alt_primary = get_primary_loc(alt_row)
            shift = get_max_shift(can_row, alt_row)

            results.append(
                {
                    "gene": gene,
                    "transcript": transcript,
                    "feature_type": alt_row["feature_type"],
                    "canonical_loc": can_primary,
                    "alternative_loc": alt_primary,
                    "shift": round(shift, 4),
                    "loc_changed": can_primary != alt_primary,
                }
            )

    if results:
        df = pd.DataFrame(results)
        df = df.sort_values("shift", ascending=False)
        df.to_csv(output_path, index=False)
        return df
    return pd.DataFrame()


def identify_loc_altering_mutations(mutations_df):
    """Identify which mutations change localization."""
    if mutations_df is None or mutations_df.empty:
        return {}

    mutations_df = mutations_df.copy()

    parsed = mutations_df["Protein_ID"].apply(parse_protein_id)
    mutations_df["gene"] = parsed.apply(lambda x: x["gene"])
    mutations_df["transcript"] = parsed.apply(lambda x: x["transcript"])
    mutations_df["feature_type"] = parsed.apply(lambda x: x["feature_type"])
    mutations_df["is_mutant"] = parsed.apply(lambda x: x["is_mutant"])
    mutations_df["mutation_info"] = parsed.apply(lambda x: x["mutation_info"])

    loc_altering = {}

    for (gene, transcript, feature_type), group in mutations_df.groupby(
        ["gene", "transcript", "feature_type"]
    ):
        if feature_type in ("canonical", "unknown"):
            continue

        refs = group[~group["is_mutant"]]
        mutants = group[group["is_mutant"]]

        if len(refs) == 0 or len(mutants) == 0:
            continue

        ref_row = refs.iloc[0]
        ref_primary = get_primary_loc(ref_row)

        altering_positions = set()
        for _, mut_row in mutants.iterrows():
            mut_primary = get_primary_loc(mut_row)
            if mut_primary != ref_primary and mut_row["mutation_info"]:
                altering_positions.add(mut_row["mutation_info"]["position"])

        if altering_positions:
            key = (gene, transcript, feature_type)
            loc_altering[key] = altering_positions

    return loc_altering


def extract_position_from_variant_id(var_id):
    """Extract genomic position from variant ID."""
    parts = var_id.split("_")
    if len(parts) >= 4:
        pos_candidate = parts[2]
        if pos_candidate.isdigit() and len(pos_candidate) >= 6:
            return pos_candidate
    return None


def count_loc_altering_missense(feature_row, loc_altering_map, source_prefix):
    """Count how many missense variants in a feature are localization-altering."""
    gene = feature_row["gene_name"]
    transcript = feature_row["transcript_id"]
    feature_type = feature_row["feature_type"]

    key = (gene, transcript, feature_type)
    altering_positions = loc_altering_map.get(key, set())

    if not altering_positions:
        return 0

    # Get missense IDs using the appropriate source prefix
    ids_col = f"ids_{source_prefix}_missense_variant"
    missense_ids = feature_row.get(ids_col, "")
    if pd.isna(missense_ids) or missense_ids == "":
        return 0

    count = 0
    for var_id in re.split(r"[;,\n]", str(missense_ids)):
        var_id = var_id.strip()
        if not var_id:
            continue
        pos = extract_position_from_variant_id(var_id)
        if pos and pos in altering_positions:
            count += 1

    return count


def safe_int(val):
    """Convert value to int, handling NaN and None."""
    if pd.isna(val) or val is None:
        return 0
    return int(val)


def generate_deleterious_burden(mane_df, loc_altering_map, output_path, source_prefix):
    """Generate summary_deleterious_burden.csv."""
    if mane_df is None or mane_df.empty:
        return pd.DataFrame()

    results = []

    for _, row in mane_df.iterrows():
        # Get counts using the appropriate source prefix
        count_nonsense = safe_int(row.get(f"count_{source_prefix}_nonsense_variant", 0))
        count_frameshift = safe_int(
            row.get(f"count_{source_prefix}_frameshift_variant", 0)
        )
        count_inframe_del = safe_int(
            row.get(f"count_{source_prefix}_inframe_deletion", 0)
        )
        count_inframe_ins = safe_int(
            row.get(f"count_{source_prefix}_inframe_insertion", 0)
        )
        count_missense = safe_int(row.get(f"count_{source_prefix}_missense_variant", 0))
        count_synonymous = safe_int(
            row.get(f"count_{source_prefix}_synonymous_variant", 0)
        )

        count_loc_altering = count_loc_altering_missense(
            row, loc_altering_map, source_prefix
        )
        count_other_missense = max(0, count_missense - count_loc_altering)

        total_deleterious = (
            count_nonsense
            + count_frameshift
            + count_inframe_del
            + count_inframe_ins
            + count_loc_altering
        )

        if count_synonymous > 0:
            del_to_syn_ratio = round(total_deleterious / count_synonymous, 3)
        else:
            del_to_syn_ratio = float("inf") if total_deleterious > 0 else 0

        mane_status = "Yes" if row.get("MANE", "") == "Yes" else "No"

        results.append(
            {
                "gene": row["gene_name"],
                "transcript": row["transcript_id"],
                "MANE": mane_status,
                "feature_id": row.get("feature_id", ""),
                "feature_type": row.get("feature_type", ""),
                "count_nonsense": count_nonsense,
                "count_frameshift": count_frameshift,
                "count_inframe_del": count_inframe_del,
                "count_inframe_ins": count_inframe_ins,
                "count_loc_altering_missense": count_loc_altering,
                "count_other_missense": count_other_missense,
                "total_deleterious": total_deleterious,
                "count_synonymous": count_synonymous,
                "deleterious_to_synonymous_ratio": del_to_syn_ratio,
            }
        )

    if results:
        df = pd.DataFrame(results)
        df = df[(df["total_deleterious"] > 0) | (df["count_synonymous"] > 0)]
        df = df.sort_values("total_deleterious", ascending=False)
        df.to_csv(output_path, index=False)
        return df
    return pd.DataFrame()


# ============================================================================
# Detailed Localization Analysis (from analyze_localization.py)
# ============================================================================


def analyze_pairs_detailed(pairs_df):
    """Analyze canonical vs alternative isoform localization differences in detail.

    Returns DataFrame with detailed shift statistics for each gene-isoform pair.
    """
    results = []

    # Parse all protein IDs
    parsed = pairs_df["Protein_ID"].apply(parse_protein_id)
    pairs_df = pairs_df.copy()
    pairs_df["gene"] = parsed.apply(lambda x: x["gene"])
    pairs_df["feature_type"] = parsed.apply(lambda x: x["feature_type"])
    pairs_df["is_canonical"] = parsed.apply(lambda x: x["is_canonical"])
    pairs_df["is_mutant"] = parsed.apply(lambda x: x["is_mutant"])

    # Group by gene
    for gene, group in pairs_df.groupby("gene"):
        canonical = group[group["is_canonical"]]
        extensions = group[
            (group["feature_type"] == "extension") & (~group["is_mutant"])
        ]

        if len(canonical) == 0 or len(extensions) == 0:
            continue

        canonical_row = canonical.iloc[0]
        canonical_loc = canonical_row.get("Localizations", "")
        canonical_primary = get_primary_loc(canonical_row)

        for _, ext_row in extensions.iterrows():
            ext_loc = ext_row.get("Localizations", "")
            ext_primary = get_primary_loc(ext_row)

            shift_info = calculate_localization_shift(canonical_row, ext_row)

            # Determine if localization changed
            loc_changed = canonical_loc != ext_loc
            primary_changed = canonical_primary != ext_primary

            results.append(
                {
                    "gene": gene,
                    "canonical_id": canonical_row["Protein_ID"],
                    "extension_id": ext_row["Protein_ID"],
                    "canonical_localization": canonical_loc,
                    "extension_localization": ext_loc,
                    "canonical_primary": canonical_primary,
                    "extension_primary": ext_primary,
                    "localization_changed": loc_changed,
                    "primary_changed": primary_changed,
                    "total_shift": round(shift_info["total_shift"], 4),
                    "max_shift_compartment": shift_info["max_shift_comp"],
                    "max_shift_value": round(shift_info["max_shift_value"], 4),
                }
            )

    return pd.DataFrame(results)


def analyze_mutations_detailed(mutations_df, pairs_df=None):
    """Analyze mutation effects on localization in detail.

    Returns DataFrame with detailed statistics for each mutation.
    """
    if mutations_df is None or mutations_df.empty:
        return pd.DataFrame()

    results = []

    # Parse all protein IDs
    parsed = mutations_df["Protein_ID"].apply(parse_protein_id)
    mutations_df = mutations_df.copy()
    mutations_df["gene"] = parsed.apply(lambda x: x["gene"])
    mutations_df["feature_type"] = parsed.apply(lambda x: x["feature_type"])
    mutations_df["is_mutant"] = parsed.apply(lambda x: x["is_mutant"])
    mutations_df["mutation_info"] = parsed.apply(lambda x: x["mutation_info"])

    # Group by gene and feature_type
    for (gene, feature_type), group in mutations_df.groupby(["gene", "feature_type"]):
        if feature_type == "canonical":
            continue

        # Find reference (non-mutant) and mutant sequences
        references = group[~group["is_mutant"]]
        mutants = group[group["is_mutant"]]

        if len(references) == 0 or len(mutants) == 0:
            continue

        ref_row = references.iloc[0]
        ref_loc = ref_row.get("Localizations", "")
        ref_primary = get_primary_loc(ref_row)

        for _, mut_row in mutants.iterrows():
            mut_loc = mut_row.get("Localizations", "")
            mut_primary = get_primary_loc(mut_row)

            shift_info = calculate_localization_shift(ref_row, mut_row)

            # Determine if localization changed
            loc_changed = ref_loc != mut_loc
            primary_changed = ref_primary != mut_primary

            # Extract mutation string
            mutation_str = ""
            if mut_row["mutation_info"]:
                mutation_str = (
                    f"{mut_row['mutation_info']['position']}_"
                    f"{mut_row['mutation_info']['change']}_"
                    f"{mut_row['mutation_info']['aa_change']}"
                )

            results.append(
                {
                    "gene": gene,
                    "feature_type": feature_type,
                    "reference_id": ref_row["Protein_ID"],
                    "mutant_id": mut_row["Protein_ID"],
                    "mutation": mutation_str,
                    "reference_localization": ref_loc,
                    "mutant_localization": mut_loc,
                    "reference_primary": ref_primary,
                    "mutant_primary": mut_primary,
                    "localization_changed": loc_changed,
                    "primary_changed": primary_changed,
                    "total_shift": round(shift_info["total_shift"], 4),
                    "max_shift_compartment": shift_info["max_shift_comp"],
                    "max_shift_value": round(shift_info["max_shift_value"], 4),
                }
            )

    return pd.DataFrame(results)


def process_ranked_summaries(dataset, results_dir, verbose=False):
    """Generate ranked summary CSVs for localization changes and deleterious burden."""
    mut_dir = results_dir / "mutations"
    loc_dir = results_dir / "localization"

    print("\n" + "=" * 80)
    print("GENERATING RANKED SUMMARIES")
    print("=" * 80)

    # Find and load mutation data
    mane_file, expected_sources = find_mutation_file(mut_dir)
    mane_df = None
    sources_to_process = []

    if mane_file and mane_file.exists():
        mane_df = pd.read_csv(mane_file)
        print(f"\nLoaded mutation data: {len(mane_df)} features from {mane_file.name}")

        # Detect available sources in the data
        available_sources = detect_source_prefixes(mane_df)
        if available_sources:
            print(f"Available sources in data: {', '.join(available_sources)}")

        # Determine which sources to process (use all available)
        if expected_sources:
            sources_to_process = expected_sources
        elif available_sources:
            sources_to_process = available_sources

        if not sources_to_process:
            print("Warning: Could not detect mutation source prefix")
    else:
        print(f"\nWarning: No mutation file found in {mut_dir}")

    # Process both Accurate and Fast modes
    modes = ["Accurate", "Fast"]
    all_results = []

    for source_prefix in sources_to_process if sources_to_process else [None]:
        if source_prefix and len(sources_to_process) > 1:
            print(f"\n{'-' * 80}")
            print(f"Processing source: {source_prefix.upper()}")
            print(f"{'-' * 80}")

        for mode in modes:
            summary_dir = results_dir / "summary" / mode.lower()
            summary_dir.mkdir(parents=True, exist_ok=True)

            pairs_file = loc_dir / f"protein_sequences_pairs_{mode}_results.csv"
            mutations_loc_file = loc_dir / f"protein_sequences_mutations_{mode}_results.csv"

            if not pairs_file.exists():
                if verbose:
                    print(f"  Skipping {mode} mode - no localization data")
                continue

            print(f"\n  Processing {mode} mode...")

            result = {
                "mode": mode,
                "source": source_prefix,
                "loc_changes_count": 0,
                "loc_changed_count": 0,
                "del_burden_count": 0,
                "total_deleterious": 0,
                "total_synonymous": 0,
            }

            # Generate localization changes
            loc_changes_out = summary_dir / "summary_localization_changes.csv"
            pairs_df = pd.read_csv(pairs_file)
            loc_changes_df = generate_localization_changes(pairs_df, loc_changes_out)

            if len(loc_changes_df) > 0:
                changed = loc_changes_df[loc_changes_df["loc_changed"]]
                result["loc_changes_count"] = len(loc_changes_df)
                result["loc_changed_count"] = len(changed)
                print(
                    f"    Localization changes: {len(loc_changes_df)} features, {len(changed)} with change"
                )

            # Generate detailed pairs analysis
            pairs_detailed_out = summary_dir / f"differential_localization_pairs_{mode}.csv"
            pairs_detailed_df = analyze_pairs_detailed(pairs_df)
            if len(pairs_detailed_df) > 0:
                pairs_detailed_df.to_csv(pairs_detailed_out, index=False)
                primary_changed = pairs_detailed_df[pairs_detailed_df["primary_changed"]]
                print(
                    f"    Detailed pairs analysis: {len(pairs_detailed_df)} pairs, "
                    f"{len(primary_changed)} with primary change"
                )

            # Identify localization-altering mutations
            loc_altering_map = {}
            mutations_loc_df = None
            if mutations_loc_file.exists():
                mutations_loc_df = pd.read_csv(mutations_loc_file)
                loc_altering_map = identify_loc_altering_mutations(mutations_loc_df)
                print(f"    Loc-altering mutations: {len(loc_altering_map)} features")

                # Generate detailed mutations analysis
                mutations_detailed_out = summary_dir / f"differential_localization_mutations_{mode}.csv"
                mutations_detailed_df = analyze_mutations_detailed(mutations_loc_df, pairs_df)
                if len(mutations_detailed_df) > 0:
                    mutations_detailed_df.to_csv(mutations_detailed_out, index=False)
                    primary_changed_mut = mutations_detailed_df[mutations_detailed_df["primary_changed"]]
                    print(
                        f"    Detailed mutation analysis: {len(mutations_detailed_df)} mutations, "
                        f"{len(primary_changed_mut)} with primary change"
                    )

            # Generate deleterious burden
            if mane_df is not None and source_prefix is not None:
                del_burden_out = summary_dir / "summary_deleterious_burden.csv"
                del_burden_df = generate_deleterious_burden(
                    mane_df, loc_altering_map, del_burden_out, source_prefix
                )

                if len(del_burden_df) > 0:
                    result["del_burden_count"] = len(del_burden_df)
                    result["total_deleterious"] = int(del_burden_df["total_deleterious"].sum())
                    result["total_synonymous"] = int(del_burden_df["count_synonymous"].sum())
                    print(
                        f"    Deleterious burden: {len(del_burden_df)} features, "
                        f"{result['total_deleterious']} deleterious variants"
                    )

            print(f"    Saved to: {summary_dir}/")
            all_results.append(result)

    # Print summary
    if all_results:
        print(f"\n{'=' * 80}")
        print("RANKED SUMMARY RESULTS")
        print(f"{'=' * 80}")

        for result in all_results:
            if result["loc_changes_count"] > 0 or result["del_burden_count"] > 0:
                source_label = f" [{result['source']}]" if result["source"] else ""
                print(f"\n{result['mode']} mode{source_label}:")
                print(
                    f"  Localization changes: {result['loc_changes_count']} features "
                    f"({result['loc_changed_count']} changed)"
                )
                print(f"  Deleterious burden: {result['del_burden_count']} features")
                print(
                    f"  Total deleterious: {result['total_deleterious']}, "
                    f"Total synonymous: {result['total_synonymous']}"
                )


# ============================================================================
# Main Function
# ============================================================================


def main():
    """Main analysis function that processes the specified dataset."""
    # Get dataset from command line arg, environment variable, or default
    if len(sys.argv) > 1:
        dataset = sys.argv[1]
    else:
        dataset = os.environ.get("DATASET", "hela")

    print("SwissIsoform Pipeline Results Summary")
    print("=" * 70)
    print(f"Dataset: {dataset}")

    results_dir = Path("../results") / dataset

    # Initialize the analyzer
    analyzer = SummaryAnalyzer()

    print(f"\n{'=' * 20} ANALYZING {dataset.upper()} DATASET {'=' * 20}")

    # Check if this dataset has any data
    if not analyzer.dataset_has_data(dataset):
        print(f"Error: {dataset} dataset has no data available")
        sys.exit(1)

    # Run standard summary analysis (text summaries, basic CSVs)
    try:
        analyzer.analyze_dataset(dataset)
        print(f"\nStandard summary analysis completed successfully for {dataset}!")
    except Exception as e:
        print(f"Error in standard summary analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # Generate ranked summaries (CSV files for further analysis)
    try:
        process_ranked_summaries(dataset, results_dir, verbose=False)
        print(f"\nRanked summary generation completed successfully!")
    except Exception as e:
        print(f"Error generating ranked summaries: {e}")
        import traceback
        traceback.print_exc()
        # Don't exit - standard summaries were successful

    print(f"\n{'=' * 70}")
    print(f"All results saved to: ../results/{dataset}/summary/")


if __name__ == "__main__":
    main()
