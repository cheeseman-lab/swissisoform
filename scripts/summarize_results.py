#!/usr/bin/env python3
"""Summary analysis of SwissIsoform pipeline results.

This script analyzes mutation analysis results and localization predictions
to provide comprehensive summaries and identify genes with interesting
localization changes for a specific dataset.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import warnings
import argparse

# Suppress pandas warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def load_mutation_results(dataset):
    """Load mutation analysis results for a dataset."""
    base_path = Path(f"../results/{dataset}/mutations")

    gene_results = None
    pair_results = None

    gene_file = base_path / "gene_level_results.csv"
    pair_file = base_path / "truncation_level_results.csv"

    if gene_file.exists():
        gene_results = pd.read_csv(gene_file)
        print(f"Loaded {len(gene_results)} gene-level mutation results for {dataset}")

    if pair_file.exists():
        pair_results = pd.read_csv(pair_file)
        print(
            f"Loaded {len(pair_results)} transcript-truncation pair results for {dataset}"
        )

    return gene_results, pair_results


def load_localization_results(dataset):
    """Load localization prediction results for a dataset."""
    base_path = Path(f"../results/{dataset}/localization")

    results = {}

    # Define the files to look for
    file_patterns = [
        ("pairs_accurate", "protein_sequences_pairs_Accurate_results.csv"),
        ("pairs_fast", "protein_sequences_pairs_Fast_results.csv"),
        ("mutations_accurate", "protein_sequences_mutations_Accurate_results.csv"),
        ("mutations_fast", "protein_sequences_mutations_Fast_results.csv"),
    ]

    for key, filename in file_patterns:
        file_path = base_path / filename
        if file_path.exists():
            try:
                df = pd.read_csv(file_path)

                # Standardize DeepLoc column names to our expected format
                column_mapping = {
                    "Protein_ID": "Sequence_ID",
                    "Localizations": "Localisation",
                }

                df = df.rename(columns=column_mapping)

                results[key] = df
                print(f"Loaded {len(df)} {key} localization predictions for {dataset}")

            except Exception as e:
                print(f"Error loading {filename}: {e}")
                continue
        else:
            print(f"Missing {key} localization predictions for {dataset}")

    return results


def analyze_mutations(dataset, gene_results, pair_results):
    """Analyze mutation results for a specific dataset."""
    print(f"\n=== ANALYZING MUTATIONS FOR {dataset.upper()} DATASET ===")

    summary_lines = []
    summary_lines.append(f"MUTATION ANALYSIS SUMMARY - {dataset.upper()} DATASET")
    summary_lines.append("=" * 60)

    if gene_results is None:
        summary_lines.append("No gene-level mutation results available")
        return summary_lines

    # Basic statistics
    total_genes = len(gene_results)
    successful_genes = len(gene_results[gene_results["status"] == "success"])
    summary_lines.append(f"Total genes processed: {total_genes}")
    summary_lines.append(f"Genes with successful analysis: {successful_genes}")

    if pair_results is not None and not pair_results.empty:
        # Look for mutation columns in pair results
        mutation_columns = [
            col for col in pair_results.columns if col.startswith("mutations_")
        ]

        if mutation_columns:
            summary_lines.append(f"\nMutation breakdown by type:")

            # Count genes with each mutation type
            gene_mutation_counts = {}

            for col in mutation_columns:
                # Extract mutation type from column name
                mutation_type = col.replace("mutations_", "").replace("_", " ").title()

                # Count genes that have at least one mutation of this type
                genes_with_mutation = pair_results[pair_results[col] > 0][
                    "gene_name"
                ].nunique()
                total_mutations = pair_results[col].sum()

                gene_mutation_counts[mutation_type] = genes_with_mutation
                summary_lines.append(
                    f"  {mutation_type}: {genes_with_mutation} genes ({total_mutations} total mutations)"
                )

            # Special combined category: frameshift OR nonsense
            frameshift_col = "mutations_frameshift_variant"
            nonsense_col = "mutations_nonsense_variant"

            if (
                frameshift_col in pair_results.columns
                and nonsense_col in pair_results.columns
            ):
                genes_with_truncating = pair_results[
                    (pair_results[frameshift_col] > 0)
                    | (pair_results[nonsense_col] > 0)
                ]["gene_name"].nunique()
                summary_lines.append(
                    f"  Truncating (Frameshift OR Nonsense): {genes_with_truncating} genes"
                )

            # Total mutations
            total_mutations_all = (
                pair_results["mutation_count_total"].sum()
                if "mutation_count_total" in pair_results.columns
                else 0
            )
            summary_lines.append(
                f"\nTotal mutations across all transcript-truncation pairs: {total_mutations_all}"
            )
        else:
            summary_lines.append("No mutation data found in pair results")
    else:
        summary_lines.append("No pair-level mutation data available")

    return summary_lines


def parse_sequence_id(seq_id):
    """Parse a sequence ID to extract gene, transcript, and variant information."""
    if pd.isna(seq_id) or not seq_id:
        return None, None, None

    seq_id = str(seq_id).strip()
    parts = seq_id.split("_")

    if len(parts) >= 3:
        gene = parts[0]
        transcript = parts[1]
        variant = "_".join(parts[2:])  # Rest is variant info
        return gene, transcript, variant
    elif len(parts) == 2:
        # Sometimes might just be gene_variant
        gene = parts[0]
        transcript = "unknown"
        variant = parts[1]
        return gene, transcript, variant
    else:
        # Can't parse
        return None, None, None


def parse_deeploc_localization(localization):
    """Parse DeepLoc localization which can be multiple locations separated by |."""
    if pd.isna(localization) or not localization:
        return "Unknown"

    localization = str(localization).strip()

    # If multiple localizations, take the first one as primary
    if "|" in localization:
        primary_loc = localization.split("|")[0].strip()
        return primary_loc

    return localization


def analyze_localizations(dataset, loc_results):
    """Analyze localization predictions for a specific dataset."""
    print(f"\n=== ANALYZING LOCALIZATIONS FOR {dataset.upper()} DATASET ===")

    summary_lines = []
    summary_lines.append(f"LOCALIZATION ANALYSIS SUMMARY - {dataset.upper()} DATASET")
    summary_lines.append("=" * 60)

    all_localization_comparisons = []

    if not loc_results:
        summary_lines.append("No localization data available")
        return summary_lines, all_localization_comparisons

    # Prefer accurate results over fast results
    pairs_data = None
    mutations_data = None

    if "pairs_accurate" in loc_results:
        pairs_data = loc_results["pairs_accurate"]
        model_type = "Accurate"
    elif "pairs_fast" in loc_results:
        pairs_data = loc_results["pairs_fast"]
        model_type = "Fast"

    if "mutations_accurate" in loc_results:
        mutations_data = loc_results["mutations_accurate"]
    elif "mutations_fast" in loc_results:
        mutations_data = loc_results["mutations_fast"]

    if pairs_data is None:
        summary_lines.append("No localization data available")
        return summary_lines, all_localization_comparisons

    summary_lines.append(f"Using {model_type} model predictions")
    summary_lines.append(f"Total sequences analyzed: {len(pairs_data)}")

    # Check if we have the required columns
    if "Sequence_ID" not in pairs_data.columns:
        summary_lines.append("Error: No Sequence_ID column found")
        return summary_lines, all_localization_comparisons

    if "Localisation" not in pairs_data.columns:
        summary_lines.append("Error: No Localisation column found")
        return summary_lines, all_localization_comparisons

    # Parse sequence IDs and analyze localizations
    pairs_parsed = []
    for _, row in pairs_data.iterrows():
        # Use the standardized column names
        seq_id = row.get("Sequence_ID", "")
        localization = parse_deeploc_localization(row.get("Localisation", ""))

        if seq_id and localization:
            gene, transcript, variant = parse_sequence_id(seq_id)
            if gene and transcript and variant:
                pairs_parsed.append(
                    {
                        "dataset": dataset,
                        "gene": gene,
                        "transcript": transcript,
                        "variant": variant,
                        "sequence_id": seq_id,
                        "prediction": localization,
                        "confidence": 0,  # DeepLoc doesn't provide a single confidence score
                        "source": "pairs",
                    }
                )

    pairs_df = pd.DataFrame(pairs_parsed)

    if not pairs_df.empty:
        # Group by gene and analyze canonical vs truncated
        gene_localization_changes = []

        for gene in pairs_df["gene"].unique():
            gene_data = pairs_df[pairs_df["gene"] == gene]

            # Find canonical and truncated sequences
            canonical_seqs = gene_data[gene_data["variant"] == "canonical"]
            truncated_seqs = gene_data[gene_data["variant"].str.startswith("trunc")]

            if not canonical_seqs.empty and not truncated_seqs.empty:
                canonical_loc = canonical_seqs.iloc[0]["prediction"]

                for _, trunc_seq in truncated_seqs.iterrows():
                    trunc_loc = trunc_seq["prediction"]

                    if canonical_loc != trunc_loc:
                        gene_localization_changes.append(
                            {
                                "dataset": dataset,
                                "gene": gene,
                                "transcript": trunc_seq["transcript"],
                                "canonical_localization": canonical_loc,
                                "truncated_localization": trunc_loc,
                                "truncated_variant": trunc_seq["variant"],
                                "localization_change_type": "canonical_vs_truncated",
                                "model_type": model_type,
                            }
                        )

        summary_lines.append(
            f"Genes with canonical vs truncated localization changes: {len(gene_localization_changes)}"
        )
        all_localization_comparisons.extend(gene_localization_changes)

        # Analyze mutations if available
        if mutations_data is not None:
            mutations_parsed = []
            for _, row in mutations_data.iterrows():
                # Use the standardized column names
                seq_id = row.get("Sequence_ID", "")
                localization = parse_deeploc_localization(row.get("Localisation", ""))

                if seq_id and localization:
                    gene, transcript, variant = parse_sequence_id(seq_id)
                    if gene and transcript and variant:
                        mutations_parsed.append(
                            {
                                "dataset": dataset,
                                "gene": gene,
                                "transcript": transcript,
                                "variant": variant,
                                "sequence_id": seq_id,
                                "prediction": localization,
                                "confidence": 0,  # DeepLoc doesn't provide a single confidence score
                                "source": "mutations",
                            }
                        )

            mutations_df = pd.DataFrame(mutations_parsed)

            if not mutations_df.empty:
                summary_lines.append(
                    f"Mutation sequences analyzed: {len(mutations_df)}"
                )

                # Find genes with missense mutations causing localization changes
                missense_changes = []

                for gene in mutations_df["gene"].unique():
                    gene_mut_data = mutations_df[mutations_df["gene"] == gene]

                    # Find canonical sequence in pairs data
                    gene_pairs_data = pairs_df[pairs_df["gene"] == gene]
                    canonical_pairs = gene_pairs_data[
                        gene_pairs_data["variant"] == "canonical"
                    ]

                    if not canonical_pairs.empty:
                        canonical_loc = canonical_pairs.iloc[0]["prediction"]

                        # Look for missense mutations with different localizations
                        missense_variants = gene_mut_data[
                            gene_mut_data["variant"].str.contains("mut", na=False)
                        ]

                        for _, mut_seq in missense_variants.iterrows():
                            mut_loc = mut_seq["prediction"]

                            if canonical_loc != mut_loc:
                                missense_changes.append(
                                    {
                                        "dataset": dataset,
                                        "gene": gene,
                                        "transcript": mut_seq["transcript"],
                                        "canonical_localization": canonical_loc,
                                        "mutated_localization": mut_loc,
                                        "mutated_variant": mut_seq["variant"],
                                        "localization_change_type": "canonical_vs_missense",
                                        "model_type": model_type,
                                    }
                                )

                summary_lines.append(
                    f"Genes with missense mutation localization changes: {len(missense_changes)}"
                )
                all_localization_comparisons.extend(missense_changes)

        # Most common localizations
        top_localizations = pairs_df["prediction"].value_counts().head(5)
        summary_lines.append(f"\nTop 5 predicted localizations:")
        for loc, count in top_localizations.items():
            summary_lines.append(f"  {loc}: {count} sequences")

    return summary_lines, all_localization_comparisons


def create_detailed_localization_analysis(dataset, loc_results):
    """Create a detailed analysis of all localization predictions for a specific dataset."""
    print(
        f"\n=== CREATING DETAILED LOCALIZATION ANALYSIS FOR {dataset.upper()} DATASET ==="
    )

    detailed_results = []

    if not loc_results:
        return pd.DataFrame(detailed_results)

    # Process pairs data
    pairs_data = None
    mutations_data = None
    model_type = "Unknown"

    if "pairs_accurate" in loc_results:
        pairs_data = loc_results["pairs_accurate"]
        model_type = "Accurate"
    elif "pairs_fast" in loc_results:
        pairs_data = loc_results["pairs_fast"]
        model_type = "Fast"

    if "mutations_accurate" in loc_results:
        mutations_data = loc_results["mutations_accurate"]
    elif "mutations_fast" in loc_results:
        mutations_data = loc_results["mutations_fast"]

    # Process pairs data
    if pairs_data is not None:
        for _, row in pairs_data.iterrows():
            # Use the standardized column names
            seq_id = row.get("Sequence_ID", "")
            localization = parse_deeploc_localization(row.get("Localisation", ""))

            if seq_id and localization:
                gene, transcript, variant = parse_sequence_id(seq_id)
                if gene and transcript and variant:
                    detailed_results.append(
                        {
                            "dataset": dataset,
                            "gene": gene,
                            "transcript": transcript,
                            "variant": variant,
                            "sequence_id": seq_id,
                            "prediction": localization,
                            "confidence": 0,  # DeepLoc doesn't provide a single confidence score
                            "sequence_type": "canonical"
                            if variant == "canonical"
                            else "truncated",
                            "model_type": model_type,
                            "analysis_type": "pairs",
                        }
                    )

    # Process mutations data
    if mutations_data is not None:
        for _, row in mutations_data.iterrows():
            # Use the standardized column names
            seq_id = row.get("Sequence_ID", "")
            localization = parse_deeploc_localization(row.get("Localisation", ""))

            if seq_id and localization:
                gene, transcript, variant = parse_sequence_id(seq_id)
                if gene and transcript and variant:
                    sequence_type = (
                        "canonical"
                        if variant == "canonical"
                        else ("mutated" if "mut" in variant else "truncated")
                    )
                    detailed_results.append(
                        {
                            "dataset": dataset,
                            "gene": gene,
                            "transcript": transcript,
                            "variant": variant,
                            "sequence_id": seq_id,
                            "prediction": localization,
                            "confidence": 0,  # DeepLoc doesn't provide a single confidence score
                            "sequence_type": sequence_type,
                            "model_type": model_type,
                            "analysis_type": "mutations",
                        }
                    )

    return pd.DataFrame(detailed_results)


def main():
    """Main analysis function."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Analyze SwissIsoform pipeline results for a specific dataset"
    )
    parser.add_argument(
        "--dataset",
        required=True,
        choices=["reduced", "full"],
        help="Dataset to analyze (reduced or full)",
    )

    args = parser.parse_args()
    dataset = args.dataset

    print(f"SwissIsoform Pipeline Results Summary - {dataset.upper()} Dataset")
    print("=" * 70)

    # Create summary directory for this dataset
    summary_dir = Path(f"../results/{dataset}/summary")
    summary_dir.mkdir(parents=True, exist_ok=True)

    # Load mutation results for this dataset
    print(f"\nLoading mutation analysis results for {dataset} dataset...")
    gene_results, pair_results = load_mutation_results(dataset)

    # Load localization results for this dataset
    print(f"\nLoading localization prediction results for {dataset} dataset...")
    loc_results = load_localization_results(dataset)

    # Analyze mutations
    mutation_summary = analyze_mutations(dataset, gene_results, pair_results)

    # Analyze localizations
    localization_summary, localization_comparisons = analyze_localizations(
        dataset, loc_results
    )

    # Create detailed localization analysis
    detailed_localization_df = create_detailed_localization_analysis(
        dataset, loc_results
    )

    # Save summaries
    print(f"\n=== SAVING RESULTS FOR {dataset.upper()} DATASET ===")

    # Save mutation summary
    with open(summary_dir / "mutation_summary.txt", "w") as f:
        f.write("\n".join(mutation_summary))
    print("Saved mutation summary")

    # Save localization summary
    with open(summary_dir / "localization_summary.txt", "w") as f:
        f.write("\n".join(localization_summary))
    print("Saved localization summary")

    # Save genes with localization changes
    if localization_comparisons:
        localization_changes_df = pd.DataFrame(localization_comparisons)
        localization_changes_df.to_csv(
            summary_dir / "genes_with_localization_changes.csv", index=False
        )
        print(f"Saved {len(localization_changes_df)} localization changes")
    else:
        # Create empty file
        pd.DataFrame().to_csv(
            summary_dir / "genes_with_localization_changes.csv", index=False
        )
        print("No localization changes found")

    # Save detailed localization analysis
    if not detailed_localization_df.empty:
        detailed_localization_df.to_csv(
            summary_dir / "detailed_localization_analysis.csv", index=False
        )
        print(
            f"Saved detailed analysis of {len(detailed_localization_df)} localization predictions"
        )
    else:
        # Create empty file
        pd.DataFrame().to_csv(
            summary_dir / "detailed_localization_analysis.csv", index=False
        )
        print("No detailed localization data available")

    print(f"\nSummary analysis for {dataset} dataset completed!")
    print(f"Results saved to: {summary_dir}")


if __name__ == "__main__":
    main()
