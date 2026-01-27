#!/usr/bin/env python3
"""Analyze differential mutation burden between datasets using MANE-filtered data."""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse


def load_mane_data(dataset: str, results_dir: str = "../results") -> pd.DataFrame:
    """Load MANE-filtered isoform level results."""
    path = Path(results_dir) / dataset / "mutations" / "isoform_level_results_mane.csv"
    df = pd.read_csv(path)
    df["dataset"] = dataset
    return df


def main():
    """Main entry point for differential mutation analysis."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=["hela_bch", "hela_msk"],
        help="Datasets to compare",
    )
    parser.add_argument("--results-dir", default="../results")
    parser.add_argument("--top-n", type=int, default=30)
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    output_lines = []

    def out(line=""):
        output_lines.append(line)

    out(f"{'=' * 90}")
    out(f"DIFFERENTIAL MUTATION ANALYSIS: {' vs '.join(args.datasets).upper()}")
    out(f"Using MANE-filtered variants only")
    out(f"{'=' * 90}")

    # Load all datasets
    dfs = {}
    for dataset in args.datasets:
        dfs[dataset] = load_mane_data(dataset, args.results_dir)
        out(f"\n{dataset}: {len(dfs[dataset])} features loaded")

    # =========================================================================
    # PART 1: Overall mutation counts by dataset
    # =========================================================================
    out(f"\n{'─' * 90}")
    out("PART 1: OVERALL MUTATION COUNTS")
    out(f"{'─' * 90}")

    for name, df in dfs.items():
        total_missense = df["count_custom_missense_variant"].sum()
        total_nonsense = df["count_custom_nonsense_variant"].sum()
        total_frameshift = df["count_custom_frameshift_variant"].sum()
        total_synonymous = df["count_custom_synonymous_variant"].sum()
        features_with_mutations = (df["total_mutations"] > 0).sum()

        out(f"\n{name}:")
        out(f"  Features with mutations: {features_with_mutations}")
        out(f"  Total missense variants: {int(total_missense)}")
        out(f"  Total nonsense variants: {int(total_nonsense)}")
        out(f"  Total frameshift variants: {int(total_frameshift)}")
        out(f"  Total synonymous variants: {int(total_synonymous)}")

    # =========================================================================
    # PART 2: Gene-level comparison
    # =========================================================================
    out(f"\n{'─' * 90}")
    out("PART 2: GENE-LEVEL MUTATION BURDEN COMPARISON")
    out(f"{'─' * 90}")

    # Aggregate by gene
    gene_summaries = {}
    for name, df in dfs.items():
        gene_df = (
            df.groupby("gene_name")
            .agg(
                {
                    "count_custom_missense_variant": "sum",
                    "count_custom_nonsense_variant": "sum",
                    "count_custom_frameshift_variant": "sum",
                    "count_custom_synonymous_variant": "sum",
                    "total_mutations": "sum",
                }
            )
            .reset_index()
        )
        gene_df.columns = [
            "gene",
            f"{name}_missense",
            f"{name}_nonsense",
            f"{name}_frameshift",
            f"{name}_synonymous",
            f"{name}_total",
        ]
        gene_summaries[name] = gene_df

    # Merge gene summaries
    merged = gene_summaries[args.datasets[0]]
    for name in args.datasets[1:]:
        merged = merged.merge(gene_summaries[name], on="gene", how="outer")
    merged = merged.fillna(0)

    # Calculate differential
    d1, d2 = args.datasets[0], args.datasets[1]
    merged["diff_missense"] = merged[f"{d1}_missense"] - merged[f"{d2}_missense"]
    merged["diff_total"] = merged[f"{d1}_total"] - merged[f"{d2}_total"]

    # Genes with more mutations in dataset 1
    out(f"\nGenes with MORE missense mutations in {d1} vs {d2} (top {args.top_n}):")
    out(f"{'Gene':<15} {d1:>12} {d2:>12} {'Diff':>10}")
    out(f"{'-' * 15} {'-' * 12} {'-' * 12} {'-' * 10}")

    more_in_d1 = merged[merged["diff_missense"] > 0].nlargest(
        args.top_n, "diff_missense"
    )
    for _, row in more_in_d1.iterrows():
        out(
            f"{row['gene']:<15} {int(row[f'{d1}_missense']):>12} {int(row[f'{d2}_missense']):>12} {int(row['diff_missense']):>+10}"
        )

    # Genes with more mutations in dataset 2
    out(f"\nGenes with MORE missense mutations in {d2} vs {d1} (top {args.top_n}):")
    out(f"{'Gene':<15} {d1:>12} {d2:>12} {'Diff':>10}")
    out(f"{'-' * 15} {'-' * 12} {'-' * 12} {'-' * 10}")

    more_in_d2 = merged[merged["diff_missense"] < 0].nsmallest(
        args.top_n, "diff_missense"
    )
    for _, row in more_in_d2.iterrows():
        out(
            f"{row['gene']:<15} {int(row[f'{d1}_missense']):>12} {int(row[f'{d2}_missense']):>12} {int(row['diff_missense']):>+10}"
        )

    # =========================================================================
    # PART 3: Unique mutations per dataset
    # =========================================================================
    out(f"\n{'─' * 90}")
    out("PART 3: DATASET-SPECIFIC MUTATIONS")
    out(f"{'─' * 90}")

    # Genes with mutations only in one dataset
    only_d1 = merged[(merged[f"{d1}_missense"] > 0) & (merged[f"{d2}_missense"] == 0)]
    only_d2 = merged[(merged[f"{d2}_missense"] > 0) & (merged[f"{d1}_missense"] == 0)]

    out(f"\nGenes with missense mutations ONLY in {d1}: {len(only_d1)}")
    if len(only_d1) > 0:
        top_only_d1 = only_d1.nlargest(min(args.top_n, len(only_d1)), f"{d1}_missense")
        out(f"Top {min(args.top_n, len(only_d1))} by count:")
        for _, row in top_only_d1.iterrows():
            out(
                f"  {row['gene']:<15} {int(row[f'{d1}_missense']):>5} missense variants"
            )

    out(f"\nGenes with missense mutations ONLY in {d2}: {len(only_d2)}")
    if len(only_d2) > 0:
        top_only_d2 = only_d2.nlargest(min(args.top_n, len(only_d2)), f"{d2}_missense")
        out(f"Top {min(args.top_n, len(only_d2))} by count:")
        for _, row in top_only_d2.iterrows():
            out(
                f"  {row['gene']:<15} {int(row[f'{d2}_missense']):>5} missense variants"
            )

    # =========================================================================
    # PART 4: Feature-level details for top differential genes
    # =========================================================================
    out(f"\n{'─' * 90}")
    out("PART 4: FEATURE-LEVEL DETAILS FOR TOP DIFFERENTIAL GENES")
    out(f"{'─' * 90}")

    # Get top 5 genes with biggest differences
    top_diff_genes = merged.nlargest(5, "diff_missense")["gene"].tolist()

    for gene in top_diff_genes:
        out(f"\n{gene}:")
        for name, df in dfs.items():
            gene_features = df[df["gene_name"] == gene]
            if len(gene_features) > 0:
                out(f"  {name}:")
                for _, feat in gene_features.iterrows():
                    if feat["count_custom_missense_variant"] > 0:
                        out(
                            f"    {feat['feature_type']}: {int(feat['count_custom_missense_variant'])} missense"
                        )
                        # Show variant IDs if available
                        if (
                            pd.notna(feat.get("ids_custom_missense_variant"))
                            and feat["ids_custom_missense_variant"]
                        ):
                            variants = str(feat["ids_custom_missense_variant"]).split(
                                ","
                            )[:5]
                            out(
                                f"      Variants: {', '.join(variants)}{'...' if len(variants) == 5 else ''}"
                            )

    # =========================================================================
    # OUTPUT
    # =========================================================================
    output_text = "\n".join(output_lines)

    if args.output:
        with open(args.output, "w") as f:
            f.write(output_text)
        print(f"Analysis written to {args.output}")
    else:
        print(output_text)

    # Save merged comparison as TSV
    output_dir = Path(args.results_dir) / "comparisons"
    output_dir.mkdir(exist_ok=True)
    merged.to_csv(
        output_dir / f"gene_mutation_comparison_{'_vs_'.join(args.datasets)}.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
