#!/usr/bin/env python3
"""Annotate isoform_level_results.csv with MANE Select information.

Adds two columns:
1. MANE (Yes/No) - Whether this transcript is MANE Select
2. transcript_notes - Context about MANE status within gene+feature_type

Usage:
    # Default paths (hela dataset)
    python annotate_mane.py

    # Custom input/output paths
    python annotate_mane.py --input results/mutations/isoform_level_results.csv --output results/mutations/isoform_level_results_mane

    # Test mode for specific gene
    python annotate_mane.py --test BRCA1

Options:
    --input FILE        Input isoform_level_results.csv file
    --output PREFIX     Output file prefix (creates .csv and .xlsx)
    --gtf FILE          GTF file with MANE annotations (default: gencode.v47)
    --test GENE_NAME    Test mode: show annotation for specific gene only
"""

import pandas as pd
import argparse
from pathlib import Path


def parse_gtf_for_mane(gtf_file):
    """Parse GTF file to extract gene → MANE Select transcript mapping."""
    print(f"Parsing GTF file for MANE Select annotations...")
    mane_transcripts = {}  # gene_name → transcript_id (without version)

    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            # Only look at transcript lines with MANE_Select tag
            if "\ttranscript\t" in line and "MANE_Select" in line:
                # Extract gene_name and transcript_id from attributes
                attrs = line.split("\t")[8]

                gene_name = None
                transcript_id = None

                for attr in attrs.split(";"):
                    attr = attr.strip()
                    if attr.startswith('gene_name "'):
                        gene_name = attr.split('"')[1]
                    elif attr.startswith('transcript_id "'):
                        transcript_id = attr.split('"')[1]
                        # Strip version number (e.g., ENST00000368474.9 → ENST00000368474)
                        transcript_id = transcript_id.rsplit(".", 1)[0]

                if gene_name and transcript_id:
                    mane_transcripts[gene_name] = transcript_id

    print(f"  Found {len(mane_transcripts)} genes with MANE Select transcripts")
    return mane_transcripts


def annotate_mane_status(df, mane_transcripts):
    """Add MANE and transcript_notes columns to dataframe."""
    # Strip version from transcript_id
    df["transcript_base"] = df["transcript_id"].str.rsplit(".", n=1).str[0]

    # Add MANE Select flag
    df["MANE"] = df.apply(
        lambda row: "Yes"
        if mane_transcripts.get(row["gene_name"]) == row["transcript_base"]
        else "No",
        axis=1,
    )

    # Generate transcript_notes for each (gene_name, feature_type) group
    def generate_notes(group):
        total_count = len(group)
        mane_count = (group["MANE"] == "Yes").sum()

        notes = []
        for idx, row in group.iterrows():
            if total_count == 1:
                if row["MANE"] == "No":
                    notes.append("Only one, no MANE Select")
                else:
                    notes.append("")  # Only one and it's MANE - no note needed
            elif mane_count == 0:
                notes.append("No MANE Select")
            elif mane_count > 1:
                # Unusual: multiple MANE transcripts for same gene+feature_type
                notes.append("Multiple MANE Select")
            else:
                # mane_count == 1 and total_count > 1
                if row["MANE"] == "Yes":
                    notes.append("")  # This is the MANE one - no note needed
                else:
                    notes.append("Non-MANE (MANE exists)")

        return notes

    # Apply notes generation per (gene_name, feature_type)
    df["transcript_notes"] = ""
    for (gene, feature), group in df.groupby(["gene_name", "feature_type"]):
        notes = generate_notes(group)
        df.loc[group.index, "transcript_notes"] = notes

    # Drop helper column
    df = df.drop(columns=["transcript_base"])

    # Reorder columns: insert MANE and transcript_notes after transcript_id
    cols = df.columns.tolist()
    transcript_id_idx = cols.index("transcript_id")

    # Remove MANE and transcript_notes from wherever they are
    cols = [c for c in cols if c not in ["MANE", "transcript_notes"]]

    # Insert after transcript_id
    cols.insert(transcript_id_idx + 1, "MANE")
    cols.insert(transcript_id_idx + 2, "transcript_notes")

    df = df[cols]

    return df


def main():
    """Annotate isoform results with MANE Select information."""
    parser = argparse.ArgumentParser(
        description="Annotate isoforms with MANE Select info"
    )
    parser.add_argument(
        "--input", type=str, help="Input isoform_level_results.csv file"
    )
    parser.add_argument(
        "--output", type=str, help="Output file prefix (creates .csv and .xlsx)"
    )
    parser.add_argument("--gtf", type=str, help="GTF file with MANE annotations")
    parser.add_argument(
        "--test", type=str, help="Test mode: show annotation for specific gene only"
    )
    args = parser.parse_args()

    # Paths - use arguments or defaults
    gtf_file = (
        Path(args.gtf)
        if args.gtf
        else Path(
            "/lab/barcheese01/mdiberna/swissisoform/data/genome_data/gencode.v47.annotation.gtf"
        )
    )
    isoform_csv = (
        Path(args.input)
        if args.input
        else Path(
            "/lab/barcheese01/mdiberna/swissisoform/results/hela/mutations/isoform_level_results.csv"
        )
    )

    if args.output:
        output_prefix = Path(args.output)
        output_csv = output_prefix.with_suffix(".csv")
        output_xlsx = output_prefix.with_suffix(".xlsx")
    else:
        output_csv = Path(
            "/lab/barcheese01/mdiberna/swissisoform/results/hela/mutations/isoform_level_results_mane.csv"
        )
        output_xlsx = output_csv.with_suffix(".xlsx")

    # Parse MANE Select transcripts from GTF
    mane_transcripts = parse_gtf_for_mane(gtf_file)

    # Load isoform data
    print(f"\nLoading isoform data from {isoform_csv}...")
    df = pd.read_csv(isoform_csv)
    print(f"  Loaded: {len(df)} isoforms from {df['gene_name'].nunique()} genes")

    # Annotate with MANE info
    print(f"\nAnnotating with MANE Select information...")
    annotated_df = annotate_mane_status(df, mane_transcripts)

    mane_count = (annotated_df["MANE"] == "Yes").sum()
    notes_count = (annotated_df["transcript_notes"] != "").sum()
    print(f"  MANE Select isoforms: {mane_count}")
    print(f"  Isoforms with notes: {notes_count}")

    # Test mode: show one gene
    if args.test:
        print(f"\n=== Test Mode: Showing annotations for gene '{args.test}' ===\n")
        gene_df = annotated_df[annotated_df["gene_name"] == args.test]

        if len(gene_df) == 0:
            print(f"  ERROR: Gene '{args.test}' not found in data")
            return

        # Show relevant columns
        display_cols = [
            "gene_name",
            "transcript_id",
            "MANE",
            "transcript_notes",
            "feature_type",
            "total_mutations",
        ]
        print(gene_df[display_cols].to_string(index=False))
        return

    # Save full annotated dataset
    print(f"\nSaving annotated data to {output_csv}...")
    annotated_df.to_csv(output_csv, index=False)

    print(f"Saving Excel file to {output_xlsx}...")
    annotated_df.to_excel(output_xlsx, index=False, engine="openpyxl")

    # Show summary statistics
    print("\n=== Summary Statistics ===")
    print(f"Total isoforms: {len(annotated_df)}")
    print(
        f"MANE Select isoforms: {mane_count} ({100 * mane_count / len(annotated_df):.1f}%)"
    )
    print(f"\nNote categories:")
    note_counts = annotated_df["transcript_notes"].value_counts()
    for note, count in note_counts.items():
        if note == "":
            print(f"  (no note): {count}")
        else:
            print(f"  '{note}': {count}")

    print("\n✅ Done!")


if __name__ == "__main__":
    main()
