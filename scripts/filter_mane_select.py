#!/usr/bin/env python3
"""
Filter isoform_level_results.csv to keep maximum one truncation + one extension per gene,
prioritizing MANE Select transcripts.

Usage:
    python filter_mane_select.py [--test GENE_NAME]

Options:
    --test GENE_NAME    Test mode: show filtering for specific gene only
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
            if line.startswith('#'):
                continue

            # Only look at transcript lines with MANE_Select tag
            if '\ttranscript\t' in line and 'MANE_Select' in line:
                # Extract gene_name and transcript_id from attributes
                attrs = line.split('\t')[8]

                gene_name = None
                transcript_id = None

                for attr in attrs.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_name "'):
                        gene_name = attr.split('"')[1]
                    elif attr.startswith('transcript_id "'):
                        transcript_id = attr.split('"')[1]
                        # Strip version number (e.g., ENST00000368474.9 → ENST00000368474)
                        transcript_id = transcript_id.rsplit('.', 1)[0]

                if gene_name and transcript_id:
                    mane_transcripts[gene_name] = transcript_id

    print(f"  Found {len(mane_transcripts)} genes with MANE Select transcripts")
    return mane_transcripts

def filter_isoforms(isoform_csv, mane_transcripts, test_gene=None):
    """Filter isoform data to keep max 1 truncation + 1 extension per gene."""
    print(f"\nLoading isoform data from {isoform_csv}...")
    df = pd.read_csv(isoform_csv)
    print(f"  Original: {len(df)} isoforms from {df['gene_name'].nunique()} genes")

    # Strip version from transcript_id in dataframe
    df['transcript_base'] = df['transcript_id'].str.rsplit('.', n=1).str[0]

    # Add MANE Select flag
    df['is_mane_select'] = df.apply(
        lambda row: mane_transcripts.get(row['gene_name']) == row['transcript_base'],
        axis=1
    )

    mane_count = df['is_mane_select'].sum()
    print(f"  {mane_count} isoforms match MANE Select transcripts")

    # Test mode: show one gene's filtering
    if test_gene:
        print(f"\n=== Test Mode: Filtering for gene '{test_gene}' ===")
        gene_df = df[df['gene_name'] == test_gene].copy()

        if len(gene_df) == 0:
            print(f"  ERROR: Gene '{test_gene}' not found in data")
            return

        print(f"\nBefore filtering ({len(gene_df)} rows):")
        print(gene_df[['gene_name', 'transcript_id', 'feature_type', 'total_mutations', 'is_mane_select']])

        # Apply filtering logic
        filtered_gene = filter_by_feature_type(gene_df)

        print(f"\nAfter filtering ({len(filtered_gene)} rows):")
        print(filtered_gene[['gene_name', 'transcript_id', 'feature_type', 'total_mutations', 'is_mane_select']])

        return

    # Full dataset filtering
    print(f"\nFiltering: max 1 truncation + 1 extension per gene...")
    filtered_df = df.groupby('gene_name', group_keys=False).apply(filter_by_feature_type)

    print(f"  Filtered: {len(filtered_df)} isoforms")
    print(f"    Extensions: {(filtered_df['feature_type'] == 'extension').sum()}")
    print(f"    Truncations: {(filtered_df['feature_type'] == 'truncation').sum()}")
    print(f"    MANE Select: {filtered_df['is_mane_select'].sum()}")
    print(f"  Eliminated: {len(df) - len(filtered_df)} rows ({100*(len(df) - len(filtered_df))/len(df):.1f}%)")

    # Drop helper columns
    filtered_df = filtered_df.drop(columns=['transcript_base', 'is_mane_select'])

    return filtered_df

def filter_by_feature_type(gene_df):
    """For one gene, select max 1 truncation + 1 extension."""
    result_rows = []

    for feature_type in ['extension', 'truncation']:
        feature_df = gene_df[gene_df['feature_type'] == feature_type]

        if len(feature_df) == 0:
            continue
        elif len(feature_df) == 1:
            result_rows.append(feature_df)
        else:
            # Multiple rows: prioritize MANE Select, then highest total_mutations, then lexicographic transcript_id
            sorted_df = feature_df.sort_values(
                by=['is_mane_select', 'total_mutations', 'transcript_id'],
                ascending=[False, False, True]
            )
            result_rows.append(sorted_df.head(1))

    if result_rows:
        return pd.concat(result_rows, ignore_index=True)
    else:
        return pd.DataFrame()

def main():
    parser = argparse.ArgumentParser(description='Filter isoforms by MANE Select')
    parser.add_argument('--test', type=str, help='Test mode: filter only this gene')
    args = parser.parse_args()

    # Paths
    gtf_file = Path('/lab/barcheese01/mdiberna/swissisoform/data/genome_data/gencode.v47.annotation.gtf')
    isoform_csv = Path('/lab/barcheese01/mdiberna/swissisoform/results/hela/mutations/isoform_level_results.csv')
    output_csv = Path('/lab/barcheese01/mdiberna/swissisoform/results/hela/mutations/isoform_level_results_minimal.csv')

    # Parse MANE Select transcripts from GTF
    mane_transcripts = parse_gtf_for_mane(gtf_file)

    # Filter isoforms
    filtered_df = filter_isoforms(isoform_csv, mane_transcripts, test_gene=args.test)

    # Save output (only in full mode)
    if not args.test and filtered_df is not None:
        print(f"\nSaving to {output_csv}...")
        filtered_df.to_csv(output_csv, index=False)
        print("✅ Done!")

if __name__ == '__main__':
    main()
