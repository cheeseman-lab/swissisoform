#!/usr/bin/env python3
"""
5' UTR Length Change Analysis - Standalone Script

This script analyzes ClinVar variants that affect 5' UTR length in MANE Select transcripts.
It identifies variants where indels/duplications cause the UTR to cross a length threshold
(default 40bp), which may affect translation efficiency.

The script automatically downloads required data files if not present.

Usage:
    python run_analysis.py [--min-size 10] [--threshold 40] [--force-download]

Output:
    results/5utr_variants_full.csv - All 5'UTR variants with calculated lengths
    results/5utr_variants_filtered.csv - Variants crossing threshold
    results/5utr_length_changes.png - Slope plot visualization
"""
from __future__ import annotations

import argparse
import gzip
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# Paths
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "data"
RESULTS_DIR = SCRIPT_DIR / "results"

# Data file paths
CLINVAR_FILE = DATA_DIR / "variant_summary.txt.gz"
GENCODE_GTF = DATA_DIR / "gencode.v47.annotation.gtf.gz"
REFSEQ_GTF = DATA_DIR / "GRCh38_latest_genomic.gtf.gz"

# Download URLs
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz"
REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz"


def download_file(url: str, output_path: Path, description: str) -> bool:
    """Download a file using wget."""
    if output_path.exists():
        print(f"  {description}: already exists")
        return True

    print(f"  Downloading {description}...")
    try:
        subprocess.run(
            ["wget", "-q", "--show-progress", "-O", str(output_path), url],
            check=True
        )
        print(f"  Downloaded: {output_path.name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  Error downloading {description}: {e}")
        return False


def download_data(force: bool = False) -> bool:
    """Download all required data files."""
    print("\n" + "=" * 60)
    print("Step 1: Downloading data files")
    print("=" * 60)

    DATA_DIR.mkdir(parents=True, exist_ok=True)

    if force:
        for f in [CLINVAR_FILE, GENCODE_GTF, REFSEQ_GTF]:
            if f.exists():
                f.unlink()

    success = True
    success &= download_file(CLINVAR_URL, CLINVAR_FILE, "ClinVar variants")
    success &= download_file(GENCODE_URL, GENCODE_GTF, "GENCODE GTF")
    success &= download_file(REFSEQ_URL, REFSEQ_GTF, "RefSeq GTF")

    return success


def build_mane_select_mapping() -> tuple[set, dict, dict]:
    """
    Build MANE Select RefSeq transcript set and gene→UTR genomic interval mapping.

    Returns:
        mane_refseq: Set of MANE Select RefSeq IDs (e.g., "NM_000037")
        utr_lengths: Dict of gene_name → 5'UTR length (in bp)
        utr_intervals: Dict of gene_name → list of (chrom, start, end, strand) tuples
    """
    print("\n" + "=" * 60)
    print("Step 2: Building MANE Select transcript mapping")
    print("=" * 60)

    # First: Get MANE Select Ensembl IDs from GENCODE
    print("  Parsing GENCODE GTF for MANE Select transcripts...")
    mane_ensembl = {}  # ensembl_id -> gene_name
    gene_strands = {}
    cds_starts = {}
    utr_regions = {}  # gene_name -> list of (chrom, start, end)

    open_func = gzip.open if str(GENCODE_GTF).endswith('.gz') else open
    with open_func(GENCODE_GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if 'tag "MANE_Select"' not in line:
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attrs = fields[8]
            chrom = fields[0].replace("chr", "")  # Normalize to match ClinVar

            gene_match = re.search(r'gene_name "([^"]+)"', attrs)
            transcript_match = re.search(r'transcript_id "([^"]+)"', attrs)

            if not gene_match:
                continue

            gene_name = gene_match.group(1)
            strand = fields[6]
            gene_strands[gene_name] = strand

            if feature_type == "transcript" and transcript_match:
                ensembl_id = transcript_match.group(1).split(".")[0]
                mane_ensembl[ensembl_id] = gene_name

            elif feature_type == "start_codon":
                start = int(fields[3])
                end = int(fields[4])
                if gene_name not in cds_starts:
                    cds_starts[gene_name] = start if strand == "+" else end

            elif "UTR" in feature_type:
                start = int(fields[3])
                end = int(fields[4])
                if gene_name not in utr_regions:
                    utr_regions[gene_name] = []
                utr_regions[gene_name].append((chrom, start, end))

    print(f"    Found {len(mane_ensembl)} MANE Select transcripts")

    # Filter to 5' UTR regions only and calculate lengths
    utr_lengths = {}
    utr_intervals = {}  # gene_name -> [(chrom, start, end, strand), ...]

    for gene_name, regions in utr_regions.items():
        if gene_name not in cds_starts or gene_name not in gene_strands:
            continue

        cds_start = cds_starts[gene_name]
        strand = gene_strands[gene_name]
        total_5utr = 0
        gene_5utr_intervals = []

        for chrom, start, end in regions:
            is_5utr = False
            if strand == "+":
                if end < cds_start:  # 5' UTR is before CDS on + strand
                    is_5utr = True
            else:
                if start > cds_start:  # 5' UTR is after CDS on - strand
                    is_5utr = True

            if is_5utr:
                total_5utr += end - start + 1
                gene_5utr_intervals.append((chrom, start, end, strand))

        if total_5utr > 0:
            utr_lengths[gene_name] = total_5utr
            utr_intervals[gene_name] = gene_5utr_intervals

    print(f"    Calculated UTR lengths for {len(utr_lengths)} genes")
    print(f"    Stored genomic intervals for {len(utr_intervals)} genes")

    # Second: Get Ensembl→RefSeq mapping from RefSeq GTF
    print("  Parsing RefSeq GTF for Ensembl cross-references...")
    ensembl_to_refseq = {}

    open_func = gzip.open if str(REFSEQ_GTF).endswith('.gz') else open
    with open_func(REFSEQ_GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue

            attrs = fields[8]
            refseq_match = re.search(r'transcript_id "([^"]+)"', attrs)
            ensembl_match = re.search(r'db_xref "Ensembl:([^"]+)"', attrs)

            if refseq_match and ensembl_match:
                refseq_id = refseq_match.group(1).split(".")[0]
                ensembl_id = ensembl_match.group(1).split(".")[0]

                if refseq_id.startswith("NM_"):
                    ensembl_to_refseq[ensembl_id] = refseq_id

    print(f"    Found {len(ensembl_to_refseq)} Ensembl→RefSeq mappings")

    # Map MANE Ensembl IDs to RefSeq
    mane_refseq = set()
    for ensembl_id in mane_ensembl.keys():
        if ensembl_id in ensembl_to_refseq:
            mane_refseq.add(ensembl_to_refseq[ensembl_id])

    print(f"    Mapped to {len(mane_refseq)} MANE Select RefSeq IDs")

    return mane_refseq, utr_lengths, utr_intervals


def check_utr_containment(
    chrom: str, var_start: int, var_stop: int, utr_intervals: list
) -> tuple[bool, bool, int]:
    """
    Check if variant start and/or stop positions are within 5' UTR region(s).

    Args:
        chrom: Chromosome of the variant (without 'chr' prefix)
        var_start: Genomic start position of variant
        var_stop: Genomic stop position of variant
        utr_intervals: List of (chrom, start, end, strand) tuples for this gene's 5' UTR

    Returns:
        (start_in_utr, stop_in_utr, variant_size): Flags for each position and size in bp
    """
    if not utr_intervals:
        return (False, False, 0)

    # Check if start and stop positions fall within any UTR interval
    start_in_utr = False
    stop_in_utr = False

    for interval_chrom, interval_start, interval_end, strand in utr_intervals:
        if chrom != interval_chrom:
            continue
        # Check if var_start is within this interval
        if interval_start <= var_start <= interval_end:
            start_in_utr = True
        # Check if var_stop is within this interval
        if interval_start <= var_stop <= interval_end:
            stop_in_utr = True

    variant_size = var_stop - var_start + 1
    return (start_in_utr, stop_in_utr, variant_size)


def extract_5utr_variants(mane_refseq: set, utr_intervals: dict) -> pd.DataFrame:
    """
    Extract 5' UTR indel/duplication variants for MANE Select transcripts.

    Uses genomic coordinates to check UTR containment. Returns all variants
    where at least the START position is in the 5' UTR, with flags indicating
    whether start and/or stop are in UTR.
    """
    print("\n" + "=" * 60)
    print("Step 3: Extracting 5' UTR variants from ClinVar")
    print("=" * 60)

    variants = []
    variant_types = {"Deletion", "Duplication", "Insertion"}

    # Build gene name -> UTR intervals lookup (for genes we have UTR data for)
    genes_with_utr = set(utr_intervals.keys())

    # Pattern to extract RefSeq ID and HGVS from name (for reference, not filtering)
    name_pattern = re.compile(r'(NM_\d+)\.?\d*\([^)]+\):([^\s]+)')

    total_indels = 0
    mane_indels = 0
    start_in_utr_count = 0
    both_in_utr_count = 0

    with gzip.open(CLINVAR_FILE, "rt") as f:
        header = f.readline().strip().split("\t")

        # Get column indices
        col_idx = {name: header.index(name) for name in [
            "Type", "Name", "GeneSymbol", "Assembly", "VariationID",
            "Chromosome", "Start", "Stop", "ClinicalSignificance"
        ]}

        for line in f:
            fields = line.strip().split("\t")

            # Filter: GRCh38 only
            if fields[col_idx["Assembly"]] != "GRCh38":
                continue

            # Filter: Indels only
            if fields[col_idx["Type"]] not in variant_types:
                continue

            total_indels += 1

            gene_name = fields[col_idx["GeneSymbol"]]
            chrom = fields[col_idx["Chromosome"]]
            name = fields[col_idx["Name"]]

            # Get genomic coordinates
            try:
                var_start = int(fields[col_idx["Start"]])
                var_stop = int(fields[col_idx["Stop"]])
            except ValueError:
                continue

            # Check if gene has UTR data
            if gene_name not in genes_with_utr:
                continue

            mane_indels += 1

            # Check UTR containment for start and stop positions
            gene_utr_intervals = utr_intervals[gene_name]
            start_in_utr, stop_in_utr, variant_size = check_utr_containment(
                chrom, var_start, var_stop, gene_utr_intervals
            )

            # Include if at least START is in UTR (more permissive filter)
            if not start_in_utr:
                continue

            start_in_utr_count += 1
            if start_in_utr and stop_in_utr:
                both_in_utr_count += 1

            # Extract RefSeq ID and HGVS if available (for reference)
            match = name_pattern.search(name)
            refseq_id = match.group(1) if match else ""
            hgvsc = match.group(2) if match else ""

            variants.append({
                "variant_id": fields[col_idx["VariationID"]],
                "gene_name": gene_name,
                "refseq_id": refseq_id,
                "chromosome": chrom,
                "start": var_start,
                "stop": var_stop,
                "genomic_size": variant_size,
                "start_in_utr": start_in_utr,
                "stop_in_utr": stop_in_utr,
                "both_in_utr": start_in_utr and stop_in_utr,
                "hgvsc": hgvsc,
                "variant_type": fields[col_idx["Type"]].lower(),
                "title": name,
                "clinical_significance": fields[col_idx["ClinicalSignificance"]],
            })

    df = pd.DataFrame(variants)
    if not df.empty:
        df = df.drop_duplicates(subset=["variant_id"])

    print(f"  Total GRCh38 indels scanned: {total_indels}")
    print(f"  Indels in genes with UTR data: {mane_indels}")
    print(f"  Start in 5' UTR: {start_in_utr_count}")
    print(f"  Both start AND stop in 5' UTR: {both_in_utr_count}")
    print(f"  Unique variants (start in UTR): {len(df)}")
    for vtype in ["deletion", "insertion", "duplication"]:
        count = len(df[df["variant_type"] == vtype]) if not df.empty else 0
        print(f"    {vtype.capitalize()}s: {count}")

    return df


def parse_size_from_hgvs(hgvsc: str, variant_type: str) -> int:
    """
    Parse size change from HGVS notation.

    Only counts the portion affecting the 5'UTR. For variants that cross
    into the CDS (e.g., c.-14_10del), only the UTR portion is counted.
    """
    if not hgvsc or not isinstance(hgvsc, str):
        return 0

    # Range pattern: c.-X_-Ydel, c.-X_-Ydup, c.-X_Ydel
    range_match = re.match(r'c\.(-?\d+)_(-?\d+)(del|dup|ins)', hgvsc)
    if range_match:
        start = int(range_match.group(1))
        end = int(range_match.group(2))
        action = range_match.group(3)

        # Calculate UTR-specific size
        # Positions < 0 are in 5'UTR, positions >= 1 are in CDS
        if start < 0 and end < 0:
            # Entirely in 5'UTR
            size = abs(end - start) + 1
        elif start < 0 and end >= 1:
            # Crosses from 5'UTR into CDS - only count UTR portion
            # UTR portion is from start to -1
            size = abs(start)  # e.g., c.-14_10del -> 14bp in UTR
        else:
            # Entirely in CDS (shouldn't happen for 5'UTR variants)
            size = 0

        if action == "del":
            return -size
        elif action == "dup":
            return size
        elif action == "ins":
            ins_match = re.search(r'ins([ACGT]+)', hgvsc)
            return len(ins_match.group(1)) if ins_match else size

    # Single position: c.-138del, c.-14dup
    single_match = re.match(r'c\.(-?\d+)(del|dup|ins)', hgvsc)
    if single_match:
        action = single_match.group(2)
        seq_match = re.search(r'(del|dup|ins)([ACGT]+)', hgvsc)
        if seq_match:
            seq_len = len(seq_match.group(2))
            return -seq_len if action == "del" else seq_len
        return -1 if action == "del" else 1

    # Repeat notation: c.-820_-817G(4_5)
    repeat_match = re.match(r'c\.(-?\d+)_(-?\d+)([ACGT])\((\d+)_(\d+)\)', hgvsc)
    if repeat_match:
        from_count = int(repeat_match.group(4))
        to_count = int(repeat_match.group(5))
        return to_count - from_count

    return -1 if variant_type == "deletion" else 1


def variant_in_utr(hgvsc: str, utr_length: float) -> bool:
    """Check if variant start position is actually within the UTR."""
    if pd.isna(utr_length):
        return False

    # Get the start position from HGVS (first negative number)
    match = re.match(r'c\.(-?\d+)', hgvsc)
    if match:
        start_pos = int(match.group(1))
        # Position must be between -1 and -utr_length to be in UTR
        # c.-1 is last base of UTR (adjacent to start codon)
        # c.-utr_len is first base of UTR
        if start_pos < 0:
            return abs(start_pos) <= utr_length
    return True


def parse_insertion_size(hgvsc: str) -> int:
    """Parse inserted sequence length from HGVS notation."""
    if not hgvsc or not isinstance(hgvsc, str):
        return 1

    # Look for inserted sequence: insACGT
    match = re.search(r'ins([ACGT]+)', hgvsc)
    if match:
        return len(match.group(1))

    return 1  # Default to 1bp if can't parse


def parse_duplication_size(hgvsc: str) -> int:
    """Parse duplicated region size from HGVS notation."""
    if not hgvsc or not isinstance(hgvsc, str):
        return 1

    # Range duplication: c.-X_-Ydup or c.-X_Ydup
    range_match = re.match(r'c\.(-?\d+)_(-?\d+)dup', hgvsc)
    if range_match:
        start = int(range_match.group(1))
        end = int(range_match.group(2))
        # Both in UTR (negative positions)
        if start < 0 and end < 0:
            return abs(start) - abs(end) + 1
        # Crosses into CDS - only count UTR portion
        elif start < 0 and end >= 1:
            return abs(start)
        else:
            return abs(end - start) + 1

    # Duplication with sequence: dupACGT
    seq_match = re.search(r'dup([ACGT]+)', hgvsc)
    if seq_match:
        return len(seq_match.group(1))

    # Single position duplication: c.-14dup
    single_match = re.match(r'c\.(-?\d+)dup', hgvsc)
    if single_match:
        return 1

    return 1  # Default to 1bp if can't parse


def calculate_utr_changes(df: pd.DataFrame, utr_lengths: dict) -> pd.DataFrame:
    """
    Calculate size changes and mutant UTR lengths.

    Uses a hybrid approach:
    - Deletions: Use genomic coordinates (accurate)
    - Insertions: Parse HGVS for inserted sequence length
    - Duplications: Parse HGVS for duplicated range
    """
    print("\n" + "=" * 60)
    print("Step 4: Calculating UTR length changes")
    print("=" * 60)

    if df.empty:
        print("  No variants to process")
        return df

    df = df.copy()

    # Add WT UTR lengths
    df["wt_utr_length"] = df["gene_name"].map(utr_lengths)

    # Filter variants without UTR data (shouldn't happen with new logic, but safety check)
    df = df[df["wt_utr_length"].notna()].copy()
    print(f"  Variants with UTR data: {len(df)}")

    # Calculate size changes using hybrid approach
    def calc_size_change(row):
        vtype = row["variant_type"]
        if vtype == "deletion":
            # Deletions: genomic coordinates are accurate
            return -row["genomic_size"]
        elif vtype == "insertion":
            # Insertions: parse HGVS for inserted sequence length
            return parse_insertion_size(row["hgvsc"])
        elif vtype == "duplication":
            # Duplications: parse HGVS for duplicated range
            return parse_duplication_size(row["hgvsc"])
        else:
            return 0

    df["size_change"] = df.apply(calc_size_change, axis=1)

    # Calculate mutant UTR length
    df["mutant_utr_length"] = df["wt_utr_length"] + df["size_change"]
    df.loc[df["mutant_utr_length"] < 0, "mutant_utr_length"] = 0

    print(f"  Mean WT UTR: {df['wt_utr_length'].mean():.1f} bp")
    print(f"  Mean size change: {df['size_change'].mean():.1f} bp")

    return df


def crosses_into_cds(hgvsc: str) -> bool:
    """Check if variant crosses from UTR into CDS (affects start codon)."""
    match = re.match(r'c\.(-?\d+)_(-?\d+)', hgvsc)
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        return start < 0 and end >= 1
    return False


def crosses_threshold(row: pd.Series, threshold: int) -> bool:
    """Check if a variant crosses a specific threshold."""
    wt = row["wt_utr_length"]
    mut = row["mutant_utr_length"]
    vtype = row["variant_type"]

    if vtype == "deletion":
        # Deletion crosses if: WT > threshold AND Mut <= threshold
        return wt > threshold and mut <= threshold
    else:
        # Insertion/duplication crosses if: WT < threshold AND Mut >= threshold
        return wt < threshold and mut >= threshold


def filter_variants(df: pd.DataFrame, min_size: int, thresholds: list[int]) -> tuple[pd.DataFrame, dict]:
    """
    Filter by min size and check multiple threshold crossings.

    Args:
        df: Input dataframe
        min_size: Minimum size change to include
        thresholds: List of thresholds to check (e.g., [10, 20, 30, 40, 50])

    Returns:
        df_sized: Variants with size >= min_size, with threshold columns added
        crossing_dfs: Dict of threshold -> dataframe of variants crossing that threshold
    """
    print("\n" + "=" * 60)
    print(f"Step 5: Filtering (min size: {min_size}bp)")
    print("=" * 60)

    # Filter by minimum size
    df_sized = df[df["size_change"].abs() >= min_size].copy()
    print(f"  Variants with size >= {min_size}bp: {len(df_sized)}")

    # Check each threshold and add columns
    crossed_thresholds = []
    crossing_dfs = {}

    print(f"\n  Checking thresholds: {thresholds}")
    for threshold in thresholds:
        col_name = f"crosses_{threshold}bp"
        df_sized[col_name] = df_sized.apply(lambda row: crosses_threshold(row, threshold), axis=1)

        # Get variants crossing this threshold
        df_crossing = df_sized[df_sized[col_name]].copy()
        crossing_dfs[threshold] = df_crossing

        # Count by type
        del_count = len(df_crossing[df_crossing["variant_type"] == "deletion"])
        ins_count = len(df_crossing[df_crossing["variant_type"] == "insertion"])
        dup_count = len(df_crossing[df_crossing["variant_type"] == "duplication"])

        print(f"    {threshold}bp: {len(df_crossing)} variants (del: {del_count}, ins: {ins_count}, dup: {dup_count})")

    # Add summary column showing all thresholds crossed
    def get_crossed_thresholds(row):
        crossed = []
        for threshold in thresholds:
            if row[f"crosses_{threshold}bp"]:
                crossed.append(str(threshold))
        return ",".join(crossed) if crossed else ""

    df_sized["thresholds_crossed"] = df_sized.apply(get_crossed_thresholds, axis=1)

    return df_sized, crossing_dfs


def format_output(df: pd.DataFrame) -> pd.DataFrame:
    """Format output DataFrame for clean CSV export."""
    df = df.copy()

    # Add ClinVar URL
    df["clinvar_url"] = df["variant_id"].apply(
        lambda x: f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{x}/"
    )

    # Convert numeric columns to integers
    df["wt_utr_length"] = df["wt_utr_length"].astype(int)
    df["mutant_utr_length"] = df["mutant_utr_length"].astype(int)
    df["size_change"] = df["size_change"].astype(int)

    # Create clean clinical significance categories for sorting
    def clin_sig_order(sig):
        sig_lower = str(sig).lower()
        if "pathogenic" in sig_lower and "likely" not in sig_lower:
            return 0
        elif "likely pathogenic" in sig_lower:
            return 1
        elif "uncertain" in sig_lower:
            return 2
        elif "conflicting" in sig_lower:
            return 3
        elif "likely benign" in sig_lower:
            return 4
        elif "benign" in sig_lower:
            return 5
        return 6

    # Sort by size_change (smallest/most negative to largest/most positive)
    df = df.sort_values("size_change")

    # Add genomic coordinates column for clarity
    df["genomic_coords"] = df.apply(
        lambda row: f"chr{row['chromosome']}:{row['start']}-{row['stop']}", axis=1
    )

    # Base columns
    output_cols = [
        "gene_name",
        "hgvsc",
        "variant_type",
        "wt_utr_length",
        "mutant_utr_length",
        "size_change",
    ]

    # Add threshold columns if present
    threshold_cols = [col for col in df.columns if col.startswith("crosses_")]
    if threshold_cols:
        output_cols.extend(sorted(threshold_cols, key=lambda x: int(x.split("_")[1].replace("bp", ""))))

    # Add thresholds_crossed summary if present
    if "thresholds_crossed" in df.columns:
        output_cols.append("thresholds_crossed")

    # Add remaining columns
    output_cols.extend([
        "genomic_coords",
        "clinical_significance",
        "clinvar_url",
        "variant_id",
        "refseq_id",
        "chromosome",
        "start",
        "stop",
    ])

    # Only include columns that exist
    output_cols = [col for col in output_cols if col in df.columns]

    return df[output_cols]


def create_plot(df: pd.DataFrame, threshold: int, output_file: Path):
    """Create slope plot showing WT vs mutant UTR lengths."""
    print("\n" + "=" * 60)
    print("Step 6: Generating plot")
    print("=" * 60)

    if len(df) == 0:
        print("  No variants to plot!")
        return

    fig, ax = plt.subplots(figsize=(10, 8))

    colors = {"deletion": "#E74C3C", "insertion": "#27AE60", "duplication": "#3498DB"}

    for _, row in df.iterrows():
        color = colors.get(row["variant_type"], "gray")
        ax.plot([0, 1], [row["wt_utr_length"], row["mutant_utr_length"]],
                color=color, alpha=0.7, linewidth=1.5)
        ax.scatter([0, 1], [row["wt_utr_length"], row["mutant_utr_length"]],
                   color=color, s=50, alpha=0.7)

    ax.axhline(y=threshold, color="black", linestyle="--", linewidth=2)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["WT UTR Length", "Mutant UTR Length"], fontsize=12)
    ax.set_ylabel("5' UTR Length (bp)", fontsize=12)
    ax.set_title(f"5' UTR Length Changes Crossing {threshold}bp Threshold\n({len(df)} variants)",
                 fontsize=14)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color="#E74C3C", linewidth=2, label="Deletion"),
        Line2D([0], [0], color="#27AE60", linewidth=2, label="Insertion"),
        Line2D([0], [0], color="#3498DB", linewidth=2, label="Duplication"),
        Line2D([0], [0], color="black", linestyle="--", linewidth=2, label=f"Threshold ({threshold}bp)")
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    ax.set_xlim(-0.2, 1.2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"  Saved: {output_file}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze 5' UTR length changes in ClinVar variants"
    )
    parser.add_argument(
        "--min-size", type=int, default=10,
        help="Minimum variant size in bp (default: 10)"
    )
    parser.add_argument(
        "--thresholds", type=str, default="10,20,30,40,50",
        help="Comma-separated UTR length thresholds in bp (default: 10,20,30,40,50)"
    )
    parser.add_argument(
        "--force-download", action="store_true",
        help="Force re-download of data files"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Parse thresholds
    thresholds = [int(t.strip()) for t in args.thresholds.split(",")]

    print("=" * 60)
    print("5' UTR Length Change Analysis (Genomic Coordinate Validation)")
    print("=" * 60)
    print(f"  Min size filter: {args.min_size} bp")
    print(f"  Thresholds: {thresholds} bp")
    print("  Filtering: Both variant start AND stop must be in 5' UTR")

    # Create output directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Download data
    if not download_data(force=args.force_download):
        print("\nError: Failed to download required data files")
        sys.exit(1)

    # Step 2: Build MANE Select mapping
    mane_refseq, utr_lengths, utr_intervals = build_mane_select_mapping()

    # Step 3: Extract 5' UTR variants (using genomic coordinate validation)
    df_all = extract_5utr_variants(mane_refseq, utr_intervals)

    if df_all.empty:
        print("\nNo 5' UTR variants found!")
        sys.exit(1)

    # Save intermediate: all variants where START is in UTR
    start_in_utr_output = RESULTS_DIR / "5utr_variants_start_in_utr.csv"
    df_all.to_csv(start_in_utr_output, index=False)
    print(f"\n  Saved start-in-UTR dataset: {start_in_utr_output} ({len(df_all)} variants)")

    # Filter to variants where BOTH start AND stop are in UTR
    df_both = df_all[df_all["both_in_utr"]].copy()
    both_in_utr_output = RESULTS_DIR / "5utr_variants_both_in_utr.csv"
    df_both.to_csv(both_in_utr_output, index=False)
    print(f"  Saved both-in-UTR dataset: {both_in_utr_output} ({len(df_both)} variants)")

    if df_both.empty:
        print("\nNo variants with both start AND stop in UTR!")
        sys.exit(1)

    # Step 4: Calculate UTR changes (only for variants fully in UTR)
    df = calculate_utr_changes(df_both, utr_lengths)

    # Save full results with UTR calculations
    full_output = RESULTS_DIR / "5utr_variants_full.csv"
    df.to_csv(full_output, index=False)
    print(f"  Saved full dataset with UTR calculations: {full_output}")

    # Step 5: Filter variants with multiple thresholds
    df_sized, crossing_dfs = filter_variants(df, args.min_size, thresholds)

    # Save size-filtered results (with all threshold columns)
    sized_output = RESULTS_DIR / f"5utr_variants_min{args.min_size}bp.csv"
    df_sized_formatted = format_output(df_sized) if not df_sized.empty else df_sized
    df_sized_formatted.to_csv(sized_output, index=False)
    print(f"  Saved size-filtered dataset: {sized_output} ({len(df_sized)} variants)")

    # Save results for each threshold crossing
    saved_files = []
    for threshold in thresholds:
        df_crossing = crossing_dfs[threshold]
        if not df_crossing.empty:
            crossing_output = RESULTS_DIR / f"5utr_variants_min{args.min_size}bp_cross{threshold}bp.csv"
            df_crossing_formatted = format_output(df_crossing)
            df_crossing_formatted.to_csv(crossing_output, index=False)
            saved_files.append((threshold, crossing_output, len(df_crossing)))
            print(f"  Saved {threshold}bp threshold dataset: {crossing_output}")

            # Create plot for this threshold
            plot_output = RESULTS_DIR / f"5utr_length_changes_min{args.min_size}bp_cross{threshold}bp.png"
            create_plot(df_crossing, threshold, plot_output)

    # Final summary
    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    print(f"\nOutput files in {RESULTS_DIR}/:")
    print(f"  - 5utr_variants_start_in_utr.csv ({len(df_all)} variants)")
    print(f"  - 5utr_variants_both_in_utr.csv ({len(df_both)} variants)")
    print(f"  - 5utr_variants_full.csv ({len(df)} variants)")
    print(f"  - {sized_output.name} ({len(df_sized)} variants)")
    for threshold, filepath, count in saved_files:
        print(f"  - {filepath.name} ({count} variants)")
        print(f"  - 5utr_length_changes_min{args.min_size}bp_cross{threshold}bp.png")

    # Summary table for collaborator - filtering funnel
    print("\n" + "=" * 60)
    print("FILTERING FUNNEL")
    print("=" * 60)

    print(f"\n[Step 1] START position in 5' UTR: {len(df_all)}")
    print(f"         (Variant start coordinate within GENCODE 5' UTR region)")
    print(f"         Deletions:    {len(df_all[df_all['variant_type'] == 'deletion'])}")
    print(f"         Insertions:   {len(df_all[df_all['variant_type'] == 'insertion'])}")
    print(f"         Duplications: {len(df_all[df_all['variant_type'] == 'duplication'])}")

    print(f"\n[Step 2] BOTH start AND stop in 5' UTR: {len(df_both)}")
    print(f"         (Both coordinates within GENCODE 5' UTR regions)")
    print(f"         Deletions:    {len(df_both[df_both['variant_type'] == 'deletion'])}")
    print(f"         Insertions:   {len(df_both[df_both['variant_type'] == 'insertion'])}")
    print(f"         Duplications: {len(df_both[df_both['variant_type'] == 'duplication'])}")

    print(f"\n[Step 3] Size change >= {args.min_size}bp: {len(df_sized)}")
    print(f"         Deletions:    {len(df_sized[df_sized['variant_type'] == 'deletion'])}")
    print(f"         Insertions:   {len(df_sized[df_sized['variant_type'] == 'insertion'])}")
    print(f"         Duplications: {len(df_sized[df_sized['variant_type'] == 'duplication'])}")

    print(f"\n[Step 4] Threshold crossings:")
    for threshold in thresholds:
        df_crossing = crossing_dfs[threshold]
        del_count = len(df_crossing[df_crossing['variant_type'] == 'deletion'])
        ins_count = len(df_crossing[df_crossing['variant_type'] == 'insertion'])
        dup_count = len(df_crossing[df_crossing['variant_type'] == 'duplication'])
        print(f"         {threshold}bp: {len(df_crossing)} variants (del: {del_count}, ins: {ins_count}, dup: {dup_count})")

    # Show variants that cross any threshold
    any_crossing = df_sized[df_sized["thresholds_crossed"] != ""]
    if not any_crossing.empty:
        print(f"\nVariants crossing any threshold ({len(any_crossing)}):")
        display_cols = ["gene_name", "hgvsc", "variant_type", "wt_utr_length",
                        "mutant_utr_length", "thresholds_crossed", "clinical_significance"]
        print(format_output(any_crossing)[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
