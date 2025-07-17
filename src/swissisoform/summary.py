#!/usr/bin/env python3
"""Summary analysis backend for SwissIsoform pipeline results."""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import warnings

# Suppress pandas warnings
warnings.filterwarnings("ignore", category=FutureWarning)


class SummaryAnalyzer:
    """Class to handle summary analysis of SwissIsoform pipeline results."""

    def __init__(self):
        """Initialize the analyzer."""
        pass

    def dataset_has_data(self, dataset):
        """Check if a dataset has any data available for analysis."""
        # Check for mutation data
        mutation_gene_file = Path(
            f"../results/{dataset}/mutations/gene_level_results.csv"
        )

        # Check for localization data
        loc_files = [
            Path(
                f"../results/{dataset}/localization/protein_sequences_pairs_Accurate_results.csv"
            ),
            Path(
                f"../results/{dataset}/localization/protein_sequences_pairs_Fast_results.csv"
            ),
        ]

        has_mutation_data = mutation_gene_file.exists()
        has_localization_data = any(f.exists() for f in loc_files)

        return has_mutation_data or has_localization_data

    def load_mutation_results(self, dataset):
        """Load mutation analysis results for a dataset."""
        base_path = Path(f"../results/{dataset}/mutations")

        gene_results = None
        pair_results = None

        gene_file = base_path / "gene_level_results.csv"
        pair_file = base_path / "truncation_level_results.csv"

        if gene_file.exists():
            gene_results = pd.read_csv(gene_file)
            print(
                f"Loaded {len(gene_results)} gene-level mutation results for {dataset}"
            )

        if pair_file.exists():
            pair_results = pd.read_csv(pair_file)
            print(
                f"Loaded {len(pair_results)} transcript-truncation pair results for {dataset}"
            )

        return gene_results, pair_results

    def load_protein_sequences(self, dataset):
        """Load protein sequence data to get variant metadata."""
        base_path = Path(f"../results/{dataset}/proteins")

        pairs_file = base_path / "protein_sequences_pairs.csv"
        mutations_file = base_path / "protein_sequences_with_mutations.csv"

        protein_data = {}

        if pairs_file.exists():
            pairs_df = pd.read_csv(pairs_file)
            print(f"Loaded {len(pairs_df)} protein sequence pairs for {dataset}")
            # Index by sequence identifier
            for _, row in pairs_df.iterrows():
                seq_id = f"{row['gene']}_{row['transcript_id']}_{row['variant_id']}"
                protein_data[seq_id] = row.to_dict()

        if mutations_file.exists():
            mutations_df = pd.read_csv(mutations_file)
            print(
                f"Loaded {len(mutations_df)} protein sequence mutations for {dataset}"
            )
            # Index by sequence identifier
            for _, row in mutations_df.iterrows():
                seq_id = f"{row['gene']}_{row['transcript_id']}_{row['variant_id']}"
                protein_data[seq_id] = row.to_dict()

        return protein_data

    def load_localization_results(self, dataset):
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
                    print(
                        f"Loaded {len(df)} {key} localization predictions for {dataset}"
                    )

                except Exception as e:
                    print(f"Error loading {filename}: {e}")
                    continue
            else:
                print(f"Missing {key} localization predictions for {dataset}")

        return results

    def analyze_mutations(self, dataset, gene_results, pair_results):
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
                    mutation_type = (
                        col.replace("mutations_", "").replace("_", " ").title()
                    )

                    # Count genes that have at least one mutation of this type
                    genes_with_mutation = pair_results[pair_results[col] > 0][
                        "gene_name"
                    ].nunique()
                    total_mutations = pair_results[col].sum()

                    gene_mutation_counts[mutation_type] = genes_with_mutation
                    summary_lines.append(
                        f"  {mutation_type}: {genes_with_mutation} genes ({total_mutations} total mutations)"
                    )

                # Special combined category: frameshift OR nonsense (changed from "Truncating")
                frameshift_col = "mutations_frameshift_variant"
                nonsense_col = "mutations_nonsense_variant"

                if (
                    frameshift_col in pair_results.columns
                    and nonsense_col in pair_results.columns
                ):
                    genes_with_frameshift_or_nonsense = pair_results[
                        (pair_results[frameshift_col] > 0)
                        | (pair_results[nonsense_col] > 0)
                    ]["gene_name"].nunique()

                    # Calculate total frameshift + nonsense mutations
                    total_frameshift_or_nonsense_mutations = (
                        pair_results[frameshift_col].sum()
                        + pair_results[nonsense_col].sum()
                    )

                    summary_lines.append(
                        f"  Frameshift OR Nonsense: {genes_with_frameshift_or_nonsense} genes ({total_frameshift_or_nonsense_mutations} total mutations)"
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

    def parse_sequence_id(self, seq_id):
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

    def normalize_localization(self, localization):
        """Normalize localization by sorting multiple locations."""
        if pd.isna(localization) or not localization:
            return "Unknown"

        localization = str(localization).strip()

        # If multiple localizations, sort them for consistent comparison
        if "|" in localization:
            locations = [loc.strip() for loc in localization.split("|")]
            locations.sort()  # Sort alphabetically for consistency
            return "|".join(locations)

        return localization

    def are_localizations_different(self, loc1, loc2):
        """Check if two localizations are different, accounting for multi-location predictions."""
        if pd.isna(loc1) or pd.isna(loc2):
            return True

        # Normalize both localizations
        norm_loc1 = self.normalize_localization(loc1)
        norm_loc2 = self.normalize_localization(loc2)

        return norm_loc1 != norm_loc2

    def analyze_localizations(self, dataset, loc_results, protein_data):
        """Analyze localization predictions for a specific dataset."""
        print(f"\n=== ANALYZING LOCALIZATIONS FOR {dataset.upper()} DATASET ===")

        summary_lines = []
        summary_lines.append(
            f"LOCALIZATION ANALYSIS SUMMARY - {dataset.upper()} DATASET"
        )
        summary_lines.append("=" * 60)

        all_localization_comparisons = []

        if not loc_results:
            summary_lines.append("No localization data available")
            return summary_lines, all_localization_comparisons

        # Prefer accurate results over fast results, but ensure consistency
        pairs_data = None
        mutations_data = None
        model_type = "Unknown"

        # First, check what's available and prioritize Accurate if both pairs and mutations exist
        has_pairs_accurate = "pairs_accurate" in loc_results
        has_pairs_fast = "pairs_fast" in loc_results
        has_mutations_accurate = "mutations_accurate" in loc_results
        has_mutations_fast = "mutations_fast" in loc_results

        # Prioritize Accurate model if available for both, otherwise use Fast
        if has_pairs_accurate and (has_mutations_accurate or not has_mutations_fast):
            pairs_data = loc_results["pairs_accurate"]
            model_type = "Accurate"
            if has_mutations_accurate:
                mutations_data = loc_results["mutations_accurate"]
            # If no mutations data at all, that's fine
        elif has_pairs_fast:
            pairs_data = loc_results["pairs_fast"]
            model_type = "Fast"
            if has_mutations_fast:
                mutations_data = loc_results["mutations_fast"]

        print(f"Using {model_type} model for localization analysis")
        if pairs_data is not None and mutations_data is not None:
            print(f"Both pairs and mutations data available with {model_type} model")
        elif pairs_data is not None:
            print(f"Only pairs data available with {model_type} model")
        elif mutations_data is not None:
            print(f"Only mutations data available with {model_type} model")

        if pairs_data is None:
            summary_lines.append("No localization data available")
            return summary_lines, all_localization_comparisons

        summary_lines.append(f"Using {model_type} model predictions")

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
            seq_id = row.get("Sequence_ID", "")
            localization = row.get("Localisation", "")

            if seq_id and localization:
                gene, transcript, variant = self.parse_sequence_id(seq_id)
                if gene and transcript and variant:
                    # Get additional metadata from protein data
                    protein_info = protein_data.get(seq_id, {})

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
                            "protein_info": protein_info,
                        }
                    )

        pairs_df = pd.DataFrame(pairs_parsed)

        if not pairs_df.empty:
            # Count total gene pairs analyzed
            total_genes = pairs_df["gene"].nunique()

            # Group by gene and analyze canonical vs truncated
            gene_localization_changes = []
            localization_change_details = []

            for gene in pairs_df["gene"].unique():
                gene_data = pairs_df[pairs_df["gene"] == gene]

                # Find canonical and truncated sequences
                canonical_seqs = gene_data[gene_data["variant"] == "canonical"]
                truncated_seqs = gene_data[gene_data["variant"].str.startswith("trunc")]

                if not canonical_seqs.empty and not truncated_seqs.empty:
                    canonical_loc = canonical_seqs.iloc[0]["prediction"]
                    canonical_info = canonical_seqs.iloc[0]["protein_info"]

                    for _, trunc_seq in truncated_seqs.iterrows():
                        trunc_loc = trunc_seq["prediction"]
                        trunc_info = trunc_seq["protein_info"]

                        if self.are_localizations_different(canonical_loc, trunc_loc):
                            # Merge protein info for comprehensive record
                            change_record = {
                                "dataset": dataset,
                                "gene": gene,
                                "transcript": trunc_seq["transcript"],
                                "canonical_localization": canonical_loc,
                                "truncated_localization": trunc_loc,
                                "truncated_variant": trunc_seq["variant"],
                                "localization_change_type": "canonical_vs_truncated",
                                "model_type": model_type,
                                # Add protein sequence info
                                "canonical_length": canonical_info.get("length", ""),
                                "truncated_length": trunc_info.get("length", ""),
                                "variant_type": trunc_info.get("variant_type", ""),
                                "mutation_position": trunc_info.get(
                                    "mutation_position", ""
                                ),
                                "mutation_change": trunc_info.get(
                                    "mutation_change", ""
                                ),
                                "aa_change": trunc_info.get("aa_change", ""),
                                "hgvsc": trunc_info.get("hgvsc", ""),
                                "hgvsp": trunc_info.get("hgvsp", ""),
                                "mutation_impact": trunc_info.get(
                                    "mutation_impact", ""
                                ),
                                "mutation_source": trunc_info.get(
                                    "mutation_source", ""
                                ),
                                "clinvar_variant_id": trunc_info.get(
                                    "clinvar_variant_id", ""
                                ),
                            }

                            gene_localization_changes.append(change_record)

                            # Store change details for summary
                            localization_change_details.append(
                                {
                                    "gene": gene,
                                    "from": canonical_loc,
                                    "to": trunc_loc,
                                    "type": "canonical_vs_truncated",
                                }
                            )

            # Count unique genes with localization changes
            unique_genes_with_changes = len(
                set(change["gene"] for change in gene_localization_changes)
            )

            summary_lines.append(f"Total gene pairs analyzed: {total_genes}")
            summary_lines.append(
                f"Genes with canonical vs truncated localization changes: {unique_genes_with_changes} out of {total_genes}"
            )

            all_localization_comparisons.extend(gene_localization_changes)

            # Analyze mutations if available
            if mutations_data is not None:
                mutations_parsed = []
                for _, row in mutations_data.iterrows():
                    seq_id = row.get("Sequence_ID", "")
                    localization = row.get("Localisation", "")

                    if seq_id and localization:
                        gene, transcript, variant = self.parse_sequence_id(seq_id)
                        if gene and transcript and variant:
                            protein_info = protein_data.get(seq_id, {})

                            mutations_parsed.append(
                                {
                                    "dataset": dataset,
                                    "gene": gene,
                                    "transcript": transcript,
                                    "variant": variant,
                                    "sequence_id": seq_id,
                                    "prediction": localization,
                                    "confidence": 0,
                                    "source": "mutations",
                                    "protein_info": protein_info,
                                }
                            )

                mutations_df = pd.DataFrame(mutations_parsed)

                if not mutations_df.empty:
                    # Find genes with missense mutations causing localization changes
                    missense_changes = []
                    missense_change_details = []

                    for gene in mutations_df["gene"].unique():
                        gene_mut_data = mutations_df[mutations_df["gene"] == gene]

                        # Find canonical sequence in pairs data
                        gene_pairs_data = pairs_df[pairs_df["gene"] == gene]
                        canonical_pairs = gene_pairs_data[
                            gene_pairs_data["variant"] == "canonical"
                        ]

                        if not canonical_pairs.empty:
                            canonical_loc = canonical_pairs.iloc[0]["prediction"]
                            canonical_info = canonical_pairs.iloc[0]["protein_info"]

                            # Look for missense mutations with different localizations
                            missense_variants = gene_mut_data[
                                gene_mut_data["variant"].str.contains("mut", na=False)
                            ]

                            for _, mut_seq in missense_variants.iterrows():
                                mut_loc = mut_seq["prediction"]
                                mut_info = mut_seq["protein_info"]

                                if self.are_localizations_different(
                                    canonical_loc, mut_loc
                                ):
                                    change_record = {
                                        "dataset": dataset,
                                        "gene": gene,
                                        "transcript": mut_seq["transcript"],
                                        "canonical_localization": canonical_loc,
                                        "mutated_localization": mut_loc,
                                        "mutated_variant": mut_seq["variant"],
                                        "localization_change_type": "canonical_vs_missense",
                                        "model_type": model_type,
                                        # Add protein sequence info
                                        "canonical_length": canonical_info.get(
                                            "length", ""
                                        ),
                                        "mutated_length": mut_info.get("length", ""),
                                        "variant_type": mut_info.get(
                                            "variant_type", ""
                                        ),
                                        "mutation_position": mut_info.get(
                                            "mutation_position", ""
                                        ),
                                        "mutation_change": mut_info.get(
                                            "mutation_change", ""
                                        ),
                                        "aa_change": mut_info.get("aa_change", ""),
                                        "hgvsc": mut_info.get("hgvsc", ""),
                                        "hgvsp": mut_info.get("hgvsp", ""),
                                        "mutation_impact": mut_info.get(
                                            "mutation_impact", ""
                                        ),
                                        "mutation_source": mut_info.get(
                                            "mutation_source", ""
                                        ),
                                        "clinvar_variant_id": mut_info.get(
                                            "clinvar_variant_id", ""
                                        ),
                                    }

                                    missense_changes.append(change_record)

                                    # Store change details for summary
                                    missense_change_details.append(
                                        {
                                            "gene": gene,
                                            "from": canonical_loc,
                                            "to": mut_loc,
                                            "type": "canonical_vs_missense",
                                        }
                                    )

                    unique_genes_with_missense_changes = len(
                        set(change["gene"] for change in missense_changes)
                    )
                    summary_lines.append(
                        f"Genes with missense mutation localization changes: {unique_genes_with_missense_changes} out of {total_genes}"
                    )
                    all_localization_comparisons.extend(missense_changes)

            # Show the actual localization changes instead of just top localizations
            if localization_change_details or (
                mutations_data is not None and missense_change_details
            ):
                summary_lines.append(f"\nLocalization changes observed:")

                # Combine all changes and count them
                all_changes = localization_change_details + (
                    missense_change_details if mutations_data is not None else []
                )

                # Count transitions
                transition_counts = {}
                for change in all_changes:
                    transition = f"{change['from']} → {change['to']}"
                    transition_counts[transition] = (
                        transition_counts.get(transition, 0) + 1
                    )

                # Show top transitions
                sorted_transitions = sorted(
                    transition_counts.items(), key=lambda x: x[1], reverse=True
                )
                for transition, count in sorted_transitions[:5]:  # Show top 5
                    summary_lines.append(f"  {transition}: {count} changes")
            else:
                summary_lines.append(f"\nNo localization changes observed")

        return summary_lines, all_localization_comparisons

    def create_detailed_localization_analysis(self, dataset, loc_results, protein_data):
        """Create a detailed analysis of all localization predictions for a specific dataset."""
        print(
            f"\n=== CREATING DETAILED LOCALIZATION ANALYSIS FOR {dataset.upper()} DATASET ==="
        )

        detailed_results = []

        if not loc_results:
            return pd.DataFrame(detailed_results)

        # Process pairs data with consistent model selection
        pairs_data = None
        mutations_data = None
        model_type = "Unknown"

        # First, check what's available and prioritize Accurate if both pairs and mutations exist
        has_pairs_accurate = "pairs_accurate" in loc_results
        has_pairs_fast = "pairs_fast" in loc_results
        has_mutations_accurate = "mutations_accurate" in loc_results
        has_mutations_fast = "mutations_fast" in loc_results

        # Prioritize Accurate model if available for both, otherwise use Fast
        if has_pairs_accurate and (has_mutations_accurate or not has_mutations_fast):
            pairs_data = loc_results["pairs_accurate"]
            model_type = "Accurate"
            if has_mutations_accurate:
                mutations_data = loc_results["mutations_accurate"]
        elif has_pairs_fast:
            pairs_data = loc_results["pairs_fast"]
            model_type = "Fast"
            if has_mutations_fast:
                mutations_data = loc_results["mutations_fast"]

        # Process pairs data
        if pairs_data is not None:
            for _, row in pairs_data.iterrows():
                seq_id = row.get("Sequence_ID", "")
                localization = row.get("Localisation", "")

                if seq_id and localization:
                    gene, transcript, variant = self.parse_sequence_id(seq_id)
                    if gene and transcript and variant:
                        protein_info = protein_data.get(seq_id, {})

                        result_record = {
                            "dataset": dataset,
                            "gene": gene,
                            "transcript": transcript,
                            "variant": variant,
                            "sequence_id": seq_id,
                            "prediction": localization,
                            "confidence": 0,
                            "sequence_type": "canonical"
                            if variant == "canonical"
                            else "truncated",
                            "model_type": model_type,
                            "analysis_type": "pairs",
                            # Add protein sequence metadata
                            "sequence_length": protein_info.get("length", ""),
                            "variant_type": protein_info.get("variant_type", ""),
                            "mutation_position": protein_info.get(
                                "mutation_position", ""
                            ),
                            "mutation_change": protein_info.get("mutation_change", ""),
                            "aa_change": protein_info.get("aa_change", ""),
                            "hgvsc": protein_info.get("hgvsc", ""),
                            "hgvsp": protein_info.get("hgvsp", ""),
                            "mutation_impact": protein_info.get("mutation_impact", ""),
                            "mutation_source": protein_info.get("mutation_source", ""),
                            "clinvar_variant_id": protein_info.get(
                                "clinvar_variant_id", ""
                            ),
                        }

                        detailed_results.append(result_record)

        # Process mutations data
        if mutations_data is not None:
            for _, row in mutations_data.iterrows():
                seq_id = row.get("Sequence_ID", "")
                localization = row.get("Localisation", "")

                if seq_id and localization:
                    gene, transcript, variant = self.parse_sequence_id(seq_id)
                    if gene and transcript and variant:
                        protein_info = protein_data.get(seq_id, {})

                        sequence_type = (
                            "canonical"
                            if variant == "canonical"
                            else ("mutated" if "mut" in variant else "truncated")
                        )

                        result_record = {
                            "dataset": dataset,
                            "gene": gene,
                            "transcript": transcript,
                            "variant": variant,
                            "sequence_id": seq_id,
                            "prediction": localization,
                            "confidence": 0,
                            "sequence_type": sequence_type,
                            "model_type": model_type,
                            "analysis_type": "mutations",
                            # Add protein sequence metadata
                            "sequence_length": protein_info.get("length", ""),
                            "variant_type": protein_info.get("variant_type", ""),
                            "mutation_position": protein_info.get(
                                "mutation_position", ""
                            ),
                            "mutation_change": protein_info.get("mutation_change", ""),
                            "aa_change": protein_info.get("aa_change", ""),
                            "hgvsc": protein_info.get("hgvsc", ""),
                            "hgvsp": protein_info.get("hgvsp", ""),
                            "mutation_impact": protein_info.get("mutation_impact", ""),
                            "mutation_source": protein_info.get("mutation_source", ""),
                            "clinvar_variant_id": protein_info.get(
                                "clinvar_variant_id", ""
                            ),
                        }

                        detailed_results.append(result_record)

        return pd.DataFrame(detailed_results)

    def create_gene_summary(
        self, dataset, loc_results, protein_data, pair_results=None
    ):
        """Create a gene-level summary showing prioritized targets."""
        print(f"\n=== CREATING GENE-LEVEL SUMMARY FOR {dataset.upper()} DATASET ===")

        if not loc_results:
            return pd.DataFrame()

        # Get the primary localization data with consistent model selection
        pairs_data = None
        mutations_data = None
        model_type = "Unknown"

        # First, check what's available and prioritize Accurate if both pairs and mutations exist
        has_pairs_accurate = "pairs_accurate" in loc_results
        has_pairs_fast = "pairs_fast" in loc_results
        has_mutations_accurate = "mutations_accurate" in loc_results
        has_mutations_fast = "mutations_fast" in loc_results

        # Prioritize Accurate model if available for both, otherwise use Fast
        if has_pairs_accurate and (has_mutations_accurate or not has_mutations_fast):
            pairs_data = loc_results["pairs_accurate"]
            model_type = "Accurate"
            if has_mutations_accurate:
                mutations_data = loc_results["mutations_accurate"]
        elif has_pairs_fast:
            pairs_data = loc_results["pairs_fast"]
            model_type = "Fast"
            if has_mutations_fast:
                mutations_data = loc_results["mutations_fast"]

        if pairs_data is None:
            return pd.DataFrame()

        # Parse pairs data
        pairs_parsed = []
        for _, row in pairs_data.iterrows():
            seq_id = row.get("Sequence_ID", "")
            localization = row.get("Localisation", "")

            if seq_id and localization:
                gene, transcript, variant = self.parse_sequence_id(seq_id)
                if gene and transcript and variant:
                    pairs_parsed.append(
                        {
                            "gene": gene,
                            "transcript": transcript,
                            "variant": variant,
                            "prediction": localization,
                            "seq_id": seq_id,
                        }
                    )

        pairs_df = pd.DataFrame(pairs_parsed)

        # Parse mutations data if available
        mutations_parsed = []
        if mutations_data is not None:
            for _, row in mutations_data.iterrows():
                seq_id = row.get("Sequence_ID", "")
                localization = row.get("Localisation", "")

                if seq_id and localization:
                    gene, transcript, variant = self.parse_sequence_id(seq_id)
                    if gene and transcript and variant:
                        mutations_parsed.append(
                            {
                                "gene": gene,
                                "transcript": transcript,
                                "variant": variant,
                                "prediction": localization,
                                "seq_id": seq_id,
                            }
                        )

        mutations_df = pd.DataFrame(mutations_parsed)

        # Create gene-level summary
        gene_summaries = []

        for gene in pairs_df["gene"].unique():
            gene_pairs = pairs_df[pairs_df["gene"] == gene]

            # Find canonical sequence
            canonical_seqs = gene_pairs[gene_pairs["variant"] == "canonical"]
            if canonical_seqs.empty:
                continue

            canonical_loc = canonical_seqs.iloc[0]["prediction"]

            # Count truncating isoforms with localization changes (changed from "variants")
            truncated_seqs = gene_pairs[gene_pairs["variant"].str.startswith("trunc")]
            truncating_changes = 0
            for _, trunc_seq in truncated_seqs.iterrows():
                if self.are_localizations_different(
                    canonical_loc, trunc_seq["prediction"]
                ):
                    truncating_changes += 1

            # Count missense isoforms with localization changes (changed from "variants")
            missense_changes = 0
            if not mutations_df.empty:
                gene_mutations = mutations_df[mutations_df["gene"] == gene]
                missense_variants = gene_mutations[
                    gene_mutations["variant"].str.contains("mut", na=False)
                ]

                for _, mut_seq in missense_variants.iterrows():
                    if self.are_localizations_different(
                        canonical_loc, mut_seq["prediction"]
                    ):
                        missense_changes += 1

            # Get mutation counts from pair_results if available
            frameshift_mutations = 0
            nonsense_mutations = 0
            missense_mutations = 0
            total_mutations = 0
            frameshift_variants = []
            nonsense_variants = []
            missense_variants = []

            if pair_results is not None:
                gene_mutation_data = pair_results[pair_results["gene_name"] == gene]
                if not gene_mutation_data.empty:
                    # Sum up mutations for this gene across all transcript-truncation pairs
                    if "mutations_frameshift_variant" in gene_mutation_data.columns:
                        frameshift_mutations = gene_mutation_data[
                            "mutations_frameshift_variant"
                        ].sum()
                    if "mutations_nonsense_variant" in gene_mutation_data.columns:
                        nonsense_mutations = gene_mutation_data[
                            "mutations_nonsense_variant"
                        ].sum()
                    if "mutations_missense_variant" in gene_mutation_data.columns:
                        missense_mutations = gene_mutation_data[
                            "mutations_missense_variant"
                        ].sum()
                    if "mutation_count_total" in gene_mutation_data.columns:
                        total_mutations = gene_mutation_data[
                            "mutation_count_total"
                        ].sum()

                    # Collect variant IDs for each mutation type
                    for _, row in gene_mutation_data.iterrows():
                        # Extract frameshift variant IDs
                        if (
                            "clinvar_ids_frameshift_variant" in row.index
                            and pd.notna(row["clinvar_ids_frameshift_variant"])
                            and row["clinvar_ids_frameshift_variant"] != ""
                        ):
                            frameshift_variants.extend(
                                str(row["clinvar_ids_frameshift_variant"]).split(",")
                            )

                        # Extract nonsense variant IDs
                        if (
                            "clinvar_ids_nonsense_variant" in row.index
                            and pd.notna(row["clinvar_ids_nonsense_variant"])
                            and row["clinvar_ids_nonsense_variant"] != ""
                        ):
                            nonsense_variants.extend(
                                str(row["clinvar_ids_nonsense_variant"]).split(",")
                            )

                        # Extract missense variant IDs
                        if (
                            "clinvar_ids_missense_variant" in row.index
                            and pd.notna(row["clinvar_ids_missense_variant"])
                            and row["clinvar_ids_missense_variant"] != ""
                        ):
                            missense_variants.extend(
                                str(row["clinvar_ids_missense_variant"]).split(",")
                            )

                    # Clean up variant lists (remove empty strings and duplicates)
                    frameshift_variants = list(
                        set([v.strip() for v in frameshift_variants if v.strip()])
                    )
                    nonsense_variants = list(
                        set([v.strip() for v in nonsense_variants if v.strip()])
                    )
                    missense_variants = list(
                        set([v.strip() for v in missense_variants if v.strip()])
                    )

            # Calculate frameshift OR nonsense mutations (changed from "truncating")
            frameshift_or_nonsense_mutations = frameshift_mutations + nonsense_mutations
            frameshift_or_nonsense_variants = list(
                set(frameshift_variants + nonsense_variants)
            )

            # Get gene metadata from protein data
            canonical_seq_id = canonical_seqs.iloc[0]["seq_id"]
            gene_info = protein_data.get(canonical_seq_id, {})

            gene_summary = {
                "dataset": dataset,
                "gene": gene,
                "canonical_localization": canonical_loc,
                # Localization analysis results (truncating = isoforms, missense = variants)
                "total_truncating_isoforms": len(truncated_seqs),
                "truncating_isoforms_with_localization_change": truncating_changes,
                "total_missense_variants": len(missense_variants)
                if not mutations_df.empty
                else 0,
                "missense_variants_with_localization_change": missense_changes,
                "total_variants_with_localization_change": truncating_changes
                + missense_changes,
                # Mutation analysis results (from clinical data)
                "total_frameshift_mutations": frameshift_mutations,
                "total_nonsense_mutations": nonsense_mutations,
                "total_frameshift_or_nonsense_mutations": frameshift_or_nonsense_mutations,  # Changed from "truncating"
                "total_missense_mutations": missense_mutations,
                "total_mutations_all_types": total_mutations,
                # Variant IDs for each mutation type
                "frameshift_variant_ids": ",".join(frameshift_variants)
                if frameshift_variants
                else "",
                "nonsense_variant_ids": ",".join(nonsense_variants)
                if nonsense_variants
                else "",
                "frameshift_or_nonsense_variant_ids": ",".join(
                    frameshift_or_nonsense_variants
                )
                if frameshift_or_nonsense_variants
                else "",  # Changed from "truncating"
                "missense_variant_ids": ",".join(missense_variants)
                if missense_variants
                else "",
                "model_type": model_type,
                # Add gene metadata
                "transcript_id": gene_info.get("transcript_id", ""),
                "canonical_sequence_length": gene_info.get("length", ""),
            }

            gene_summaries.append(gene_summary)

        gene_summary_df = pd.DataFrame(gene_summaries)

        # Sort by total variants with localization changes (prioritized targets)
        if not gene_summary_df.empty:
            gene_summary_df = gene_summary_df.sort_values(
                "total_variants_with_localization_change", ascending=False
            )

        return gene_summary_df

    def analyze_dataset(self, dataset):
        """Analyze a complete dataset and save all results."""
        print(f"\nAnalyzing {dataset} dataset...")

        # Create summary directory for this dataset
        summary_dir = Path(f"../results/{dataset}/summary")
        summary_dir.mkdir(parents=True, exist_ok=True)

        # Load all data for this dataset
        print(f"\nLoading data for {dataset} dataset...")
        gene_results, pair_results = self.load_mutation_results(dataset)
        protein_data = self.load_protein_sequences(dataset)
        loc_results = self.load_localization_results(dataset)

        # Analyze mutations
        mutation_summary = self.analyze_mutations(dataset, gene_results, pair_results)

        # Analyze localizations
        localization_summary, localization_comparisons = self.analyze_localizations(
            dataset, loc_results, protein_data
        )

        # Create detailed localization analysis
        detailed_localization_df = self.create_detailed_localization_analysis(
            dataset, loc_results, protein_data
        )

        # Create gene-level summary
        gene_summary_df = self.create_gene_summary(
            dataset, loc_results, protein_data, pair_results
        )

        # Save all results
        print(f"\n=== SAVING RESULTS FOR {dataset.upper()} DATASET ===")

        # Save mutation summary
        with open(summary_dir / "mutation_summary.txt", "w") as f:
            f.write("\n".join(mutation_summary))
        print("Saved mutation summary")

        # Save localization summary
        with open(summary_dir / "localization_summary.txt", "w") as f:
            f.write("\n".join(localization_summary))
        print("Saved localization summary")

        # Save genes with localization changes (only variants with changes)
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

        # Save detailed localization analysis (all variants assessed)
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

        # Save gene-level summary (prioritized targets)
        if not gene_summary_df.empty:
            gene_summary_df.to_csv(summary_dir / "gene_level_summary.csv", index=False)
            print(f"Saved gene-level summary for {len(gene_summary_df)} genes")
        else:
            # Create empty file
            pd.DataFrame().to_csv(summary_dir / "gene_level_summary.csv", index=False)
            print("No gene-level summary data available")

        print(f"\nSummary analysis for {dataset} dataset completed!")
        print(f"Results saved to: {summary_dir}")

        # Print brief summary of gene-level results
        if not gene_summary_df.empty:
            print(f"\n=== GENE-LEVEL SUMMARY FOR {dataset.upper()} DATASET ===")

            # Top genes by total variants with localization changes
            top_genes = gene_summary_df.head(10)
            print(f"Top 10 genes by variants with localization changes:")
            for _, gene in top_genes.iterrows():
                print(
                    f"  {gene['gene']}: {gene['total_variants_with_localization_change']} variants with loc changes "
                    f"({gene['truncating_isoforms_with_localization_change']} truncating, "
                    f"{gene['missense_variants_with_localization_change']} missense) | "
                    f"Clinical mutations: {gene['total_mutations_all_types']} total "
                    f"({gene['total_frameshift_or_nonsense_mutations']} frameshift/nonsense, {gene['total_missense_mutations']} missense)"
                )

            # Summary statistics
            total_genes_with_changes = len(
                gene_summary_df[
                    gene_summary_df["total_variants_with_localization_change"] > 0
                ]
            )
            total_genes = len(gene_summary_df)
            total_genes_with_mutations = len(
                gene_summary_df[gene_summary_df["total_mutations_all_types"] > 0]
            )
            print(
                f"\nSummary: {total_genes_with_changes} out of {total_genes} genes have variants with localization changes"
            )
            print(
                f"Summary: {total_genes_with_mutations} out of {total_genes} genes have clinical mutations in dataset"
            )
