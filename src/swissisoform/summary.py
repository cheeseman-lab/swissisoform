#!/usr/bin/env python3
"""Summary analysis backend for SwissIsoform pipeline results."""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from collections import defaultdict
import warnings

# Suppress pandas warnings
warnings.filterwarnings("ignore", category=FutureWarning)

logger = logging.getLogger(__name__)


class SummaryAnalyzer:
    """Class to handle summary analysis of SwissIsoform pipeline results."""

    def __init__(self):
        """Initialize the analyzer."""
        pass

    def detect_mutation_sources(self, pair_results):
        """Detect available mutation sources from column names.

        Scans column names to find available sources by looking for patterns:
        - count_{source} (e.g., count_custom, count_clinvar)
        - ids_{source}_* (e.g., ids_custom_missense_variant)
        - {source}_allele_count_* (gnomad specific)
        - {source}_sample_count_* (cosmic specific)

        Args:
            pair_results (pd.DataFrame): Isoform-level mutation results DataFrame.

        Returns:
            List[str]: List of detected source names (e.g., ['custom', 'clinvar']).
        """
        if pair_results is None or pair_results.empty:
            return []

        known_sources = ["custom", "clinvar", "gnomad", "cosmic"]
        detected_sources = []

        columns = pair_results.columns.tolist()

        for source in known_sources:
            # Check for count_{source} column
            if f"count_{source}" in columns:
                detected_sources.append(source)
                continue

            # Check for ids_{source}_* pattern
            if any(col.startswith(f"ids_{source}_") for col in columns):
                detected_sources.append(source)
                continue

            # Check for gnomad-specific pattern: gnomad_allele_count_*
            if source == "gnomad" and any(
                col.startswith("gnomad_allele_count_") for col in columns
            ):
                detected_sources.append(source)
                continue

            # Check for cosmic-specific pattern: cosmic_sample_count_*
            if source == "cosmic" and any(
                col.startswith("cosmic_sample_count_") for col in columns
            ):
                detected_sources.append(source)
                continue

        return detected_sources

    def dataset_has_data(self, dataset, source="gnomad"):
        """Check if a dataset has any data available for analysis.

        Args:
            dataset (str): Name of the dataset to check.
            source (str): Mutation source (default: gnomad).

        Returns:
            bool: True if the dataset has mutation or localization data available.
        """
        # Check for mutation data (source-specific)
        mutation_isoform_file = Path(
            f"../results/{dataset}/{source}/mutations/isoform_level_results.csv"
        )

        # Check for localization data (source-specific mutations)
        loc_files = [
            Path(
                f"../results/{dataset}/{source}/localization/protein_sequences_mutations_Accurate_results.csv"
            ),
            Path(
                f"../results/{dataset}/{source}/localization/protein_sequences_mutations_Fast_results.csv"
            ),
        ]

        has_mutation_data = mutation_isoform_file.exists()
        has_localization_data = any(f.exists() for f in loc_files)

        return has_mutation_data or has_localization_data

    def get_available_models(self, dataset, source="gnomad"):
        """Determine which models have data available for a dataset and source.

        Args:
            dataset (str): Name of the dataset to check.
            source (str): Mutation source (default: gnomad).

        Returns:
            List[str]: List of available model names (e.g., ['Accurate', 'Fast']).
        """
        loc_results = self.load_localization_results(dataset, source)

        available_models = []

        # Check for Accurate model
        if "pairs_accurate" in loc_results:
            available_models.append("Accurate")

        # Check for Fast model
        if "pairs_fast" in loc_results:
            available_models.append("Fast")

        return available_models

    def load_mutation_results(self, dataset, source="gnomad"):
        """Load mutation analysis results for a dataset and source.

        Args:
            dataset (str): Name of the dataset to load results for.
            source (str): Mutation source (default: gnomad).

        Returns:
            Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]: Tuple containing
                gene-level results DataFrame and pair-level results DataFrame.
        """
        base_path = Path(f"../results/{dataset}/{source}/mutations")

        gene_results = None
        pair_results = None

        gene_file = base_path / "gene_level_results.csv"
        pair_file = base_path / "isoform_level_results.csv"

        if gene_file.exists():
            gene_results = pd.read_csv(gene_file)
            logger.info(
                f"Loaded {len(gene_results)} gene-level mutation results for {dataset}/{source}"
            )

        if pair_file.exists():
            pair_results = pd.read_csv(pair_file)
            logger.info(
                f"Loaded {len(pair_results)} transcript-truncation pair results for {dataset}/{source}"
            )

        return gene_results, pair_results

    def load_protein_sequences(self, dataset, source="gnomad"):
        """Load protein sequence data to get variant metadata.

        Args:
            dataset (str): Name of the dataset to load protein sequences for.
            source (str): Mutation source (default: gnomad).

        Returns:
            Dict[str, Dict]: Dictionary mapping sequence IDs to protein metadata.
        """
        # Base pairs are in default folder
        base_pairs_path = Path(f"../results/{dataset}/default/proteins")
        # Mutations are in source-specific folder
        source_path = Path(f"../results/{dataset}/{source}/proteins")

        pairs_file = base_pairs_path / "protein_sequences_pairs.csv"
        mutations_file = source_path / "protein_sequences_with_mutations.csv"

        protein_data = {}

        if pairs_file.exists():
            pairs_df = pd.read_csv(pairs_file)
            logger.info(f"Loaded {len(pairs_df)} protein sequence pairs for {dataset}")
            # Index by sequence identifier
            for _, row in pairs_df.iterrows():
                seq_id = (
                    f"{row['gene_name']}_{row['transcript_id']}_{row['variant_id']}"
                )
                protein_data[seq_id] = row.to_dict()

        if mutations_file.exists():
            mutations_df = pd.read_csv(mutations_file)
            logger.info(
                f"Loaded {len(mutations_df)} protein sequence mutations for {dataset}/{source}"
            )
            # Index by sequence identifier
            for _, row in mutations_df.iterrows():
                seq_id = (
                    f"{row['gene_name']}_{row['transcript_id']}_{row['variant_id']}"
                )
                protein_data[seq_id] = row.to_dict()

        return protein_data

    def load_localization_results(self, dataset, source="gnomad"):
        """Load localization prediction results for a dataset and source.

        Args:
            dataset (str): Name of the dataset to load localization results for.
            source (str): Mutation source (default: gnomad).

        Returns:
            Dict[str, pd.DataFrame]: Dictionary mapping result types to DataFrames
                (e.g., 'pairs_accurate', 'pairs_fast', 'mutations_accurate', 'mutations_fast').
        """
        # Base pairs are in default folder
        base_pairs_path = Path(f"../results/{dataset}/default/localization")
        # Mutations are in source-specific folder
        source_path = Path(f"../results/{dataset}/{source}/localization")

        results = {}

        # Define the files to look for (pairs from default, mutations from source)
        file_patterns = [
            (
                "pairs_accurate",
                base_pairs_path / "protein_sequences_pairs_Accurate_results.csv",
            ),
            (
                "pairs_fast",
                base_pairs_path / "protein_sequences_pairs_Fast_results.csv",
            ),
            (
                "mutations_accurate",
                source_path / "protein_sequences_mutations_Accurate_results.csv",
            ),
            (
                "mutations_fast",
                source_path / "protein_sequences_mutations_Fast_results.csv",
            ),
        ]

        for key, file_path in file_patterns:
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
                    logger.info(
                        f"Loaded {len(df)} {key} localization predictions for {dataset}/{source}"
                    )

                except Exception as e:
                    logger.error(f"Error loading {file_path.name}: {e}")
                    continue
            else:
                logger.warning(
                    f"Missing {key} localization predictions for {dataset}/{source}"
                )

        return results

    def analyze_mutations(self, dataset, gene_results, pair_results):
        """Analyze mutation results for a specific dataset.

        Args:
            dataset (str): Name of the dataset being analyzed.
            gene_results (Optional[pd.DataFrame]): Gene-level mutation results DataFrame.
            pair_results (Optional[pd.DataFrame]): Isoform-level results DataFrame.

        Returns:
            List[str]: List of summary lines describing the mutation analysis.
        """
        logger.info(f"\n=== ANALYZING MUTATIONS FOR {dataset.upper()} DATASET ===")

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
            # Detect available sources
            detected_sources = self.detect_mutation_sources(pair_results)
            if detected_sources:
                summary_lines.append(
                    f"\nMutation sources detected: {', '.join(detected_sources)}"
                )

            # Look for total count columns (count_{variant_type})
            variant_types = [
                "missense_variant",
                "nonsense_variant",
                "frameshift_variant",
                "inframe_deletion",
                "inframe_insertion",
                "synonymous_variant",
            ]

            count_columns = [
                f"count_{vt}"
                for vt in variant_types
                if f"count_{vt}" in pair_results.columns
            ]

            if count_columns:
                summary_lines.append(f"\nMutation breakdown by type (all sources):")

                for col in count_columns:
                    # Extract mutation type from column name
                    mutation_type = col.replace("count_", "").replace("_", " ").title()

                    # Count genes that have at least one mutation of this type
                    genes_with_mutation = pair_results[pair_results[col] > 0][
                        "gene_name"
                    ].nunique()
                    total_mutations = int(pair_results[col].sum())

                    summary_lines.append(
                        f"  {mutation_type}: {genes_with_mutation} genes ({total_mutations} total mutations)"
                    )

                # Special combined category: frameshift OR nonsense
                frameshift_col = "count_frameshift_variant"
                nonsense_col = "count_nonsense_variant"

                if (
                    frameshift_col in pair_results.columns
                    and nonsense_col in pair_results.columns
                ):
                    genes_with_frameshift_or_nonsense = pair_results[
                        (pair_results[frameshift_col] > 0)
                        | (pair_results[nonsense_col] > 0)
                    ]["gene_name"].nunique()

                    total_frameshift_or_nonsense_mutations = int(
                        pair_results[frameshift_col].sum()
                        + pair_results[nonsense_col].sum()
                    )

                    summary_lines.append(
                        f"  Frameshift OR Nonsense: {genes_with_frameshift_or_nonsense} genes ({total_frameshift_or_nonsense_mutations} total mutations)"
                    )

                # Total mutations
                total_mutations_all = int(
                    pair_results["total_mutations"].sum()
                    if "total_mutations" in pair_results.columns
                    else 0
                )
                summary_lines.append(
                    f"\nTotal mutations across all isoform pairs: {total_mutations_all}"
                )

                # Per-source breakdown if multiple sources
                if len(detected_sources) > 1:
                    summary_lines.append(f"\nPer-source breakdown:")
                    for source in detected_sources:
                        source_count_col = f"count_{source}"
                        if source_count_col in pair_results.columns:
                            source_total = int(pair_results[source_count_col].sum())
                            genes_with_source = pair_results[
                                pair_results[source_count_col] > 0
                            ]["gene_name"].nunique()
                            summary_lines.append(
                                f"  {source.upper()}: {genes_with_source} genes ({source_total} total mutations)"
                            )
                elif len(detected_sources) == 1:
                    source = detected_sources[0]
                    summary_lines.append(f"\nSource: {source.upper()}")
            else:
                summary_lines.append(
                    "No mutation count columns found in isoform results"
                )
        else:
            summary_lines.append("No isoform-level mutation data available")

        return summary_lines

    def parse_sequence_id(self, seq_id):
        """Parse a sequence ID to extract gene, transcript, and variant information.

        Args:
            seq_id (str): Sequence identifier to parse.

        Returns:
            Tuple[Optional[str], Optional[str], Optional[str]]: Tuple containing gene name,
                transcript ID, and variant identifier.
        """
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
        """Normalize localization by sorting multiple locations.

        Args:
            localization (str): Localization string (may contain multiple locations separated by '|').

        Returns:
            str: Normalized localization string with sorted location names.
        """
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
        """Check if two localizations are different, accounting for multi-location predictions.

        Args:
            loc1 (str): First localization string.
            loc2 (str): Second localization string.

        Returns:
            bool: True if the normalizations are different.
        """
        if pd.isna(loc1) or pd.isna(loc2):
            return True

        # Normalize both localizations
        norm_loc1 = self.normalize_localization(loc1)
        norm_loc2 = self.normalize_localization(loc2)

        return norm_loc1 != norm_loc2

    def analyze_localizations_for_model(
        self, dataset, loc_results, protein_data, model_type
    ):
        """Analyze localization predictions for a specific dataset and model."""
        logger.info(
            f"\n=== ANALYZING LOCALIZATIONS FOR {dataset.upper()} DATASET - {model_type.upper()} MODEL ==="
        )

        summary_lines = []
        summary_lines.append(
            f"LOCALIZATION ANALYSIS SUMMARY - {dataset.upper()} DATASET - {model_type.upper()} MODEL"
        )
        summary_lines.append("=" * 60)

        all_localization_comparisons = []

        if not loc_results:
            summary_lines.append("No localization data available")
            return summary_lines, all_localization_comparisons

        # Get data for the specified model
        pairs_data = None
        mutations_data = None

        if model_type.lower() == "accurate":
            pairs_data = loc_results.get("pairs_accurate")
            mutations_data = loc_results.get("mutations_accurate")
        elif model_type.lower() == "fast":
            pairs_data = loc_results.get("pairs_fast")
            mutations_data = loc_results.get("mutations_fast")

        if pairs_data is None:
            summary_lines.append(f"No {model_type} model localization data available")
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

                # Find canonical and non-canonical sequences (truncations + extensions)
                canonical_seqs = gene_data[gene_data["variant"] == "canonical"]
                truncated_seqs = gene_data[gene_data["variant"] != "canonical"]

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

    def create_detailed_localization_analysis_for_model(
        self, dataset, loc_results, protein_data, model_type
    ):
        """Create a detailed analysis of all localization predictions for a specific dataset and model."""
        logger.info(
            f"\n=== CREATING DETAILED LOCALIZATION ANALYSIS FOR {dataset.upper()} DATASET - {model_type.upper()} MODEL ==="
        )

        detailed_results = []

        if not loc_results:
            return pd.DataFrame(detailed_results)

        # Get data for the specified model
        pairs_data = None
        mutations_data = None

        if model_type.lower() == "accurate":
            pairs_data = loc_results.get("pairs_accurate")
            mutations_data = loc_results.get("mutations_accurate")
        elif model_type.lower() == "fast":
            pairs_data = loc_results.get("pairs_fast")
            mutations_data = loc_results.get("mutations_fast")

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

    def create_gene_summary_for_model(
        self, dataset, loc_results, protein_data, model_type, pair_results=None
    ):
        """Create a gene-level summary showing prioritized targets for a specific model."""
        logger.info(
            f"\n=== CREATING GENE-LEVEL SUMMARY FOR {dataset.upper()} DATASET - {model_type.upper()} MODEL ==="
        )

        if not loc_results:
            return pd.DataFrame()

        # Get data for the specified model
        pairs_data = None
        mutations_data = None

        if model_type.lower() == "accurate":
            pairs_data = loc_results.get("pairs_accurate")
            mutations_data = loc_results.get("mutations_accurate")
        elif model_type.lower() == "fast":
            pairs_data = loc_results.get("pairs_fast")
            mutations_data = loc_results.get("mutations_fast")

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

            # Count non-canonical isoforms with localization changes (truncations + extensions)
            truncated_seqs = gene_pairs[gene_pairs["variant"] != "canonical"]
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
            source_counts = {}  # Track per-source mutation counts

            if pair_results is not None:
                gene_mutation_data = pair_results[pair_results["gene_name"] == gene]
                if not gene_mutation_data.empty:
                    # Sum up mutations for this gene using count_{variant_type} columns
                    if "count_frameshift_variant" in gene_mutation_data.columns:
                        frameshift_mutations = int(
                            gene_mutation_data["count_frameshift_variant"].sum()
                        )
                    if "count_nonsense_variant" in gene_mutation_data.columns:
                        nonsense_mutations = int(
                            gene_mutation_data["count_nonsense_variant"].sum()
                        )
                    if "count_missense_variant" in gene_mutation_data.columns:
                        missense_mutations = int(
                            gene_mutation_data["count_missense_variant"].sum()
                        )
                    if "total_mutations" in gene_mutation_data.columns:
                        total_mutations = int(
                            gene_mutation_data["total_mutations"].sum()
                        )

                    # Detect available sources and collect variant IDs
                    detected_sources = self.detect_mutation_sources(pair_results)

                    for source in detected_sources:
                        source_counts[source] = {
                            "total": 0,
                            "frameshift": 0,
                            "nonsense": 0,
                            "missense": 0,
                        }

                        # Get source total
                        source_count_col = f"count_{source}"
                        if source_count_col in gene_mutation_data.columns:
                            source_counts[source]["total"] = int(
                                gene_mutation_data[source_count_col].sum()
                            )

                        # Get per-variant-type counts for this source
                        for vt, vt_key in [
                            ("frameshift_variant", "frameshift"),
                            ("nonsense_variant", "nonsense"),
                            ("missense_variant", "missense"),
                        ]:
                            col = f"count_{source}_{vt}"
                            if col in gene_mutation_data.columns:
                                source_counts[source][vt_key] = int(
                                    gene_mutation_data[col].sum()
                                )

                    # Collect variant IDs from all detected sources
                    for _, row in gene_mutation_data.iterrows():
                        for source in detected_sources:
                            # Extract frameshift variant IDs
                            col = f"ids_{source}_frameshift_variant"
                            if (
                                col in row.index
                                and pd.notna(row[col])
                                and row[col] != ""
                            ):
                                frameshift_variants.extend(str(row[col]).split(","))

                            # Extract nonsense variant IDs
                            col = f"ids_{source}_nonsense_variant"
                            if (
                                col in row.index
                                and pd.notna(row[col])
                                and row[col] != ""
                            ):
                                nonsense_variants.extend(str(row[col]).split(","))

                            # Extract missense variant IDs
                            col = f"ids_{source}_missense_variant"
                            if (
                                col in row.index
                                and pd.notna(row[col])
                                and row[col] != ""
                            ):
                                missense_variants.extend(str(row[col]).split(","))

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
                "total_frameshift_or_nonsense_mutations": frameshift_or_nonsense_mutations,
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
                else "",
                "missense_variant_ids": ",".join(missense_variants)
                if missense_variants
                else "",
                "model_type": model_type,
                # Add gene metadata
                "transcript_id": gene_info.get("transcript_id", ""),
                "canonical_sequence_length": gene_info.get("length", ""),
            }

            # Add per-source mutation counts dynamically
            for source, counts in source_counts.items():
                gene_summary[f"{source}_total_mutations"] = counts["total"]
                gene_summary[f"{source}_frameshift_mutations"] = counts["frameshift"]
                gene_summary[f"{source}_nonsense_mutations"] = counts["nonsense"]
                gene_summary[f"{source}_missense_mutations"] = counts["missense"]

            gene_summaries.append(gene_summary)

        gene_summary_df = pd.DataFrame(gene_summaries)

        # Sort by total variants with localization changes (prioritized targets)
        if not gene_summary_df.empty:
            gene_summary_df = gene_summary_df.sort_values(
                "total_variants_with_localization_change", ascending=False
            )

        return gene_summary_df

    def analyze_dataset(self, dataset, source="gnomad"):
        """Analyze a complete dataset and source, saving all results for both models when available.

        Args:
            dataset (str): Name of the dataset to analyze.
            source (str): Mutation source (default: gnomad).
        """
        logger.info(f"\nAnalyzing {dataset}/{source}...")

        # Create summary directory for this dataset/source
        summary_dir = Path(f"../results/{dataset}/{source}/summary")
        summary_dir.mkdir(parents=True, exist_ok=True)

        # Load shared data for this dataset/source
        logger.info(f"\nLoading data for {dataset}/{source}...")
        gene_results, pair_results = self.load_mutation_results(dataset, source)
        protein_data = self.load_protein_sequences(dataset, source)
        loc_results = self.load_localization_results(dataset, source)

        # Analyze mutations (same for both models)
        mutation_summary = self.analyze_mutations(dataset, gene_results, pair_results)

        # Get available models
        available_models = self.get_available_models(dataset, source)

        if not available_models:
            logger.warning(f"No localization models available for {dataset}/{source}")
            return

        logger.info(
            f"Available models for {dataset}/{source}: {', '.join(available_models)}"
        )

        # Analyze each available model separately
        for model_type in available_models:
            logger.info(
                f"\n=== ANALYZING {model_type.upper()} MODEL FOR {dataset.upper()} DATASET ==="
            )

            # Create model-specific subdirectory
            model_summary_dir = summary_dir / model_type.lower()
            model_summary_dir.mkdir(parents=True, exist_ok=True)

            # Analyze localizations for this model
            localization_summary, localization_comparisons = (
                self.analyze_localizations_for_model(
                    dataset, loc_results, protein_data, model_type
                )
            )

            # Create detailed localization analysis for this model
            detailed_localization_df = (
                self.create_detailed_localization_analysis_for_model(
                    dataset, loc_results, protein_data, model_type
                )
            )

            # Create gene-level summary for this model
            gene_summary_df = self.create_gene_summary_for_model(
                dataset, loc_results, protein_data, model_type, pair_results
            )

            # Save model-specific results
            logger.info(
                f"\n=== SAVING {model_type.upper()} MODEL RESULTS FOR {dataset.upper()} DATASET ==="
            )

            # Save localization summary
            with open(model_summary_dir / "localization_summary.txt", "w") as f:
                f.write("\n".join(localization_summary))
            logger.info(f"Saved {model_type} localization summary")

            # Save genes with localization changes (only variants with changes)
            if localization_comparisons:
                localization_changes_df = pd.DataFrame(localization_comparisons)
                localization_changes_df.to_csv(
                    model_summary_dir / "genes_with_localization_changes.csv",
                    index=False,
                )
                logger.info(
                    f"Saved {len(localization_changes_df)} {model_type} localization changes"
                )
            else:
                # Create empty file
                pd.DataFrame().to_csv(
                    model_summary_dir / "genes_with_localization_changes.csv",
                    index=False,
                )
                logger.info(f"No {model_type} localization changes found")

            # Save detailed localization analysis (all variants assessed)
            if not detailed_localization_df.empty:
                detailed_localization_df.to_csv(
                    model_summary_dir / "detailed_localization_analysis.csv",
                    index=False,
                )
                logger.info(
                    f"Saved detailed {model_type} analysis of {len(detailed_localization_df)} localization predictions"
                )
            else:
                # Create empty file
                pd.DataFrame().to_csv(
                    model_summary_dir / "detailed_localization_analysis.csv",
                    index=False,
                )
                logger.info(f"No detailed {model_type} localization data available")

            # Save gene-level summary (prioritized targets)
            if not gene_summary_df.empty:
                gene_summary_df.to_csv(
                    model_summary_dir / "gene_level_summary.csv", index=False
                )
                logger.info(
                    f"Saved {model_type} gene-level summary for {len(gene_summary_df)} genes"
                )
            else:
                # Create empty file
                pd.DataFrame().to_csv(
                    model_summary_dir / "gene_level_summary.csv", index=False
                )
                logger.info(f"No {model_type} gene-level summary data available")

            # Print brief summary of gene-level results for this model
            if not gene_summary_df.empty:
                logger.info(
                    f"\n=== GENE-LEVEL SUMMARY FOR {dataset.upper()} DATASET - {model_type.upper()} MODEL ==="
                )

                # Top genes by total variants with localization changes
                top_genes = gene_summary_df.head(10)
                logger.info(
                    f"Top 10 genes by variants with localization changes ({model_type} model):"
                )
                for _, gene in top_genes.iterrows():
                    logger.info(
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
                logger.info(
                    f"\n{model_type} model summary: {total_genes_with_changes} out of {total_genes} genes have variants with localization changes"
                )
                logger.info(
                    f"{model_type} model summary: {total_genes_with_mutations} out of {total_genes} genes have clinical mutations in dataset"
                )

        # Save mutation summary to main summary directory (shared across models)
        with open(summary_dir / "mutation_summary.txt", "w") as f:
            f.write("\n".join(mutation_summary))
        logger.info("Saved mutation summary (shared)")

        logger.info(f"\nSummary analysis for {dataset} dataset completed!")
        logger.info(f"Results saved to: {summary_dir}")
        for model_type in available_models:
            logger.info(
                f"  {model_type} model results: {summary_dir / model_type.lower()}/"
            )

    # Keep original methods for backward compatibility (they now use the first available model)
    def analyze_localizations(self, dataset, loc_results, protein_data):
        """Backward compatibility wrapper - uses first available model."""
        available_models = []
        if "pairs_accurate" in loc_results:
            available_models.append("Accurate")
        if "pairs_fast" in loc_results:
            available_models.append("Fast")

        if not available_models:
            return [], []

        # Use first available model for backward compatibility
        return self.analyze_localizations_for_model(
            dataset, loc_results, protein_data, available_models[0]
        )

    def create_detailed_localization_analysis(self, dataset, loc_results, protein_data):
        """Backward compatibility wrapper - uses first available model."""
        available_models = []
        if "pairs_accurate" in loc_results:
            available_models.append("Accurate")
        if "pairs_fast" in loc_results:
            available_models.append("Fast")

        if not available_models:
            return pd.DataFrame()

        # Use first available model for backward compatibility
        return self.create_detailed_localization_analysis_for_model(
            dataset, loc_results, protein_data, available_models[0]
        )

    def create_gene_summary(
        self, dataset, loc_results, protein_data, pair_results=None
    ):
        """Backward compatibility wrapper - uses first available model."""
        available_models = []
        if "pairs_accurate" in loc_results:
            available_models.append("Accurate")
        if "pairs_fast" in loc_results:
            available_models.append("Fast")

        if not available_models:
            return pd.DataFrame()

        # Use first available model for backward compatibility
        return self.create_gene_summary_for_model(
            dataset, loc_results, protein_data, available_models[0], pair_results
        )
