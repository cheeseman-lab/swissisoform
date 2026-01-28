#!/usr/bin/env python3
"""Summary analysis backend for SwissIsoform pipeline results."""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from collections import defaultdict
import warnings
import re

# Suppress pandas warnings
warnings.filterwarnings("ignore", category=FutureWarning)

logger = logging.getLogger(__name__)

# Localization compartment columns from DeepLoc output
COMPARTMENT_COLS = [
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

# Signal/membrane type columns from DeepLoc output
SIGNAL_COLS = [
    "Peripheral",
    "Transmembrane",
    "Lipid anchor",
    "Soluble",
]

# All probability columns
ALL_PROB_COLS = COMPARTMENT_COLS + SIGNAL_COLS


class SummaryAnalyzer:
    """Class to handle summary analysis of SwissIsoform pipeline results."""

    def __init__(self):
        """Initialize the analyzer."""
        self.compartment_cols = COMPARTMENT_COLS
        self.signal_cols = SIGNAL_COLS
        self.all_prob_cols = ALL_PROB_COLS

    # =========================================================================
    # Probability-based Analysis Helper Methods
    # =========================================================================

    def compute_shift_magnitude(self, row1, row2, prob_cols=None):
        """Compute maximum probability change across all compartments.

        Args:
            row1: First row (Series) with probability columns.
            row2: Second row (Series) with probability columns.
            prob_cols: List of probability columns to compare (default: COMPARTMENT_COLS).

        Returns:
            float: Maximum absolute probability shift.
        """
        if prob_cols is None:
            prob_cols = self.compartment_cols

        max_shift = 0.0
        for col in prob_cols:
            if col in row1.index and col in row2.index:
                shift = abs(float(row2[col]) - float(row1[col]))
                if shift > max_shift:
                    max_shift = shift
        return max_shift

    def get_top_localizations(self, row, prob_cols=None, n=3):
        """Return top N localizations with their probabilities.

        Args:
            row: Row (Series) with probability columns.
            prob_cols: List of probability columns to check (default: COMPARTMENT_COLS).
            n: Number of top localizations to return.

        Returns:
            List[Tuple[str, float]]: List of (compartment, probability) tuples.
        """
        if prob_cols is None:
            prob_cols = self.compartment_cols

        probs = {}
        for col in prob_cols:
            if col in row.index:
                val = row[col]
                if pd.notna(val):
                    probs[col] = float(val)

        sorted_probs = sorted(probs.items(), key=lambda x: x[1], reverse=True)
        return sorted_probs[:n]

    def get_primary_confidence(self, row, prob_cols=None):
        """Get confidence of primary prediction (max probability).

        Args:
            row: Row (Series) with probability columns.
            prob_cols: List of probability columns to check (default: COMPARTMENT_COLS).

        Returns:
            float: Maximum probability value (confidence of primary localization).
        """
        if prob_cols is None:
            prob_cols = self.compartment_cols

        max_prob = 0.0
        for col in prob_cols:
            if col in row.index:
                val = row[col]
                if pd.notna(val) and float(val) > max_prob:
                    max_prob = float(val)
        return max_prob

    def get_primary_localization(self, row, prob_cols=None):
        """Get the compartment with highest probability.

        Args:
            row: Row (Series) with probability columns.
            prob_cols: List of probability columns to check (default: COMPARTMENT_COLS).

        Returns:
            str: Name of compartment with highest probability.
        """
        if prob_cols is None:
            prob_cols = self.compartment_cols

        top = self.get_top_localizations(row, prob_cols, n=1)
        return top[0][0] if top else "Unknown"

    def format_top_probs(self, row, prob_cols=None, n=3):
        """Format top N probabilities as a string.

        Args:
            row: Row (Series) with probability columns.
            prob_cols: List of probability columns to check.
            n: Number of top localizations to format.

        Returns:
            str: Formatted string like "Nucleus(0.72)|Cytoplasm(0.18)|ER(0.05)"
        """
        top = self.get_top_localizations(row, prob_cols, n)
        parts = []
        for loc, prob in top:
            # Abbreviate some long names
            abbrev = {
                "Endoplasmic reticulum": "ER",
                "Lysosome/Vacuole": "Lyso/Vac",
                "Cell membrane": "CellMem",
            }
            loc_name = abbrev.get(loc, loc)
            parts.append(f"{loc_name}({prob:.2f})")
        return "|".join(parts)

    def parse_protein_id(self, protein_id):
        """Parse a Protein_ID to extract gene, transcript, feature_type, and mutation info.

        Args:
            protein_id (str): Protein identifier string.

        Returns:
            dict: Dictionary with gene, transcript, feature_type, feature_id, is_mutant,
                  is_canonical, mutation_info keys.
        """
        if pd.isna(protein_id) or not protein_id:
            return {
                "gene": None,
                "transcript": None,
                "feature_type": None,
                "feature_id": None,
                "is_mutant": False,
                "is_canonical": False,
                "mutation_info": None,
            }

        protein_id = str(protein_id).strip()
        parts = protein_id.split("_")
        gene = parts[0]
        transcript = parts[1] if len(parts) > 1 else None

        is_mutant = "_mutated_" in protein_id
        is_canonical = protein_id.endswith("_canonical")

        if is_canonical:
            feature_type = "canonical"
            feature_id = protein_id
        elif "_extension_" in protein_id:
            feature_type = "extension"
            # Extract feature_id (everything after gene_transcript_)
            match = re.search(r"(extension_[A-Z]+_\d+)", protein_id)
            feature_id = match.group(1) if match else protein_id
        elif "_truncation_" in protein_id:
            feature_type = "truncation"
            match = re.search(r"(truncation_[A-Z]+_\d+_\d+)", protein_id)
            feature_id = match.group(1) if match else protein_id
        else:
            feature_type = "unknown"
            feature_id = protein_id

        mutation_info = None
        if is_mutant:
            # Pattern: _mutated_{position}_{change}_{aa_change}
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
            "feature_id": feature_id,
            "is_mutant": is_mutant,
            "is_canonical": is_canonical,
            "mutation_info": mutation_info,
        }

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
            source (str): Mutation source (default: gnomad, or 'default' for structural analysis).

        Returns:
            bool: True if the dataset has mutation or localization data available.
        """
        # For DEFAULT source, check for pairs localization data
        if source == "default":
            pairs_loc_files = [
                Path(
                    f"../results/{dataset}/default/localization/protein_sequences_pairs_Accurate_results.csv"
                ),
                Path(
                    f"../results/{dataset}/default/localization/protein_sequences_pairs_Fast_results.csv"
                ),
            ]
            return any(f.exists() for f in pairs_loc_files)

        # For mutation sources, check for mutation data
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

    def load_mutation_details(self, dataset, source):
        """Load detailed per-gene mutation CSVs from mutations/{gene}/{transcript}/ directories.

        This extracts detailed mutation metadata including:
        - clinvar_clinical_significance, clinvar_star_rating, clinvar_condition
        - gnomad_allele_frequency, gnomad_allele_count, gnomad_allele_count_hom
        - cosmic_sample_count
        - in_alt_start_site, is_alt_start_loss
        - impact, hgvsc, hgvsp

        Args:
            dataset (str): Name of the dataset.
            source (str): Mutation source (clinvar, cosmic, gnomad, custom).

        Returns:
            Dict[Tuple[str, str, str], pd.DataFrame]: Dictionary mapping
                (gene, transcript, feature_id) to DataFrame of mutations.
        """
        base_path = Path(f"../results/{dataset}/{source}/mutations")

        if not base_path.exists():
            logger.warning(f"Mutation directory not found: {base_path}")
            return {}

        mutation_details = {}

        # Iterate through gene directories
        for gene_dir in base_path.iterdir():
            if not gene_dir.is_dir() or gene_dir.name.startswith("."):
                continue

            gene_name = gene_dir.name

            # Iterate through transcript directories
            for transcript_dir in gene_dir.iterdir():
                if not transcript_dir.is_dir():
                    continue

                transcript_id = transcript_dir.name

                # Find mutation CSV files (prefer filtered if available)
                csv_files = list(transcript_dir.glob("*_mutations*.csv"))

                for csv_file in csv_files:
                    try:
                        df = pd.read_csv(csv_file)

                        if df.empty:
                            continue

                        # Extract feature_id from filename
                        # Format: TRANSCRIPT_feature_type_CODON_START_END_mutations.csv
                        filename = csv_file.stem
                        match = re.search(
                            r"(extension_[A-Z]+_\d+|truncation_[A-Z]+_\d+_\d+)",
                            filename,
                        )
                        feature_id = match.group(1) if match else filename

                        key = (gene_name, transcript_id, feature_id)

                        # If we already have data for this key, prefer filtered version
                        if key in mutation_details:
                            if "filtered" in csv_file.name:
                                mutation_details[key] = df
                        else:
                            mutation_details[key] = df

                    except Exception as e:
                        logger.debug(f"Could not load {csv_file}: {e}")
                        continue

        logger.info(
            f"Loaded mutation details for {len(mutation_details)} features from {dataset}/{source}"
        )
        return mutation_details

    def get_source_specific_rankings(self, source):
        """Return list of ranking columns available for each source.

        Args:
            source (str): Mutation source name.

        Returns:
            Dict[str, List[str]]: Dictionary with 'required' and 'optional' ranking columns.
        """
        rankings = {
            "clinvar": {
                "required": ["clinical_significance", "star_rating"],
                "optional": ["condition"],
            },
            "cosmic": {
                "required": ["sample_count"],
                "optional": [],
            },
            "gnomad": {
                "required": ["allele_frequency"],
                "optional": ["allele_count", "allele_count_hom"],
            },
            "custom": {
                "required": [],
                "optional": [],
            },
            "default": {
                "required": [],
                "optional": [],
            },
        }

        # Handle custom sources with prefixes
        source_key = source.lower()
        if source_key.startswith("custom"):
            source_key = "custom"

        return rankings.get(source_key, rankings["custom"])

    def load_localization_results(self, dataset, source="gnomad"):
        """Load localization prediction results for a dataset and source.

        Preserves ALL columns from DeepLoc output including:
        - Protein_ID (renamed to Sequence_ID)
        - Localizations (renamed to Localisation)
        - Signals, Membrane types
        - 10 compartment probabilities (Cytoplasm, Nucleus, etc.)
        - 4 signal/membrane type scores (Peripheral, Transmembrane, etc.)

        Args:
            dataset (str): Name of the dataset to load localization results for.
            source (str): Mutation source (default: gnomad).

        Returns:
            Dict[str, pd.DataFrame]: Dictionary mapping result types to DataFrames
                (e.g., 'pairs_accurate', 'pairs_fast', 'mutations_accurate', 'mutations_fast').
                Each DataFrame contains all DeepLoc columns including probability scores.
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
                    # Read all columns - don't filter any out
                    df = pd.read_csv(file_path)

                    # Standardize DeepLoc column names to our expected format
                    # Keep original names as well for compatibility
                    column_mapping = {
                        "Protein_ID": "Sequence_ID",
                        "Localizations": "Localisation",
                    }

                    df = df.rename(columns=column_mapping)

                    # Ensure probability columns are numeric
                    for col in self.all_prob_cols:
                        if col in df.columns:
                            df[col] = pd.to_numeric(df[col], errors="coerce")

                    results[key] = df
                    logger.info(
                        f"Loaded {len(df)} {key} localization predictions for {dataset}/{source} "
                        f"({len(df.columns)} columns including probabilities)"
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
                # LOCALIZATION CHANGE METRICS:
                # IMPORTANT: These have different meanings for default vs. mutation sources
                #
                # For "default" analysis (genomic context):
                #   - truncating_isoforms_with_localization_change: Alternative TIS in genome
                #   - missense_variants_with_localization_change: Not applicable (no mutations)
                #
                # For mutation sources (clinvar/cosmic/custom):
                #   - truncating_isoforms_with_localization_change: Structural features (NOT follow-up targets)
                #   - missense_variants_with_localization_change: Actual mutations (PRIMARY follow-up targets)
                #
                "total_truncating_isoforms": len(truncated_seqs),
                "truncating_isoforms_with_localization_change": truncating_changes,
                "total_missense_variants": len(missense_variants)
                if not mutations_df.empty
                else 0,
                "missense_variants_with_localization_change": missense_changes,  # ← PRIMARY METRIC for mutation sources
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

    # =========================================================================
    # Multi-View Analysis Methods (NEW)
    # =========================================================================

    def create_feature_analysis(
        self,
        dataset,
        loc_results,
        protein_data,
        model_type,
        source,
        pair_results=None,
        mutation_details=None,
        other_model_loc=None,
    ):
        """Create unified feature-level analysis table combining all metrics.

        This is the single comprehensive table that replaces multiple subdirectory files.
        One row per feature (canonical/extension/truncation) with all computed metrics.

        Args:
            dataset: Dataset name.
            loc_results: Dictionary of localization DataFrames.
            protein_data: Dictionary of protein metadata.
            model_type: "Accurate" or "Fast".
            source: Mutation source (default, clinvar, cosmic, gnomad, custom_*).
            pair_results: Optional isoform-level results DataFrame.
            mutation_details: Optional dict of per-gene mutation DataFrames.
            other_model_loc: Optional loc_results from the other model for agreement.

        Returns:
            pd.DataFrame: Unified feature analysis table.
        """
        logger.info(
            f"Creating unified feature analysis for {dataset}/{source}/{model_type}"
        )

        # Get pairs data for the specified model
        if model_type.lower() == "accurate":
            pairs_data = loc_results.get("pairs_accurate")
            mutations_data = loc_results.get("mutations_accurate")
            other_pairs = other_model_loc.get("pairs_fast") if other_model_loc else None
            other_mutations = (
                other_model_loc.get("mutations_fast") if other_model_loc else None
            )
        else:
            pairs_data = loc_results.get("pairs_fast")
            mutations_data = loc_results.get("mutations_fast")
            other_pairs = (
                other_model_loc.get("pairs_accurate") if other_model_loc else None
            )
            other_mutations = (
                other_model_loc.get("mutations_accurate") if other_model_loc else None
            )

        if pairs_data is None or pairs_data.empty:
            logger.warning(f"No pairs data available for {model_type}")
            return pd.DataFrame()

        # Build lookup for other model (for agreement)
        other_lookup = {}
        if other_pairs is not None:
            for _, row in other_pairs.iterrows():
                seq_id = row.get("Sequence_ID", "")
                if seq_id:
                    other_lookup[seq_id] = row
        if other_mutations is not None:
            for _, row in other_mutations.iterrows():
                seq_id = row.get("Sequence_ID", "")
                if seq_id:
                    other_lookup[seq_id] = row

        # Build mutations lookup for missense analysis
        mutations_lookup = {}
        if mutations_data is not None:
            for _, row in mutations_data.iterrows():
                seq_id = row.get("Sequence_ID", "")
                if seq_id:
                    info = self.parse_protein_id(seq_id)
                    # Key by gene + feature_type + mutation info for grouping
                    gene = info["gene"]
                    feature_type = info["feature_type"]
                    if gene not in mutations_lookup:
                        mutations_lookup[gene] = {}
                    if feature_type not in mutations_lookup[gene]:
                        mutations_lookup[gene][feature_type] = []
                    mutations_lookup[gene][feature_type].append((row, info))

        # Parse all protein IDs and group by gene
        parsed_info = {}
        gene_groups = {}
        for _, row in pairs_data.iterrows():
            seq_id = row.get("Sequence_ID", "")
            if seq_id:
                info = self.parse_protein_id(seq_id)
                parsed_info[seq_id] = info
                gene = info["gene"]
                if gene not in gene_groups:
                    gene_groups[gene] = {"canonical": None, "features": []}
                if info["is_canonical"]:
                    gene_groups[gene]["canonical"] = (row, info, seq_id)
                elif (
                    info["feature_type"] in ("extension", "truncation")
                    and not info["is_mutant"]
                ):
                    gene_groups[gene]["features"].append((row, info, seq_id))

        # Build pair_results lookup for mutation burden info
        pair_lookup = {}
        if pair_results is not None and not pair_results.empty:
            for _, row in pair_results.iterrows():
                gene = row.get("gene_name", "")
                feature_id = row.get("feature_id", "")
                if gene and feature_id:
                    pair_lookup[(gene, feature_id)] = row

        # Build feature analysis records
        records = []

        for gene, group in gene_groups.items():
            canonical_entry = group["canonical"]
            if canonical_entry is None:
                continue

            can_row, can_info, can_seq_id = canonical_entry
            can_loc = self.get_primary_localization(can_row)
            can_confidence = self.get_primary_confidence(can_row)
            can_top_probs = self.format_top_probs(can_row, n=3)

            # Get other model's canonical prediction for agreement
            other_can_row = other_lookup.get(can_seq_id)
            other_can_loc = (
                self.get_primary_localization(other_can_row)
                if other_can_row is not None
                else None
            )
            other_can_conf = (
                self.get_primary_confidence(other_can_row)
                if other_can_row is not None
                else None
            )

            for feat_row, feat_info, feat_seq_id in group["features"]:
                feat_loc = self.get_primary_localization(feat_row)
                feat_confidence = self.get_primary_confidence(feat_row)
                feat_top_probs = self.format_top_probs(feat_row, n=3)
                shift_magnitude = self.compute_shift_magnitude(can_row, feat_row)

                # Get other model's feature prediction
                other_feat_row = other_lookup.get(feat_seq_id)
                other_feat_loc = (
                    self.get_primary_localization(other_feat_row)
                    if other_feat_row is not None
                    else None
                )
                other_feat_conf = (
                    self.get_primary_confidence(other_feat_row)
                    if other_feat_row is not None
                    else None
                )

                # Model agreement for this feature
                model_agreement = None
                if other_feat_loc is not None:
                    model_agreement = feat_loc == other_feat_loc

                # Get feature length and other info from pair_results
                feature_length = 0
                total_mutations = 0
                mutation_density = 0.0
                alt_start_loss_count = 0
                alt_start_loss_variants = ""
                count_in_alt_start = 0
                count_in_canonical_start = 0

                # Find matching pair_results row
                pair_row = None
                for (g, fid), prow in pair_lookup.items():
                    if g == gene and feat_info["feature_id"] in str(fid):
                        pair_row = prow
                        break

                if pair_row is not None:
                    # Helper to safely convert to int (handles NaN)
                    def safe_int(val, default=0):
                        if pd.isna(val):
                            return default
                        try:
                            return int(val)
                        except (ValueError, TypeError):
                            return default

                    feature_length = safe_int(pair_row.get("feature_length_aa", 0))
                    total_mutations = safe_int(pair_row.get("total_mutations", 0))
                    alt_start_loss_count = safe_int(
                        pair_row.get("alternative_start_loss_count", 0)
                    )
                    alt_start_loss_variants = str(
                        pair_row.get("alternative_start_loss_variant_ids", "") or ""
                    )
                    count_in_alt_start = safe_int(
                        pair_row.get("count_in_alt_start_site", 0)
                    )
                    count_in_canonical_start = safe_int(
                        pair_row.get("count_in_canonical_start_site", 0)
                    )
                    if feature_length > 0:
                        mutation_density = round(
                            total_mutations / feature_length * 100, 4
                        )

                # Detect signal changes
                signal_changes = []
                for sig_col in self.signal_cols:
                    if sig_col in can_row.index and sig_col in feat_row.index:
                        can_val = (
                            float(can_row[sig_col]) if pd.notna(can_row[sig_col]) else 0
                        )
                        feat_val = (
                            float(feat_row[sig_col])
                            if pd.notna(feat_row[sig_col])
                            else 0
                        )
                        if abs(feat_val - can_val) > 0.2:
                            signal_changes.append(
                                f"{sig_col}({can_val:.2f}→{feat_val:.2f})"
                            )

                # Missense analysis (for mutation sources)
                missense_count = 0
                missense_loc_change_count = 0
                max_missense_shift = 0.0
                top_missense_variant = ""
                top_missense_hgvsp = ""

                if source != "default" and gene in mutations_lookup:
                    feat_type = feat_info["feature_type"]
                    feat_mutations = mutations_lookup[gene].get(feat_type, [])
                    # Also check for mutations on the "canonical" feature type that are related
                    can_mutations = mutations_lookup[gene].get("canonical", [])

                    for mut_row, mut_info in feat_mutations + can_mutations:
                        if mut_info["is_mutant"]:
                            missense_count += 1
                            mut_loc = self.get_primary_localization(mut_row)
                            mut_shift = self.compute_shift_magnitude(can_row, mut_row)
                            if mut_loc != can_loc:
                                missense_loc_change_count += 1
                            if mut_shift > max_missense_shift:
                                max_missense_shift = mut_shift
                                top_missense_variant = mut_row.get("Sequence_ID", "")
                                if mut_info["mutation_info"]:
                                    top_missense_hgvsp = mut_info["mutation_info"].get(
                                        "aa_change", ""
                                    )

                # Source-specific columns
                clinvar_pathogenic_count = 0
                clinvar_top_star_rating = None
                cosmic_max_sample_count = 0
                gnomad_min_af = None

                # Get source-specific data from mutation_details if available
                if mutation_details and source != "default":
                    detail_key = (
                        gene,
                        feat_info["transcript"],
                        feat_info["feature_id"],
                    )
                    if detail_key in mutation_details:
                        detail_df = mutation_details[detail_key]
                        if not detail_df.empty:
                            if source == "clinvar":
                                if "clinvar_clinical_significance" in detail_df.columns:
                                    # Convert to string first to handle mixed types
                                    clin_sig = detail_df[
                                        "clinvar_clinical_significance"
                                    ].astype(str)
                                    pathogenic_mask = clin_sig.str.contains(
                                        "Pathogenic|Likely pathogenic",
                                        case=False,
                                        na=False,
                                    )
                                    clinvar_pathogenic_count = pathogenic_mask.sum()
                                if "clinvar_star_rating" in detail_df.columns:
                                    stars = pd.to_numeric(
                                        detail_df["clinvar_star_rating"],
                                        errors="coerce",
                                    ).dropna()
                                    if not stars.empty:
                                        clinvar_top_star_rating = int(stars.max())
                            elif source == "cosmic":
                                if "cosmic_sample_count" in detail_df.columns:
                                    counts = pd.to_numeric(
                                        detail_df["cosmic_sample_count"],
                                        errors="coerce",
                                    ).dropna()
                                    if not counts.empty:
                                        cosmic_max_sample_count = int(counts.max())
                            elif source == "gnomad":
                                if "gnomad_allele_frequency" in detail_df.columns:
                                    freqs = pd.to_numeric(
                                        detail_df["gnomad_allele_frequency"],
                                        errors="coerce",
                                    ).dropna()
                                    if not freqs.empty:
                                        gnomad_min_af = float(freqs.min())

                record = {
                    # Identity
                    "gene": gene,
                    "transcript": feat_info["transcript"],
                    "feature_id": feat_info["feature_id"],
                    "feature_type": feat_info["feature_type"],
                    "feature_length_aa": feature_length,
                    # Localization comparison
                    "canonical_loc": can_loc,
                    "canonical_confidence": round(can_confidence, 3),
                    "canonical_probs": can_top_probs,
                    "feature_loc": feat_loc,
                    "feature_confidence": round(feat_confidence, 3),
                    "feature_probs": feat_top_probs,
                    "shift_magnitude": round(shift_magnitude, 4),
                    "locs_differ": can_loc != feat_loc,
                    "signal_change": "|".join(signal_changes) if signal_changes else "",
                    # Mutation burden
                    "total_mutations": total_mutations,
                    "mutation_density": mutation_density,
                    "count_in_alt_start": count_in_alt_start,
                    "count_in_canonical_start": count_in_canonical_start,
                    # Alt-start-loss
                    "alt_start_loss_count": alt_start_loss_count,
                    "alt_start_loss_variants": alt_start_loss_variants
                    if alt_start_loss_variants != "nan"
                    else "",
                    # Missense summary
                    "missense_count": missense_count,
                    "missense_loc_change_count": missense_loc_change_count,
                    "max_missense_shift": round(max_missense_shift, 4),
                    "top_missense_variant": top_missense_variant,
                    "top_missense_hgvsp": top_missense_hgvsp,
                    # Model agreement
                    "model_agreement": model_agreement,
                }

                # Add source-specific columns based on source
                if source == "clinvar":
                    record["clinvar_pathogenic_count"] = clinvar_pathogenic_count
                    record["clinvar_top_star_rating"] = clinvar_top_star_rating
                elif source == "cosmic":
                    record["cosmic_max_sample_count"] = cosmic_max_sample_count
                elif source == "gnomad":
                    record["gnomad_min_af"] = gnomad_min_af

                # Add other model columns if available
                if other_feat_loc is not None:
                    other_model = (
                        "fast" if model_type.lower() == "accurate" else "accurate"
                    )
                    record[f"{other_model}_loc"] = other_feat_loc
                    record[f"{other_model}_confidence"] = (
                        round(other_feat_conf, 3) if other_feat_conf else None
                    )

                records.append(record)

        if not records:
            logger.warning("No feature analysis records created")
            return pd.DataFrame()

        df = pd.DataFrame(records)

        # Sort by shift_magnitude descending as default
        df = df.sort_values("shift_magnitude", ascending=False).reset_index(drop=True)

        logger.info(f"Created feature analysis with {len(df)} features")
        return df

    def create_structural_analysis_views(
        self, dataset, loc_results, protein_data, model_type, pair_results=None
    ):
        """Create multiple ranked views of structural changes (DEFAULT analysis only).

        Compares canonical vs extension/truncation isoforms without mutations.
        Generates views ranked by different criteria.

        Args:
            dataset: Dataset name.
            loc_results: Dictionary of localization DataFrames.
            protein_data: Dictionary of protein metadata.
            model_type: "Accurate" or "Fast".
            pair_results: Optional isoform-level results DataFrame.

        Returns:
            Dict[str, pd.DataFrame]: Dictionary mapping view names to DataFrames.
        """
        logger.info(f"Creating structural analysis views for {dataset}/{model_type}")

        views = {}

        # Get pairs data for the specified model
        if model_type.lower() == "accurate":
            pairs_data = loc_results.get("pairs_accurate")
        else:
            pairs_data = loc_results.get("pairs_fast")

        if pairs_data is None or pairs_data.empty:
            logger.warning(f"No pairs data available for {model_type}")
            return views

        # Build comparison records
        records = []

        # Parse all protein IDs
        parsed_info = {}
        for _, row in pairs_data.iterrows():
            seq_id = row.get("Sequence_ID", "")
            if seq_id:
                parsed_info[seq_id] = self.parse_protein_id(seq_id)

        # Group by gene
        gene_groups = {}
        for _, row in pairs_data.iterrows():
            seq_id = row.get("Sequence_ID", "")
            if seq_id and seq_id in parsed_info:
                info = parsed_info[seq_id]
                gene = info["gene"]
                if gene not in gene_groups:
                    gene_groups[gene] = []
                gene_groups[gene].append((row, info))

        # Compare canonical vs alternatives for each gene
        for gene, entries in gene_groups.items():
            canonical_entry = None
            alternative_entries = []

            for row, info in entries:
                if info["is_canonical"]:
                    canonical_entry = (row, info)
                elif (
                    info["feature_type"] in ("extension", "truncation")
                    and not info["is_mutant"]
                ):
                    alternative_entries.append((row, info))

            if canonical_entry is None or not alternative_entries:
                continue

            can_row, can_info = canonical_entry
            can_loc = can_row.get("Localisation", "")
            can_confidence = self.get_primary_confidence(can_row)
            can_primary = self.get_primary_localization(can_row)
            can_top_probs = self.format_top_probs(can_row, n=3)

            for alt_row, alt_info in alternative_entries:
                alt_loc = alt_row.get("Localisation", "")
                alt_confidence = self.get_primary_confidence(alt_row)
                alt_primary = self.get_primary_localization(alt_row)
                alt_top_probs = self.format_top_probs(alt_row, n=3)

                shift_magnitude = self.compute_shift_magnitude(can_row, alt_row)

                # Get length difference from pair_results if available
                length_diff = 0
                if pair_results is not None and not pair_results.empty:
                    feature_rows = pair_results[
                        (pair_results["gene_name"] == gene)
                        & (
                            pair_results["feature_id"].str.contains(
                                alt_info["feature_id"], na=False
                            )
                        )
                    ]
                    if not feature_rows.empty:
                        length_diff = feature_rows.iloc[0].get(
                            "aa_difference_from_canonical", 0
                        )
                        if pd.isna(length_diff):
                            length_diff = 0

                # Detect signal changes
                signal_changes = []
                for sig_col in self.signal_cols:
                    if sig_col in can_row.index and sig_col in alt_row.index:
                        can_val = (
                            float(can_row[sig_col]) if pd.notna(can_row[sig_col]) else 0
                        )
                        alt_val = (
                            float(alt_row[sig_col]) if pd.notna(alt_row[sig_col]) else 0
                        )
                        if abs(alt_val - can_val) > 0.2:  # Significant change
                            signal_changes.append(
                                f"{sig_col}({can_val:.2f}→{alt_val:.2f})"
                            )

                records.append(
                    {
                        "gene": gene,
                        "transcript": alt_info["transcript"],
                        "feature_id": alt_info["feature_id"],
                        "feature_type": alt_info["feature_type"],
                        "canonical_loc": can_primary,
                        "canonical_confidence": round(can_confidence, 3),
                        "canonical_probs": can_top_probs,
                        "alternative_loc": alt_primary,
                        "alternative_confidence": round(alt_confidence, 3),
                        "alternative_probs": alt_top_probs,
                        "shift_magnitude": round(shift_magnitude, 4),
                        "length_difference_aa": int(length_diff),
                        "signal_change": "|".join(signal_changes)
                        if signal_changes
                        else "",
                        "loc_changed": can_primary != alt_primary,
                        "min_confidence": round(min(can_confidence, alt_confidence), 3),
                    }
                )

        if not records:
            logger.warning("No structural comparisons found")
            return views

        df = pd.DataFrame(records)

        # Create different ranked views
        # 1. By shift magnitude (largest biological changes first)
        views["by_shift_magnitude"] = df.sort_values(
            "shift_magnitude", ascending=False
        ).copy()

        # 2. By confidence (most reliable predictions first)
        views["by_confidence"] = df.sort_values(
            ["min_confidence", "shift_magnitude"], ascending=[False, False]
        ).copy()

        # 3. By extension length (longest extensions first)
        views["by_extension_length"] = df.sort_values(
            ["length_difference_aa", "shift_magnitude"], ascending=[False, False]
        ).copy()

        for name, view_df in views.items():
            logger.info(f"  Created {name}: {len(view_df)} rows")

        return views

    def create_missense_analysis_views(
        self,
        dataset,
        loc_results,
        protein_data,
        model_type,
        source,
        mutation_details=None,
        pair_results=None,
    ):
        """Create multiple ranked views of missense mutation effects (MUTATION analysis).

        Args:
            dataset: Dataset name.
            loc_results: Dictionary of localization DataFrames.
            protein_data: Dictionary of protein metadata.
            model_type: "Accurate" or "Fast".
            source: Mutation source (clinvar, cosmic, gnomad, custom).
            mutation_details: Optional dict of per-gene mutation DataFrames.
            pair_results: Optional isoform-level results DataFrame.

        Returns:
            Dict[str, pd.DataFrame]: Dictionary mapping view names to DataFrames.
        """
        logger.info(
            f"Creating missense analysis views for {dataset}/{source}/{model_type}"
        )

        views = {}

        # Get mutations data for the specified model
        if model_type.lower() == "accurate":
            mutations_data = loc_results.get("mutations_accurate")
            pairs_data = loc_results.get("pairs_accurate")
        else:
            mutations_data = loc_results.get("mutations_fast")
            pairs_data = loc_results.get("pairs_fast")

        if mutations_data is None or mutations_data.empty:
            logger.warning(f"No mutations data available for {model_type}")
            return views

        # Build comparison records
        records = []

        # Parse all protein IDs for mutations
        parsed_mutations = {}
        for _, row in mutations_data.iterrows():
            seq_id = row.get("Sequence_ID", "")
            if seq_id:
                parsed_mutations[seq_id] = self.parse_protein_id(seq_id)

        # Build reference localizations from pairs data (unmutated sequences)
        # Group by gene and feature_type
        reference_locs = {}
        if pairs_data is not None:
            for _, row in pairs_data.iterrows():
                seq_id = row.get("Sequence_ID", "")
                if seq_id:
                    info = self.parse_protein_id(seq_id)
                    if info["gene"] and not info["is_mutant"]:
                        key = (info["gene"], info["feature_type"])
                        reference_locs[key] = row

        # Also check mutations_data for reference sequences
        for _, row in mutations_data.iterrows():
            seq_id = row.get("Sequence_ID", "")
            if seq_id:
                info = parsed_mutations[seq_id]
                if info["gene"] and not info["is_mutant"]:
                    key = (info["gene"], info["feature_type"])
                    if key not in reference_locs:
                        reference_locs[key] = row

        # Compare mutants to their references
        for _, mut_row in mutations_data.iterrows():
            seq_id = mut_row.get("Sequence_ID", "")
            if not seq_id:
                continue

            info = parsed_mutations.get(seq_id)
            if not info or not info["is_mutant"]:
                continue

            gene = info["gene"]
            feature_type = info["feature_type"]
            mutation_info = info["mutation_info"]

            # Find reference
            ref_key = (gene, feature_type)
            if ref_key not in reference_locs:
                # Try canonical as fallback
                ref_key = (gene, "canonical")
                if ref_key not in reference_locs:
                    continue

            ref_row = reference_locs[ref_key]

            ref_confidence = self.get_primary_confidence(ref_row)
            ref_primary = self.get_primary_localization(ref_row)
            ref_top_probs = self.format_top_probs(ref_row, n=3)

            mut_confidence = self.get_primary_confidence(mut_row)
            mut_primary = self.get_primary_localization(mut_row)
            mut_top_probs = self.format_top_probs(mut_row, n=3)

            shift_magnitude = self.compute_shift_magnitude(ref_row, mut_row)

            # Extract variant_id from mutation_info or sequence_id
            variant_id = ""
            hgvsp = ""
            hgvsc = ""
            if mutation_info:
                variant_id = f"{mutation_info['position']}_{mutation_info['change']}_{mutation_info['aa_change']}"
                hgvsp = mutation_info.get("aa_change", "")

            # Get clinical metadata from mutation_details if available
            clinical_significance = ""
            star_rating = None
            sample_count = None
            allele_frequency = None

            if mutation_details and mutation_info:
                # Look up mutation details
                for key, detail_df in mutation_details.items():
                    if key[0] == gene and not detail_df.empty:
                        # Try to match by position
                        pos_matches = detail_df[
                            detail_df["position"].astype(str)
                            == str(mutation_info["position"])
                        ]
                        if not pos_matches.empty:
                            match_row = pos_matches.iloc[0]
                            clinical_significance = match_row.get(
                                "clinvar_clinical_significance", ""
                            )
                            star_rating = match_row.get("clinvar_star_rating")
                            sample_count = match_row.get("cosmic_sample_count")
                            allele_frequency = match_row.get("gnomad_allele_frequency")
                            hgvsc = match_row.get("hgvsc", "")
                            hgvsp = match_row.get("hgvsp", hgvsp)
                            variant_id = match_row.get("variant_id", variant_id)
                            break

            records.append(
                {
                    "gene": gene,
                    "transcript": info["transcript"],
                    "feature_type": feature_type,
                    "variant_id": variant_id,
                    "hgvsp": hgvsp,
                    "hgvsc": hgvsc,
                    "canonical_loc": ref_primary,
                    "canonical_confidence": round(ref_confidence, 3),
                    "canonical_probs": ref_top_probs,
                    "mutated_loc": mut_primary,
                    "mutated_confidence": round(mut_confidence, 3),
                    "mutated_probs": mut_top_probs,
                    "shift_magnitude": round(shift_magnitude, 4),
                    "loc_changed": ref_primary != mut_primary,
                    "min_confidence": round(min(ref_confidence, mut_confidence), 3),
                    "source": source,
                    "clinvar_clinical_significance": clinical_significance,
                    "clinvar_star_rating": star_rating
                    if pd.notna(star_rating)
                    else None,
                    "cosmic_sample_count": sample_count
                    if pd.notna(sample_count)
                    else None,
                    "gnomad_allele_frequency": allele_frequency
                    if pd.notna(allele_frequency)
                    else None,
                }
            )

        if not records:
            logger.warning("No missense comparisons found")
            return views

        df = pd.DataFrame(records)

        # Universal rankings (all sources)
        # 1. By shift magnitude
        views["by_shift_magnitude"] = df.sort_values(
            "shift_magnitude", ascending=False
        ).copy()

        # 2. By model agreement (both confident + largest shift)
        views["by_confidence"] = df.sort_values(
            ["min_confidence", "shift_magnitude"], ascending=[False, False]
        ).copy()

        # 3. By mutation count per gene (aggregate later if needed)
        # For now, just sort by presence of loc change
        views["by_loc_change"] = df.sort_values(
            ["loc_changed", "shift_magnitude"], ascending=[False, False]
        ).copy()

        # Source-specific rankings
        source_lower = source.lower()

        if source_lower == "clinvar" or source_lower.startswith("clinvar"):
            # Clinical significance ranking
            significance_order = {
                "Pathogenic": 0,
                "Pathogenic/Likely pathogenic": 1,
                "Likely pathogenic": 2,
                "Uncertain significance": 3,
                "Conflicting classifications of pathogenicity": 4,
                "Likely benign": 5,
                "Benign/Likely benign": 6,
                "Benign": 7,
            }

            df_clin = df.copy()
            df_clin["_sig_order"] = df_clin["clinvar_clinical_significance"].map(
                lambda x: significance_order.get(str(x), 99)
            )
            views["by_clinical_significance"] = df_clin.sort_values(
                ["_sig_order", "shift_magnitude"], ascending=[True, False]
            ).drop(columns=["_sig_order"])

            # Star rating ranking
            if df["clinvar_star_rating"].notna().any():
                views["by_star_rating"] = df.sort_values(
                    ["clinvar_star_rating", "shift_magnitude"], ascending=[False, False]
                ).copy()

        if source_lower == "cosmic" or source_lower.startswith("cosmic"):
            # Sample count ranking
            if df["cosmic_sample_count"].notna().any():
                views["by_sample_count"] = df.sort_values(
                    ["cosmic_sample_count", "shift_magnitude"], ascending=[False, False]
                ).copy()

        if source_lower == "gnomad":
            # Rarity ranking (lowest frequency = rarest)
            if df["gnomad_allele_frequency"].notna().any():
                df_gnom = df[df["gnomad_allele_frequency"].notna()].copy()
                views["by_rarity"] = df_gnom.sort_values(
                    ["gnomad_allele_frequency", "shift_magnitude"],
                    ascending=[True, False],
                )

        for name, view_df in views.items():
            logger.info(f"  Created {name}: {len(view_df)} rows")

        return views

    def create_alt_start_loss_views(
        self,
        dataset,
        pair_results,
        loc_results,
        protein_data,
        model_type,
        mutation_details=None,
    ):
        """Create multiple ranked views of alt-start-loss mutations (MUTATION analysis).

        These are mutations that destroy the alternative start codon, effectively
        switching the cell to use the canonical start.

        Args:
            dataset: Dataset name.
            pair_results: Isoform-level results DataFrame with alt_start_loss info.
            loc_results: Dictionary of localization DataFrames.
            protein_data: Dictionary of protein metadata.
            model_type: "Accurate" or "Fast".
            mutation_details: Optional dict of per-gene mutation DataFrames.

        Returns:
            Dict[str, pd.DataFrame]: Dictionary mapping view names to DataFrames.
        """
        logger.info(f"Creating alt-start-loss views for {dataset}/{model_type}")

        views = {}

        if pair_results is None or pair_results.empty:
            logger.warning("No pair_results available for alt-start-loss analysis")
            return views

        # Get pairs data for localization comparison
        if model_type.lower() == "accurate":
            pairs_data = loc_results.get("pairs_accurate")
        else:
            pairs_data = loc_results.get("pairs_fast")

        if pairs_data is None:
            logger.warning(f"No pairs data available for {model_type}")
            return views

        # Filter to features with alt-start-loss mutations
        alt_loss_features = pair_results[
            pair_results["alternative_start_loss_count"] > 0
        ]

        if alt_loss_features.empty:
            logger.info("No features with alt-start-loss mutations found")
            return views

        # Build localization lookup
        loc_lookup = {}
        for _, row in pairs_data.iterrows():
            seq_id = row.get("Sequence_ID", "")
            if seq_id:
                info = self.parse_protein_id(seq_id)
                key = (info["gene"], info["feature_type"])
                loc_lookup[key] = {
                    "loc": self.get_primary_localization(row),
                    "confidence": self.get_primary_confidence(row),
                    "probs": self.format_top_probs(row, n=3),
                    "row": row,
                }

        records = []

        for _, feat_row in alt_loss_features.iterrows():
            gene = feat_row["gene_name"]
            transcript = feat_row["transcript_id"]
            feature_id = feat_row.get("feature_id", "")
            feature_type = feat_row.get("feature_type", "")

            # Get extension/truncation localization
            ext_key = (gene, feature_type)
            ext_info = loc_lookup.get(ext_key)

            # Get canonical localization
            can_key = (gene, "canonical")
            can_info = loc_lookup.get(can_key)

            if not ext_info or not can_info:
                continue

            locs_different = ext_info["loc"] != can_info["loc"]

            # Get variant IDs
            alt_start_loss_ids = feat_row.get("alternative_start_loss_variant_ids", "")
            if pd.isna(alt_start_loss_ids):
                alt_start_loss_ids = ""

            # Compute shift if localizations differ
            shift_magnitude = 0.0
            if "row" in ext_info and "row" in can_info:
                shift_magnitude = self.compute_shift_magnitude(
                    can_info["row"], ext_info["row"]
                )

            records.append(
                {
                    "gene": gene,
                    "transcript": transcript,
                    "feature_id": feature_id,
                    "feature_type": feature_type,
                    "extension_loc": ext_info["loc"],
                    "extension_confidence": round(ext_info["confidence"], 3),
                    "extension_probs": ext_info["probs"],
                    "canonical_loc": can_info["loc"],
                    "canonical_confidence": round(can_info["confidence"], 3),
                    "canonical_probs": can_info["probs"],
                    "locs_different": locs_different,
                    "shift_magnitude": round(shift_magnitude, 4),
                    "extension_length_aa": int(
                        feat_row.get("aa_difference_from_canonical", 0) or 0
                    ),
                    "alternative_start_loss_count": int(
                        feat_row.get("alternative_start_loss_count", 0)
                    ),
                    "alt_start_loss_variant_ids": alt_start_loss_ids,
                }
            )

        if not records:
            logger.warning("No alt-start-loss records created")
            return views

        df = pd.DataFrame(records)

        # Create ranked views
        # 1. By localization difference (genes where losing extension matters)
        views["by_loc_difference"] = df.sort_values(
            ["locs_different", "shift_magnitude"], ascending=[False, False]
        ).copy()

        # 2. By mutation count
        views["by_mutation_count"] = df.sort_values(
            ["alternative_start_loss_count", "shift_magnitude"],
            ascending=[False, False],
        ).copy()

        # 3. By extension length
        views["by_extension_length"] = df.sort_values(
            ["extension_length_aa", "shift_magnitude"], ascending=[False, False]
        ).copy()

        for name, view_df in views.items():
            logger.info(f"  Created {name}: {len(view_df)} rows")

        return views

    def create_mutation_burden_views(
        self, dataset, pair_results, source, mutation_details=None
    ):
        """Create multiple ranked views of mutation burden per feature (MUTATION analysis).

        Args:
            dataset: Dataset name.
            pair_results: Isoform-level results DataFrame.
            source: Mutation source (clinvar, cosmic, gnomad, custom).
            mutation_details: Optional dict of per-gene mutation DataFrames.

        Returns:
            Dict[str, pd.DataFrame]: Dictionary mapping view names to DataFrames.
        """
        logger.info(f"Creating mutation burden views for {dataset}/{source}")

        views = {}

        if pair_results is None or pair_results.empty:
            logger.warning("No pair_results available for mutation burden analysis")
            return views

        # Filter to features with mutations
        features_with_mutations = pair_results[pair_results["total_mutations"] > 0]

        if features_with_mutations.empty:
            logger.info("No features with mutations found")
            return views

        records = []

        for _, row in features_with_mutations.iterrows():
            gene = row["gene_name"]
            transcript = row["transcript_id"]
            feature_id = row.get("feature_id", "")
            feature_type = row.get("feature_type", "")
            feature_length = (
                row.get("feature_length_aa", 0) or 1
            )  # Avoid division by zero

            # Count different mutation types
            missense = int(row.get("count_missense_variant", 0) or 0)
            nonsense = int(row.get("count_nonsense_variant", 0) or 0)
            frameshift = int(row.get("count_frameshift_variant", 0) or 0)
            total = int(row.get("total_mutations", 0) or 0)

            # Count mutations in alternative start region
            alt_start_mutations = int(row.get("count_in_alt_start_site", 0) or 0)
            can_start_mutations = int(row.get("count_in_canonical_start_site", 0) or 0)

            # Mutation density
            mutation_density = round(total / feature_length * 100, 4)  # per 100 aa

            # Count pathogenic/likely pathogenic
            pathogenic_count = 0
            likely_pathogenic_count = 0

            # Try to get from source-specific columns
            source_prefix = source.lower()
            if source_prefix.startswith("custom"):
                source_prefix = "custom"

            if f"count_{source_prefix}" in row.index:
                # We'd need mutation details to count pathogenic
                pass

            records.append(
                {
                    "gene": gene,
                    "transcript": transcript,
                    "feature_id": feature_id,
                    "feature_type": feature_type,
                    "total_mutations": total,
                    "missense_count": missense,
                    "nonsense_count": nonsense,
                    "frameshift_count": frameshift,
                    "pathogenic_count": pathogenic_count,
                    "likely_pathogenic_count": likely_pathogenic_count,
                    "count_in_alt_start_site": alt_start_mutations,
                    "count_in_canonical_start_site": can_start_mutations,
                    "feature_length_aa": int(feature_length),
                    "mutation_density": mutation_density,
                }
            )

        if not records:
            logger.warning("No mutation burden records created")
            return views

        df = pd.DataFrame(records)

        # Create ranked views
        # 1. By pathogenic count (if we have that data)
        if df["pathogenic_count"].sum() > 0:
            views["by_pathogenic_count"] = df.sort_values(
                ["pathogenic_count", "total_mutations"], ascending=[False, False]
            ).copy()

        # 2. By mutation density
        views["by_mutation_density"] = df.sort_values(
            ["mutation_density", "total_mutations"], ascending=[False, False]
        ).copy()

        # 3. By alt region mutations
        views["by_alt_region_mutations"] = df.sort_values(
            ["count_in_alt_start_site", "total_mutations"], ascending=[False, False]
        ).copy()

        # 4. By total mutations
        views["by_total_mutations"] = df.sort_values(
            "total_mutations", ascending=False
        ).copy()

        for name, view_df in views.items():
            logger.info(f"  Created {name}: {len(view_df)} rows")

        return views

    def create_model_comparison_views(self, dataset, loc_results, source):
        """Create model agreement analysis between Accurate and Fast models.

        Args:
            dataset: Dataset name.
            loc_results: Dictionary containing both accurate and fast predictions.
            source: Mutation source.

        Returns:
            Dict[str, pd.DataFrame]: Dictionary mapping view names to DataFrames.
        """
        logger.info(f"Creating model comparison views for {dataset}/{source}")

        views = {}

        # Get both model results
        accurate_pairs = loc_results.get("pairs_accurate")
        fast_pairs = loc_results.get("pairs_fast")

        if accurate_pairs is None or fast_pairs is None:
            logger.warning("Need both Accurate and Fast results for model comparison")
            return views

        # Also compare mutation results if available
        accurate_muts = loc_results.get("mutations_accurate")
        fast_muts = loc_results.get("mutations_fast")

        # Build comparison records
        records = []

        # Create lookup for fast results
        fast_lookup = {}
        for _, row in fast_pairs.iterrows():
            seq_id = row.get("Sequence_ID", "")
            if seq_id:
                fast_lookup[seq_id] = row

        if fast_muts is not None:
            for _, row in fast_muts.iterrows():
                seq_id = row.get("Sequence_ID", "")
                if seq_id:
                    fast_lookup[seq_id] = row

        # Compare accurate to fast
        for _, acc_row in accurate_pairs.iterrows():
            seq_id = acc_row.get("Sequence_ID", "")
            if not seq_id or seq_id not in fast_lookup:
                continue

            fast_row = fast_lookup[seq_id]
            info = self.parse_protein_id(seq_id)

            acc_loc = self.get_primary_localization(acc_row)
            acc_conf = self.get_primary_confidence(acc_row)
            fast_loc = self.get_primary_localization(fast_row)
            fast_conf = self.get_primary_confidence(fast_row)

            agreement = acc_loc == fast_loc
            avg_conf = (acc_conf + fast_conf) / 2
            min_conf = min(acc_conf, fast_conf)

            # Categorize confidence
            if min_conf >= 0.7:
                conf_category = "high"
            elif min_conf >= 0.5:
                conf_category = "medium"
            else:
                conf_category = "low"

            records.append(
                {
                    "gene": info["gene"],
                    "transcript": info["transcript"],
                    "feature_type": info["feature_type"],
                    "sequence_id": seq_id,
                    "accurate_loc": acc_loc,
                    "accurate_confidence": round(acc_conf, 3),
                    "fast_loc": fast_loc,
                    "fast_confidence": round(fast_conf, 3),
                    "agreement": agreement,
                    "avg_confidence": round(avg_conf, 3),
                    "min_confidence": round(min_conf, 3),
                    "confidence_category": conf_category,
                }
            )

        # Also process mutations if available
        if accurate_muts is not None:
            for _, acc_row in accurate_muts.iterrows():
                seq_id = acc_row.get("Sequence_ID", "")
                if not seq_id or seq_id not in fast_lookup:
                    continue

                # Skip if already processed from pairs
                if any(r["sequence_id"] == seq_id for r in records):
                    continue

                fast_row = fast_lookup[seq_id]
                info = self.parse_protein_id(seq_id)

                acc_loc = self.get_primary_localization(acc_row)
                acc_conf = self.get_primary_confidence(acc_row)
                fast_loc = self.get_primary_localization(fast_row)
                fast_conf = self.get_primary_confidence(fast_row)

                agreement = acc_loc == fast_loc
                avg_conf = (acc_conf + fast_conf) / 2
                min_conf = min(acc_conf, fast_conf)

                if min_conf >= 0.7:
                    conf_category = "high"
                elif min_conf >= 0.5:
                    conf_category = "medium"
                else:
                    conf_category = "low"

                records.append(
                    {
                        "gene": info["gene"],
                        "transcript": info["transcript"],
                        "feature_type": info["feature_type"],
                        "sequence_id": seq_id,
                        "accurate_loc": acc_loc,
                        "accurate_confidence": round(acc_conf, 3),
                        "fast_loc": fast_loc,
                        "fast_confidence": round(fast_conf, 3),
                        "agreement": agreement,
                        "avg_confidence": round(avg_conf, 3),
                        "min_confidence": round(min_conf, 3),
                        "confidence_category": conf_category,
                    }
                )

        if not records:
            logger.warning("No model comparison records created")
            return views

        df = pd.DataFrame(records)

        # Create different views
        # 1. All predictions with agreement status
        views["agreement_summary"] = df.sort_values(
            ["agreement", "min_confidence"], ascending=[False, False]
        ).copy()

        # 2. High confidence agreements (both >0.7)
        high_conf = df[(df["min_confidence"] >= 0.7) & (df["agreement"])].copy()
        if not high_conf.empty:
            views["high_confidence_agreement"] = high_conf.sort_values(
                "min_confidence", ascending=False
            )

        # 3. Disagreement cases
        disagreements = df[~df["agreement"]].copy()
        if not disagreements.empty:
            views["disagreement_cases"] = disagreements.sort_values(
                "avg_confidence", ascending=False
            )

        # 4. Low confidence cases
        low_conf = df[df["min_confidence"] < 0.5].copy()
        if not low_conf.empty:
            views["low_confidence_cases"] = low_conf.sort_values(
                "min_confidence", ascending=True
            )

        # Summary statistics
        total = len(df)
        agreed = df["agreement"].sum()
        high_conf_count = len(df[df["confidence_category"] == "high"])
        logger.info(
            f"  Model comparison: {total} sequences, {agreed} ({100 * agreed / total:.1f}%) agree, "
            f"{high_conf_count} high confidence"
        )

        for name, view_df in views.items():
            logger.info(f"  Created {name}: {len(view_df)} rows")

        return views

    def create_localization_confidence_profiles(self, dataset, loc_results, model_type):
        """Create full probability distribution profiles for all predictions.

        Args:
            dataset: Dataset name.
            loc_results: Dictionary of localization DataFrames.
            model_type: "Accurate" or "Fast".

        Returns:
            pd.DataFrame: DataFrame with all probability columns for each sequence.
        """
        logger.info(
            f"Creating localization confidence profiles for {dataset}/{model_type}"
        )

        # Get data for the specified model
        if model_type.lower() == "accurate":
            pairs_data = loc_results.get("pairs_accurate")
            mutations_data = loc_results.get("mutations_accurate")
        else:
            pairs_data = loc_results.get("pairs_fast")
            mutations_data = loc_results.get("mutations_fast")

        records = []

        # Process pairs
        if pairs_data is not None:
            for _, row in pairs_data.iterrows():
                seq_id = row.get("Sequence_ID", "")
                if not seq_id:
                    continue

                info = self.parse_protein_id(seq_id)

                record = {
                    "gene": info["gene"],
                    "transcript": info["transcript"],
                    "feature_type": info["feature_type"],
                    "sequence_id": seq_id,
                    "is_mutant": info["is_mutant"],
                    "primary_loc": self.get_primary_localization(row),
                    "primary_confidence": round(self.get_primary_confidence(row), 3),
                    "localisation_string": row.get("Localisation", ""),
                    "signals": row.get("Signals", ""),
                    "membrane_types": row.get("Membrane types", ""),
                }

                # Add all probability columns
                for col in self.all_prob_cols:
                    if col in row.index and pd.notna(row[col]):
                        record[col] = round(float(row[col]), 4)
                    else:
                        record[col] = None

                records.append(record)

        # Process mutations
        if mutations_data is not None:
            for _, row in mutations_data.iterrows():
                seq_id = row.get("Sequence_ID", "")
                if not seq_id:
                    continue

                # Skip if already in pairs
                if any(r["sequence_id"] == seq_id for r in records):
                    continue

                info = self.parse_protein_id(seq_id)

                record = {
                    "gene": info["gene"],
                    "transcript": info["transcript"],
                    "feature_type": info["feature_type"],
                    "sequence_id": seq_id,
                    "is_mutant": info["is_mutant"],
                    "primary_loc": self.get_primary_localization(row),
                    "primary_confidence": round(self.get_primary_confidence(row), 3),
                    "localisation_string": row.get("Localisation", ""),
                    "signals": row.get("Signals", ""),
                    "membrane_types": row.get("Membrane types", ""),
                }

                for col in self.all_prob_cols:
                    if col in row.index and pd.notna(row[col]):
                        record[col] = round(float(row[col]), 4)
                    else:
                        record[col] = None

                records.append(record)

        if records:
            df = pd.DataFrame(records)
            logger.info(f"  Created confidence profiles for {len(df)} sequences")
            return df

        return pd.DataFrame()

    def generate_summary_text(
        self, dataset, source, feature_dfs, gene_results, pair_results
    ):
        """Generate a unified text summary orienting the reader to feature_analysis.csv.

        Args:
            dataset: Dataset name (e.g., 'hela').
            source: Mutation source (e.g., 'default', 'clinvar', 'gnomad').
            feature_dfs: Dict of {model_name: feature_analysis DataFrame}.
            gene_results: Gene-level mutation results DataFrame (or None).
            pair_results: Isoform-level mutation results DataFrame (or None).

        Returns:
            List[str]: Lines of summary text.
        """
        lines = []
        lines.append(f"SWISSISOFORM SUMMARY — {dataset.upper()} / {source.upper()}")
        lines.append("=" * 60)
        lines.append("")

        # --- OVERVIEW ---
        # Use the first available model's feature_analysis for overview stats
        model_names = list(feature_dfs.keys())
        primary_model = model_names[0] if model_names else None
        primary_df = (
            feature_dfs.get(primary_model, pd.DataFrame())
            if primary_model
            else pd.DataFrame()
        )

        lines.append("OVERVIEW")

        if not primary_df.empty:
            total_features = len(primary_df)
            n_extensions = (primary_df["feature_type"] == "extension").sum()
            n_truncations = (primary_df["feature_type"] == "truncation").sum()
            lines.append(
                f"  Total features: {total_features} "
                f"({n_extensions} extensions, {n_truncations} truncations)"
            )
        else:
            lines.append("  No features available.")

        lines.append(
            f"  Models available: {', '.join(model_names) if model_names else 'None'}"
        )
        lines.append("")

        # --- LOCALIZATION SHIFTS (per model) ---
        for model_name in model_names:
            df = feature_dfs[model_name]
            if df.empty:
                continue

            lines.append(f"LOCALIZATION SHIFTS ({model_name} model)")

            total = len(df)
            n_differ = int(df["locs_differ"].sum())
            pct_differ = round(100 * n_differ / total, 1) if total > 0 else 0

            lines.append(
                f"  Features with different localization from canonical: "
                f"{n_differ} / {total} ({pct_differ}%)"
            )

            shifts = df["shift_magnitude"]
            lines.append(
                f"  Shift magnitude: mean={shifts.mean():.2f}, "
                f"median={shifts.median():.2f}, max={shifts.max():.2f}"
            )

            # Top transitions
            differ_df = df[df["locs_differ"]]
            if not differ_df.empty:
                transitions = (
                    differ_df["canonical_loc"] + " -> " + differ_df["feature_loc"]
                )
                transition_counts = transitions.value_counts().head(5)
                lines.append("  Top transitions:")
                for transition, count in transition_counts.items():
                    lines.append(f"    {transition}: {count} features")

            # Top shifted features
            top_shifted = df.nlargest(5, "shift_magnitude")
            if not top_shifted.empty:
                lines.append("  Top shifted features (by shift_magnitude):")
                for _, row in top_shifted.iterrows():
                    loc_change = (
                        f"{row['canonical_loc']} -> {row['feature_loc']}"
                        if row["locs_differ"]
                        else f"{row['canonical_loc']} (same)"
                    )
                    lines.append(
                        f"    {row['gene']} ({row['feature_type']}): "
                        f"{row['shift_magnitude']:.2f} shift, {loc_change}"
                    )

            lines.append("")

        # --- MODEL AGREEMENT ---
        if len(model_names) == 2:
            lines.append("MODEL AGREEMENT")
            # Use primary model's df which has model_agreement column
            if not primary_df.empty and "model_agreement" in primary_df.columns:
                agreement_col = primary_df["model_agreement"]
                non_null = agreement_col.dropna()
                if len(non_null) > 0:
                    n_agree = int(non_null.sum())
                    n_total = len(non_null)
                    pct_agree = round(100 * n_agree / n_total, 1)
                    lines.append(
                        f"  Both models agree: {n_agree} / {n_total} ({pct_agree}%)"
                    )
                    lines.append(f"  Disagreement cases: {n_total - n_agree}")
                else:
                    lines.append("  No model agreement data available.")
            else:
                lines.append("  No model agreement data available.")
            lines.append("")

        # --- MUTATION BURDEN (non-default sources only) ---
        if source != "default":
            lines.append("MUTATION BURDEN")

            if gene_results is not None and not gene_results.empty:
                total_genes = len(gene_results)
                successful = len(gene_results[gene_results["status"] == "success"])
                lines.append(
                    f"  Genes processed: {total_genes} ({successful} successful)"
                )

            if pair_results is not None and not pair_results.empty:
                detected_sources = self.detect_mutation_sources(pair_results)
                if detected_sources:
                    lines.append(f"  Sources detected: {', '.join(detected_sources)}")

                if "total_mutations" in pair_results.columns:
                    total_muts = int(pair_results["total_mutations"].sum())
                    lines.append(f"  Total mutations across features: {total_muts}")

                # Breakdown by type
                variant_types = [
                    ("missense_variant", "missense"),
                    ("nonsense_variant", "nonsense"),
                    ("frameshift_variant", "frameshift"),
                    ("inframe_deletion", "inframe deletion"),
                    ("inframe_insertion", "inframe insertion"),
                    ("synonymous_variant", "synonymous"),
                ]
                breakdown_parts = []
                for col_suffix, label in variant_types:
                    col = f"count_{col_suffix}"
                    if col in pair_results.columns:
                        count = int(pair_results[col].sum())
                        if count > 0:
                            breakdown_parts.append(f"{count} {label}")
                if breakdown_parts:
                    lines.append(f"  Breakdown: {', '.join(breakdown_parts)}")

                # Feature-level mutation info from feature_analysis
                if not primary_df.empty:
                    if "missense_loc_change_count" in primary_df.columns:
                        n_with_missense_loc = int(
                            (primary_df["missense_loc_change_count"] > 0).sum()
                        )
                        lines.append(
                            f"  Features with missense-induced localization changes: {n_with_missense_loc}"
                        )
                    if "alt_start_loss_count" in primary_df.columns:
                        n_with_alt_loss = int(
                            (primary_df["alt_start_loss_count"] > 0).sum()
                        )
                        lines.append(
                            f"  Features with alt-start-loss mutations: {n_with_alt_loss}"
                        )
            else:
                lines.append("  No mutation data available.")

            lines.append("")

        # --- GUIDE TO feature_analysis.csv ---
        lines.append("GUIDE TO feature_analysis.csv")
        lines.append("  Key columns for localization:")
        lines.append("    - shift_magnitude: max probability change vs canonical (0-1)")
        lines.append("    - locs_differ: True if primary localization changed")
        lines.append("    - signal_change: membrane/signal type differences")
        lines.append("    - model_agreement: True if Accurate and Fast models agree")

        if source != "default":
            lines.append("  Key columns for mutations:")
            lines.append("    - total_mutations, mutation_density")
            lines.append("    - missense_loc_change_count")
            lines.append("    - alt_start_loss_count")
            if source == "clinvar":
                lines.append("    - clinvar_pathogenic_count, clinvar_top_star_rating")
            elif source == "cosmic":
                lines.append("    - cosmic_max_sample_count")
            elif source == "gnomad":
                lines.append("    - gnomad_min_af")

        lines.append("")

        return lines

    def analyze_dataset(self, dataset, source="gnomad"):
        """Analyze a complete dataset and source, saving unified feature analysis.

        Creates a single comprehensive feature_analysis.csv per model that combines
        all metrics (structural, missense, alt-start-loss, mutation burden).

        Output structure:
        - results/{dataset}/{source}/summary/
          - summary.txt (text overview)
          - model_comparison/ (if both models available)
          - accurate/
            - feature_analysis.csv (unified analysis table)
          - fast/
            - feature_analysis.csv (unified analysis table)

        Args:
            dataset (str): Name of the dataset to analyze.
            source (str): Mutation source (default, gnomad, clinvar, cosmic, custom_*).
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

        # Load detailed mutation CSVs if available (for mutation sources)
        mutation_details = None
        if source != "default":
            mutation_details = self.load_mutation_details(dataset, source)

        # Get available models
        available_models = self.get_available_models(dataset, source)

        if not available_models:
            logger.warning(f"No localization models available for {dataset}/{source}")
            return

        logger.info(
            f"Available models for {dataset}/{source}: {', '.join(available_models)}"
        )

        # Create model comparison views if both models available
        if len(available_models) == 2:
            logger.info("\n=== CREATING MODEL COMPARISON VIEWS ===")
            comparison_views = self.create_model_comparison_views(
                dataset, loc_results, source
            )

            # Save model comparison views
            comparison_dir = summary_dir / "model_comparison"
            comparison_dir.mkdir(parents=True, exist_ok=True)

            for view_name, view_df in comparison_views.items():
                view_df.to_csv(comparison_dir / f"{view_name}.csv", index=False)
                logger.info(
                    f"  Saved model_comparison/{view_name}.csv ({len(view_df)} rows)"
                )

        # Analyze each available model separately, collecting feature DataFrames
        feature_dfs = {}

        for model_type in available_models:
            logger.info(
                f"\n=== ANALYZING {model_type.upper()} MODEL FOR {dataset.upper()}/{source.upper()} ==="
            )

            # Create model-specific subdirectory
            model_summary_dir = summary_dir / model_type.lower()
            model_summary_dir.mkdir(parents=True, exist_ok=True)

            # Determine if we have the other model for agreement column
            other_model_loc = loc_results if len(available_models) == 2 else None

            # Create unified feature analysis table
            logger.info(f"\n--- Creating unified feature analysis ({model_type}) ---")
            feature_analysis_df = self.create_feature_analysis(
                dataset=dataset,
                loc_results=loc_results,
                protein_data=protein_data,
                model_type=model_type,
                source=source,
                pair_results=pair_results,
                mutation_details=mutation_details,
                other_model_loc=other_model_loc,
            )

            # Save feature analysis
            if not feature_analysis_df.empty:
                feature_analysis_df.to_csv(
                    model_summary_dir / "feature_analysis.csv", index=False
                )
                logger.info(
                    f"  Saved feature_analysis.csv ({len(feature_analysis_df)} features)"
                )

                # Print summary statistics
                locs_differ_count = feature_analysis_df["locs_differ"].sum()
                total_features = len(feature_analysis_df)
                logger.info(
                    f"  {locs_differ_count}/{total_features} features have different localization than canonical"
                )

                if source != "default":
                    with_mutations = (feature_analysis_df["total_mutations"] > 0).sum()
                    with_alt_start_loss = (
                        feature_analysis_df["alt_start_loss_count"] > 0
                    ).sum()
                    with_missense_loc_change = (
                        feature_analysis_df["missense_loc_change_count"] > 0
                    ).sum()
                    logger.info(f"  {with_mutations} features with mutations")
                    logger.info(
                        f"  {with_alt_start_loss} features with alt-start-loss mutations"
                    )
                    logger.info(
                        f"  {with_missense_loc_change} features with missense-induced loc changes"
                    )
            else:
                # Save empty file
                pd.DataFrame().to_csv(
                    model_summary_dir / "feature_analysis.csv", index=False
                )

            # Collect for summary text generation
            feature_dfs[model_type] = feature_analysis_df

        # Generate and save unified summary text
        summary_lines = self.generate_summary_text(
            dataset, source, feature_dfs, gene_results, pair_results
        )
        with open(summary_dir / "summary.txt", "w") as f:
            f.write("\n".join(summary_lines))
        logger.info("Saved summary.txt")

        logger.info(f"\n{'=' * 60}")
        logger.info(f"Summary analysis for {dataset}/{source} completed!")
        logger.info(f"Results saved to: {summary_dir}")
        logger.info(f"{'=' * 60}")
        for model_type in available_models:
            logger.info(
                f"  {model_type} model: {summary_dir / model_type.lower()}/feature_analysis.csv"
            )
        if len(available_models) == 2:
            logger.info(f"  Model comparison: {summary_dir / 'model_comparison'}/")

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
