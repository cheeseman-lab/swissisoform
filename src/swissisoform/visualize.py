"""Genomic feature visualization module.

This module provides the GenomeVisualizer class for creating visual
representations of transcript features, alternative isoforms, and
mutation data.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from typing import Optional, List, Dict
from pathlib import Path


class GenomeVisualizer:
    """Visualize genomic features and mutations along a transcript.

    This class provides methods to create visualizations of transcripts
    with optional alternative start sites and mutation data.
    """

    def __init__(self, genome_handler):
        """Initialize the GenomeVisualizer.

        Args:
            genome_handler: Initialized GenomeHandler instance
        """
        self.genome = genome_handler

        # Feature colors
        self.feature_colors = {
            "CDS": "#cacfd2",  # light gray (as specified)
            "UTR": "#a6acaf",  # medium gray (as specified)
            "start_codon": "#f27e20",  # annotated start site (as specified)
            "stop_codon": "#e74c3c",  # black (as specified)
            "truncation": "#412D78",  # purple (as specified)
            "extension": "#2E8B57",  # green for extension
            "alternative_start": "#652d90",  # alternative (as specified)
        }

        # Mutation impact colors (updated to distinct random colors)
        self.mutation_colors = {
            "missense variant": "#1b9e77",  # bright red
            "synonymous variant": "#d95f02",  # bright blue
            "nonsense variant": "#d95f02",  # blue
            "inframe deletion": "#FF9800",  # orange
            "inframe insertion": "#9C27B0",  # purple
            "frameshift variant": "#7570b3",  # yellow
            "splice variant": "#00BCD4",  # cyan
            "start lost variant": "#795548",  # brown
            "5 prime utr variant": "#607D8B",  # blue-gray
            "3 prime utr variant": "#FF5722",  # deep orange
            "other variant": "#9E9E9E",  # medium gray
            "unknown": "#000000",  # black (kept original)
            "extension": "#2E8B57",  # green for extension
        }

        # Source-specific markers
        self.source_markers = {
            "gnomAD": "o",  # circle
            "ClinVar": "x",  # heavy asterisk
            "COSMIC": "s",  # square
        }

        # Track dimensions
        self.track_height = 0.15
        self.codon_height = 0.2

        # Y-positions for mutation tracks
        self.mutation_track_positions = {
            "gnomAD": 0.7,
            "ClinVar": 0.8,
            "COSMIC": 0.9,
        }

    def _get_mutation_color(self, impact: str) -> str:
        """Get color for mutation impact, handling various format cases.

        Args:
            impact: Impact type string

        Returns:
            Hex color code
        """
        if pd.isna(impact) or impact is None:
            return self.mutation_colors["unknown"]

        impact = str(impact).strip().lower()

        # Direct match with standardized categories
        if impact in self.mutation_colors:
            return self.mutation_colors[impact]

        return self.mutation_colors["other variant"]

    def _create_mutation_legend_groups(self, unique_impacts: List[str]) -> List[Line2D]:
        """Create legend elements for mutation impacts.

        Args:
            unique_impacts: List of unique impact types

        Returns:
            List of legend elements
        """
        # Simplified legend groups matching standardized categories
        legend_groups = {
            "Missense": ["missense variant"],
            "Synonymous": ["synonymous variant"],
            "Nonsense": ["nonsense variant"],
            "Inframe": ["inframe deletion", "inframe insertion"],
            "Frameshift": ["frameshift variant"],
            "Splice": ["splice variant"],
            "UTR": ["5 prime utr variant", "3 prime utr variant"],
            "Other": ["other variant", "unknown"],
            "Extension": ["extension"],
        }

        legend_elements = []

        # Create legend elements for each group that has variants
        for group_name, impact_types in legend_groups.items():
            group_impacts = [
                imp for imp in unique_impacts if str(imp).lower() in impact_types
            ]

            if group_impacts:
                # Use the first impact's color for the group
                color = self._get_mutation_color(group_impacts[0])
                legend_elements.append(
                    Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="w",
                        markerfacecolor=color,
                        label=group_name,
                        markersize=8,
                    )
                )

        return legend_elements

    def _preprocess_features(self, features):
        """Preprocess features to ensure proper visualization with connected elements. Makes minimal adjustments to maintain biological accuracy while enhancing visualization.

        Args:
            features: DataFrame of transcript features
        Returns:
            Preprocessed features DataFrame
        """
        # Make a copy to avoid modifying the original
        features = features.copy()

        # Constants for maximum allowed adjustments (in base pairs)
        MAX_UTR_ADJUSTMENT = 3
        MAX_CODON_ADJUSTMENT = 3
        MAX_CDS_ADJUSTMENT = 0  # No CDS modifications by default

        # Get transcript strand
        transcript_strand = "+"
        transcript_info = features[features["feature_type"] == "transcript"]
        if not transcript_info.empty:
            transcript_strand = transcript_info.iloc[0]["strand"]

        # Find relevant features
        utr_features = features[features["feature_type"] == "UTR"]
        start_codons = features[features["feature_type"] == "start_codon"]
        stop_codons = features[features["feature_type"] == "stop_codon"]
        cds_features = features[features["feature_type"] == "CDS"]

        # Classify UTRs based on their position relative to CDS
        if not utr_features.empty and not cds_features.empty:
            # Get the CDS boundaries
            cds_start = cds_features["start"].min()
            cds_end = cds_features["end"].max()

            # Create empty dataframes for 5' and 3' UTRs
            utr_5prime = pd.DataFrame()
            utr_3prime = pd.DataFrame()

            # Classify each UTR as 5' or 3' based on strand and position
            for idx, utr in utr_features.iterrows():
                if transcript_strand == "+":
                    if utr["end"] <= cds_start:
                        # This is a 5' UTR (before CDS in + strand)
                        utr_5prime = pd.concat([utr_5prime, utr.to_frame().T])
                    elif utr["start"] >= cds_end:
                        # This is a 3' UTR (after CDS in + strand)
                        utr_3prime = pd.concat([utr_3prime, utr.to_frame().T])
                else:  # negative strand
                    if utr["start"] >= cds_end:
                        # This is a 5' UTR (after CDS in - strand)
                        utr_5prime = pd.concat([utr_5prime, utr.to_frame().T])
                    elif utr["end"] <= cds_start:
                        # This is a 3' UTR (before CDS in - strand)
                        utr_3prime = pd.concat([utr_3prime, utr.to_frame().T])
        else:
            # If no CDS or UTR features, create empty dataframes
            utr_5prime = pd.DataFrame()
            utr_3prime = pd.DataFrame()

        # Handle 5'UTR and start codon connection
        if not utr_5prime.empty and not start_codons.empty:
            for utr_idx in utr_5prime.index:
                utr = features.loc[utr_idx]

                for _, start_codon in start_codons.iterrows():
                    # For positive strand
                    if transcript_strand == "+":
                        gap = start_codon["start"] - utr["end"]
                        if 0 < gap <= MAX_UTR_ADJUSTMENT:
                            features.loc[utr_idx, "end"] = start_codon["start"]
                    # For negative strand
                    elif transcript_strand == "-":
                        gap = utr["start"] - start_codon["end"]
                        if 0 < gap <= MAX_UTR_ADJUSTMENT:
                            features.loc[utr_idx, "start"] = start_codon["end"]

        # Handle start codon and CDS connection (only adjust start codon, not CDS)
        if not start_codons.empty and not cds_features.empty:
            for start_idx in start_codons.index:
                start_codon = features.loc[start_idx]

                if transcript_strand == "+":
                    # On positive strand, connect to first CDS
                    first_cds = cds_features.sort_values("start").iloc[0]
                    gap = first_cds["start"] - start_codon["end"]

                    # If there's a small gap, extend start codon (not CDS)
                    if 0 < gap <= MAX_CODON_ADJUSTMENT:
                        features.loc[start_idx, "end"] = first_cds["start"]
                else:
                    # On negative strand, connect to last CDS (by end position)
                    last_cds = cds_features.sort_values("end", ascending=False).iloc[0]
                    gap = start_codon["start"] - last_cds["end"]

                    # If there's a small gap, adjust start codon (not CDS)
                    if 0 < gap <= MAX_CODON_ADJUSTMENT:
                        features.loc[start_idx, "start"] = last_cds["end"]

        # Handle last CDS and stop codon connection (only adjust stop codon, not CDS)
        if not stop_codons.empty and not cds_features.empty:
            for stop_idx in stop_codons.index:
                stop_codon = features.loc[stop_idx]

                if transcript_strand == "+":
                    last_cds = cds_features.sort_values("end", ascending=False).iloc[0]
                    gap = stop_codon["start"] - last_cds["end"]

                    # If there's a small gap, adjust stop codon
                    if 0 < gap <= MAX_CODON_ADJUSTMENT:
                        features.loc[stop_idx, "start"] = last_cds["end"]
                else:
                    # For negative strand, connect to first CDS (by start position)
                    first_cds = cds_features.sort_values("start").iloc[0]
                    gap = first_cds["start"] - stop_codon["end"]

                    # If there's a small gap, adjust stop codon
                    if 0 < gap <= MAX_CODON_ADJUSTMENT:
                        features.loc[stop_idx, "end"] = first_cds["start"]

        # Handle stop codon and 3'UTR connection
        if not utr_3prime.empty and not stop_codons.empty:
            for utr_idx in utr_3prime.index:
                utr = features.loc[utr_idx]

                for _, stop_codon in stop_codons.iterrows():
                    # For positive strand
                    if transcript_strand == "+":
                        gap = utr["start"] - stop_codon["end"]
                        if 0 < gap <= MAX_UTR_ADJUSTMENT:
                            features.loc[utr_idx, "start"] = stop_codon["end"]
                    # For negative strand
                    elif transcript_strand == "-":
                        gap = stop_codon["start"] - utr["end"]
                        if 0 < gap <= MAX_UTR_ADJUSTMENT:
                            features.loc[utr_idx, "end"] = stop_codon["start"]

        # Log any changes made
        if hasattr(self, "logger") and self.logger:
            self.logger.debug(
                f"Feature preprocessing completed, adjustments limited to: UTR={MAX_UTR_ADJUSTMENT}bp, Codon={MAX_CODON_ADJUSTMENT}bp, CDS={MAX_CDS_ADJUSTMENT}bp"
            )

        return features

    def visualize_transcript(
        self,
        gene_name: str,
        transcript_id: str,
        alt_features: Optional[pd.DataFrame] = None,
        mutations_df: Optional[pd.DataFrame] = None,
        output_file: Optional[str] = None,
    ) -> None:
        """Visualize transcript with optional alternative start sites and mutations.

        Args:
            gene_name: Name of the gene
            transcript_id: Transcript ID to visualize
            alt_features: Alternative start features
            mutations_df: Mutations data
            output_file: Path to save the visualization
        """
        transcript_data = self.genome.get_transcript_features_with_sequence(
            transcript_id
        )
        features = transcript_data["features"]

        features = self._preprocess_features(features)

        # Get transcript bounds and create figure
        transcript_info = features[features["feature_type"] == "transcript"].iloc[0]
        transcript_start = transcript_info["start"]
        transcript_end = transcript_info["end"]
        span = transcript_end - transcript_start

        fig, ax = plt.subplots(figsize=(15, 4))

        transcript_strand = "+"  # Default to positive strand
        transcript_info = features[features["feature_type"] == "transcript"]
        if not transcript_info.empty:
            transcript_strand = transcript_info.iloc[0]["strand"]

        # Plot transcript features
        for _, feature in features.iterrows():
            width = feature["end"] - feature["start"]

            if feature["feature_type"] == "CDS":
                rect = Rectangle(
                    (feature["start"], 0.325),
                    width,
                    self.track_height,
                    facecolor=self.feature_colors["CDS"],
                )
                ax.add_patch(rect)

            elif feature["feature_type"] == "UTR":
                rect = Rectangle(
                    (feature["start"], 0.325),
                    width,
                    self.track_height,
                    facecolor=self.feature_colors["UTR"],
                )
                ax.add_patch(rect)

            elif feature["feature_type"] in ["start_codon", "stop_codon"]:
                # Position the markers based on strand direction and codon type
                if transcript_strand == "+":
                    # On positive strand:
                    # - Start codons should be at their start position
                    # - Stop codons should be at their end position
                    if feature["feature_type"] == "start_codon":
                        codon_position = feature["start"]
                    else:  # stop_codon
                        codon_position = feature["end"]
                else:  # Negative strand
                    # On negative strand:
                    # - Start codons should be at their end position
                    # - Stop codons should be at their start position
                    if feature["feature_type"] == "start_codon":
                        codon_position = feature["end"]
                    else:  # stop_codon
                        codon_position = feature["start"]

                codon_y = 0.4  # Same as the transcript line
                codon_height = self.codon_height
                # Draw a vertical line
                ax.vlines(
                    x=codon_position,
                    ymin=codon_y - (codon_height / 2),
                    ymax=codon_y + (codon_height / 2),
                    color=self.feature_colors[feature["feature_type"]],
                    linewidth=3,  # Same width as alternative start
                )
                ax.add_patch(rect)

        # Plot alternative isoforms (truncation/extension) - update in both visualization methods
        if alt_features is not None and not alt_features.empty:
            for _, alt_feature in alt_features.iterrows():
                alt_start = alt_feature["start"]
                alt_end = alt_feature["end"]
                region_type = alt_feature.get("region_type", "alternative")

                # Get transcript strand
                transcript_strand = "+"  # Default to positive strand
                if features is not None and not features.empty:
                    transcript_info = features[features["feature_type"] == "transcript"]
                    if not transcript_info.empty:
                        transcript_strand = transcript_info.iloc[0]["strand"]

                # Draw a bracket for truncation/extension
                bracket_y = 0.2  # Position for the top of the bracket
                bracket_height = 0.05  # Height of the bracket

                bracket_vertices = [
                    (alt_start, bracket_y),
                    (alt_start, bracket_y - bracket_height),
                    (alt_end, bracket_y - bracket_height),
                    (alt_end, bracket_y),
                ]
                xs, ys = zip(*bracket_vertices)
                ax.plot(
                    xs,
                    ys,
                    color=self.feature_colors.get(
                        region_type, self.feature_colors["alternative_start"]
                    ),
                    linewidth=1.0,
                    solid_capstyle="butt",
                    label="Alternative Isoform",
                )

                # Alternative start marker positioning depends on strand
                black_bar_y = 0.4
                alt_start_height = self.codon_height
                alt_start_ymin = black_bar_y - (alt_start_height / 2)
                alt_start_ymax = black_bar_y + (alt_start_height / 2)

                if region_type == "extension":
                    if transcript_strand == "+":
                        alt_start_pos = alt_start  # Start of extension region
                    else:
                        alt_start_pos = (
                            alt_end  # End of extension region (start in - strand)
                        )
                else:  # truncation
                    if transcript_strand == "+":
                        alt_start_pos = alt_end  # End of truncated region
                    else:
                        alt_start_pos = (
                            alt_start  # Start of truncated region (end in - strand)
                        )

                ax.vlines(
                    x=alt_start_pos,
                    ymin=alt_start_ymin,
                    ymax=alt_start_ymax,
                    color=self.feature_colors.get(
                        region_type, self.feature_colors["alternative_start"]
                    ),
                    linewidth=3,
                    label="Alternative Isoform",
                )

                # Place the start codon label and region type
                if region_type == "extension":
                    label_x = alt_start if transcript_strand == "+" else alt_end
                else:  # truncation
                    label_x = alt_end if transcript_strand == "+" else alt_start
                label_text = alt_feature.get("start_codon", "")
                if region_type:
                    label_text = f"{region_type.capitalize()} {label_text}".strip()
                plt.text(
                    label_x,
                    alt_start_ymax + 0.05,
                    label_text,
                    fontsize=8,
                    rotation=45,
                    ha="center",
                )

        # Plot mutations if provided
        if mutations_df is not None and not mutations_df.empty:
            # Use impact_validated if available
            impact_col = (
                "impact_validated"
                if "impact_validated" in mutations_df.columns
                else "impact"
            )
            source_legend_elements = []
            unique_impacts = mutations_df[impact_col].unique()
            mutation_legend_elements = self._create_mutation_legend_groups(
                unique_impacts
            )

            jitter_amount = 0.02
            np.random.seed(42)

            for source, base_y_pos in self.mutation_track_positions.items():
                source_data = mutations_df[mutations_df["source"] == source]

                if not source_data.empty:
                    source_legend_elements.append(
                        Line2D(
                            [0],
                            [0],
                            marker=self.source_markers[source],
                            color="w",
                            markerfacecolor="black",
                            label=source,
                            markersize=6,
                        )
                    )

                    jittered_y = base_y_pos + np.random.uniform(
                        -jitter_amount, jitter_amount, len(source_data)
                    )

                    for (_, mutation), y_pos in zip(source_data.iterrows(), jittered_y):
                        # FIX: Use impact_col consistently
                        color = self._get_mutation_color(mutation[impact_col])
                        ax.scatter(
                            mutation["position"],
                            y_pos,
                            marker=self.source_markers[source],
                            c=color,
                            s=10,
                            alpha=1,
                        )

        # Create legends - DYNAMIC VERSION
        feature_legend_elements = [
            plt.Rectangle(
                (0, 0), 1, 1, facecolor=self.feature_colors["CDS"], label="CDS"
            ),
            plt.Rectangle(
                (0, 0), 1, 1, facecolor=self.feature_colors["UTR"], label="UTR"
            ),
            Line2D(
                [0],
                [0],
                marker="|",
                color=self.feature_colors["start_codon"],
                label="Start codon",
                markersize=12,
                markeredgewidth=3,
                linestyle="None",
            ),
            Line2D(
                [0],
                [0],
                marker="|",
                color=self.feature_colors["stop_codon"],
                label="Stop codon",
                markersize=12,
                markeredgewidth=3,
                linestyle="None",
            ),
        ]

        # Add alternative isoform legend elements based on what's actually present
        if alt_features is not None and not alt_features.empty:
            region_types = alt_features["region_type"].unique()

            for region_type in region_types:
                if region_type == "truncation":
                    feature_legend_elements.extend(
                        [
                            Line2D(
                                [0],
                                [0],
                                marker="|",
                                color=self.feature_colors["truncation"],
                                label="Truncation start",
                                markersize=12,
                                markeredgewidth=3,
                                linestyle="None",
                            ),
                            Line2D(
                                [0],
                                [0],
                                color=self.feature_colors["truncation"],
                                label="Truncation region",
                                linewidth=1.5,
                            ),
                        ]
                    )
                elif region_type == "extension":
                    feature_legend_elements.extend(
                        [
                            Line2D(
                                [0],
                                [0],
                                marker="|",
                                color=self.feature_colors["extension"],
                                label="Extension start",
                                markersize=12,
                                markeredgewidth=3,
                                linestyle="None",
                            ),
                            Line2D(
                                [0],
                                [0],
                                color=self.feature_colors["extension"],
                                label="Extension region",
                                linewidth=1.5,
                            ),
                        ]
                    )

        # Add legends with improved spacing and spine handling
        if mutations_df is not None and not mutations_df.empty:
            # First legend - Features
            l1 = ax.legend(
                handles=feature_legend_elements,
                bbox_to_anchor=(1.05, 1.05),
                loc="upper left",
                borderaxespad=0.0,
                frameon=False,
                title="Transcript Features",
                prop={"size": 8},
                title_fontsize=9,
            )

            # Second legend - Sources
            ax2 = ax.twinx()
            l2 = ax2.legend(
                handles=source_legend_elements,
                bbox_to_anchor=(1.05, 0.55),
                loc="upper left",
                borderaxespad=0.0,
                frameon=False,
                title="Mutation Sources",
                prop={"size": 8},
                title_fontsize=9,
            )

            # Third legend - Mutation Types
            ax3 = ax.twinx()
            l3 = ax3.legend(
                handles=mutation_legend_elements,
                bbox_to_anchor=(1.05, 0.3),
                loc="upper left",
                borderaxespad=0.0,
                frameon=False,
                title="Mutation Types",
                prop={"size": 8},
                title_fontsize=9,
            )

            # Hide all spines and ticks for twin axes
            for axis in [ax2, ax3]:
                axis.set_yticks([])
                for spine in axis.spines.values():
                    spine.set_visible(False)
        else:
            ax.legend(
                handles=feature_legend_elements,
                bbox_to_anchor=(1.05, 1),
                loc="upper left",
                borderaxespad=0.0,
                frameon=False,
            )

        # Customize plot
        plt.title(f"{gene_name} - {transcript_id}", pad=20, y=1.05)
        ax.set_ylim(0, 1)
        ax.set_xlim(transcript_start - 10, transcript_end + 10)

        # Set tick positions
        tick_span = transcript_end - transcript_start
        if tick_span > 100000:
            tick_interval = 10000
        elif tick_span > 50000:
            tick_interval = 5000
        elif tick_span > 10000:
            tick_interval = 1000
        else:
            tick_interval = 500

        base_position = transcript_start - (transcript_start % tick_interval)
        tick_positions = range(
            base_position, int(transcript_end) + tick_interval, tick_interval
        )

        plt.xticks(tick_positions, [f"{pos:,}" for pos in tick_positions])
        plt.xticks(rotation=45)

        # Draw the transcript line after setting the ticks
        xlim = ax.get_xlim()
        ax.hlines(
            y=0.4, xmin=xlim[0], xmax=xlim[1], color="black", linewidth=1.5, zorder=0
        )

        # Remove y-axis and spines
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        plt.tight_layout()

        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(
                output_path,
                bbox_inches="tight",
                dpi=300,
                facecolor="white",
                edgecolor="none",
                format="pdf",
            )
            plt.close()
        else:
            plt.show()

    def visualize_transcript_zoomed(
        self,
        gene_name: str,
        transcript_id: str,
        alt_features: pd.DataFrame,
        mutations_df: Optional[pd.DataFrame] = None,
        output_file: Optional[str] = None,
        padding: int = 50,
    ) -> None:
        """Visualize a zoomed-in region of the transcript centered on alt_features.

        Args:
            gene_name: Name of the gene
            transcript_id: Transcript ID to visualize
            alt_features: Alternative start features defining zoom region
            mutations_df: Mutations data
            output_file: Path to save the visualization
            padding: Number of bases to add as padding on each side
        """
        if alt_features is None or alt_features.empty:
            raise ValueError("alt_features is required for zoomed visualization")

        transcript_data = self.genome.get_transcript_features_with_sequence(
            transcript_id
        )
        features = transcript_data["features"]

        if features.empty:
            raise ValueError(f"No features found for transcript {transcript_id}")

        features = self._preprocess_features(features)

        # Calculate zoom region based on alt_features
        zoom_start = alt_features["start"].min() - padding
        zoom_end = alt_features["end"].max() + padding

        # Create figure
        fig, ax = plt.subplots(figsize=(15, 4))

        transcript_strand = "+"  # Default to positive strand
        transcript_info = features[features["feature_type"] == "transcript"]
        if not transcript_info.empty:
            transcript_strand = transcript_info.iloc[0]["strand"]

        # Plot transcript features (only those in view)
        for _, feature in features.iterrows():
            # Skip features entirely outside zoom region
            if feature["end"] < zoom_start or feature["start"] > zoom_end:
                continue

            width = feature["end"] - feature["start"]

            if feature["feature_type"] == "CDS":
                rect = Rectangle(
                    (feature["start"], 0.325),
                    width,
                    self.track_height,
                    facecolor=self.feature_colors["CDS"],
                )
                ax.add_patch(rect)

            elif feature["feature_type"] == "UTR":
                rect = Rectangle(
                    (feature["start"], 0.325),
                    width,
                    self.track_height,
                    facecolor=self.feature_colors["UTR"],
                )
                ax.add_patch(rect)

            elif feature["feature_type"] in ["start_codon", "stop_codon"]:
                # Position the markers based on strand direction and codon type
                if transcript_strand == "+":
                    # On positive strand:
                    # - Start codons should be at their start position
                    # - Stop codons should be at their end position
                    if feature["feature_type"] == "start_codon":
                        codon_position = feature["start"]
                    else:  # stop_codon
                        codon_position = feature["end"]
                else:  # Negative strand
                    # On negative strand:
                    # - Start codons should be at their end position
                    # - Stop codons should be at their start position
                    if feature["feature_type"] == "start_codon":
                        codon_position = feature["end"]
                    else:  # stop_codon
                        codon_position = feature["start"]

                codon_y = 0.4  # Same as the transcript line
                codon_height = self.codon_height

                # Draw a vertical line
                ax.vlines(
                    x=codon_position,
                    ymin=codon_y - (codon_height / 2),
                    ymax=codon_y + (codon_height / 2),
                    color=self.feature_colors[feature["feature_type"]],
                    linewidth=3,  # Same width as alternative start
                )

        # Plot alternative isoform regions (truncation/extension) with unified labeling
        if alt_features is not None and not alt_features.empty:
            for _, alt_feature in alt_features.iterrows():
                alt_start = alt_feature["start"]
                alt_end = alt_feature["end"]
                region_type = alt_feature.get("region_type", "alternative")

                # Get transcript strand
                transcript_strand = "+"  # Default to positive strand
                if features is not None and not features.empty:
                    transcript_info = features[features["feature_type"] == "transcript"]
                    if not transcript_info.empty:
                        transcript_strand = transcript_info.iloc[0]["strand"]

                # Draw a bracket for truncation/extension
                bracket_y = 0.2  # Position for the top of the bracket
                bracket_height = 0.05  # Height of the bracket

                bracket_vertices = [
                    (alt_start, bracket_y),
                    (alt_start, bracket_y - bracket_height),
                    (alt_end, bracket_y - bracket_height),
                    (alt_end, bracket_y),
                ]
                xs, ys = zip(*bracket_vertices)
                ax.plot(
                    xs,
                    ys,
                    color=self.feature_colors.get(
                        region_type, self.feature_colors["alternative_start"]
                    ),
                    linewidth=1.0,
                    solid_capstyle="butt",
                    label="Alternative Isoform",
                )

                # Alternative start marker positioning depends on strand and region type
                black_bar_y = 0.4
                alt_start_height = self.codon_height
                alt_start_ymin = black_bar_y - (alt_start_height / 2)
                alt_start_ymax = black_bar_y + (alt_start_height / 2)

                # FIX: Correct positioning logic for extensions vs truncations
                if region_type == "extension":
                    # For extensions, the alternative start is at the beginning of the extension
                    if transcript_strand == "+":
                        alt_start_pos = alt_start  # Start of extension region
                    else:
                        alt_start_pos = (
                            alt_end  # End of extension region (start in - strand)
                        )
                else:  # truncation
                    # For truncations, the alternative start is at the end of the truncated region
                    if transcript_strand == "+":
                        alt_start_pos = alt_end  # End of truncated region
                    else:
                        alt_start_pos = (
                            alt_start  # Start of truncated region (end in - strand)
                        )

                ax.vlines(
                    x=alt_start_pos,
                    ymin=alt_start_ymin,
                    ymax=alt_start_ymax,
                    color=self.feature_colors.get(
                        region_type, self.feature_colors["alternative_start"]
                    ),
                    linewidth=3,
                    label="Alternative Isoform",
                )

                # Place the start codon label and region type
                # FIX: Use the same positioning logic for labels
                if region_type == "extension":
                    label_x = alt_start if transcript_strand == "+" else alt_end
                else:  # truncation
                    label_x = alt_end if transcript_strand == "+" else alt_start

                label_text = alt_feature.get("start_codon", "")
                if region_type:
                    label_text = f"{region_type.capitalize()} {label_text}".strip()
                plt.text(
                    label_x,
                    alt_start_ymax + 0.05,
                    label_text,
                    fontsize=8,
                    rotation=45,
                    ha="center",
                )

        # Plot mutations if provided (only those in zoom region)
        if mutations_df is not None and not mutations_df.empty:
            zoomed_mutations = mutations_df[
                (mutations_df["position"] >= zoom_start)
                & (mutations_df["position"] <= zoom_end)
            ].copy()

            if not zoomed_mutations.empty:
                source_legend_elements = []
                impact_col = (
                    "impact_validated"
                    if "impact_validated" in zoomed_mutations.columns
                    else "impact"
                )
                unique_impacts = mutations_df[impact_col].unique()
                mutation_legend_elements = self._create_mutation_legend_groups(
                    unique_impacts
                )

                # Set jitter parameters
                jitter_amount = 0.02
                np.random.seed(42)

                # Plot mutations by source with jitter
                for source, base_y_pos in self.mutation_track_positions.items():
                    source_data = zoomed_mutations[zoomed_mutations["source"] == source]

                    if not source_data.empty:
                        source_legend_elements.append(
                            Line2D(
                                [0],
                                [0],
                                marker=self.source_markers[source],
                                color="w",
                                markerfacecolor="black",
                                label=source,
                                markersize=6,
                            )
                        )

                        jittered_y = base_y_pos + np.random.uniform(
                            -jitter_amount, jitter_amount, len(source_data)
                        )

                        for (_, mutation), y_pos in zip(
                            source_data.iterrows(), jittered_y
                        ):
                            color = self._get_mutation_color(mutation[impact_col])
                            ax.scatter(
                                mutation["position"],
                                y_pos,
                                marker=self.source_markers[source],
                                c=color,
                                s=8,
                                alpha=0.8,
                            )

        # Create and add legends - DYNAMIC VERSION
        feature_legend_elements = [
            plt.Rectangle(
                (0, 0), 1, 1, facecolor=self.feature_colors["CDS"], label="CDS"
            ),
            plt.Rectangle(
                (0, 0), 1, 1, facecolor=self.feature_colors["UTR"], label="UTR"
            ),
            Line2D(
                [0],
                [0],
                marker="|",
                color=self.feature_colors["start_codon"],
                label="Start codon",
                markersize=15,
                markeredgewidth=3,
                linestyle="None",
            ),
            Line2D(
                [0],
                [0],
                marker="|",
                color=self.feature_colors["stop_codon"],
                label="Stop codon",
                markersize=12,
                markeredgewidth=3,
                linestyle="None",
            ),
        ]

        # Add alternative isoform legend elements based on what's actually present
        if alt_features is not None and not alt_features.empty:
            region_types = alt_features["region_type"].unique()

            for region_type in region_types:
                if region_type == "truncation":
                    feature_legend_elements.extend(
                        [
                            Line2D(
                                [0],
                                [0],
                                marker="|",
                                color=self.feature_colors["truncation"],
                                label="Truncation start",
                                markersize=12,
                                markeredgewidth=3,
                                linestyle="None",
                            ),
                            Line2D(
                                [0],
                                [0],
                                color=self.feature_colors["truncation"],
                                label="Truncation region",
                                linewidth=1.5,
                            ),
                        ]
                    )
                elif region_type == "extension":
                    feature_legend_elements.extend(
                        [
                            Line2D(
                                [0],
                                [0],
                                marker="|",
                                color=self.feature_colors["extension"],
                                label="Extension start",
                                markersize=12,
                                markeredgewidth=3,
                                linestyle="None",
                            ),
                            Line2D(
                                [0],
                                [0],
                                color=self.feature_colors["extension"],
                                label="Extension region",
                                linewidth=1.5,
                            ),
                        ]
                    )

        # Add legends with improved spacing
        if (
            mutations_df is not None
            and not mutations_df.empty
            and not zoomed_mutations.empty
        ):
            l1 = ax.legend(
                handles=feature_legend_elements,
                bbox_to_anchor=(1.05, 1.05),
                loc="upper left",
                borderaxespad=0.0,
                frameon=False,
                title="Transcript Features",
                prop={"size": 8},
                title_fontsize=9,
            )

            ax2 = ax.twinx()
            l2 = ax2.legend(
                handles=source_legend_elements,
                bbox_to_anchor=(1.05, 0.55),
                loc="upper left",
                borderaxespad=0.0,
                frameon=False,
                title="Mutation Sources",
                prop={"size": 8},
                title_fontsize=9,
            )

            ax3 = ax.twinx()
            l3 = ax3.legend(
                handles=mutation_legend_elements,
                bbox_to_anchor=(1.05, 0.3),
                loc="upper left",
                borderaxespad=0.0,
                frameon=False,
                title="Mutation Types",
                prop={"size": 8},
                title_fontsize=9,
            )

            # Hide spines for twin axes
            for axis in [ax2, ax3]:
                axis.set_yticks([])
                for spine in axis.spines.values():
                    spine.set_visible(False)
        else:
            ax.legend(
                handles=feature_legend_elements,
                bbox_to_anchor=(1.05, 1),
                loc="upper left",
                borderaxespad=0.0,
                frameon=False,
            )

        # Customize plot
        plt.title(f"{gene_name} - {transcript_id} (Zoomed)", pad=20, y=1.05)
        ax.set_ylim(0, 1)
        ax.set_xlim(zoom_start, zoom_end)

        # Set tick positions for zoomed region
        tick_interval = 100  # Smaller interval for zoomed view
        base_position = zoom_start - (zoom_start % tick_interval)
        tick_positions = range(
            int(base_position), int(zoom_end) + tick_interval, tick_interval
        )

        plt.xticks(tick_positions, [f"{pos:,}" for pos in tick_positions])
        plt.xticks(rotation=45)

        # Draw the transcript line after setting the ticks
        xlim = ax.get_xlim()
        ax.hlines(
            y=0.4, xmin=xlim[0], xmax=xlim[1], color="black", linewidth=1.5, zorder=0
        )

        # Remove spines
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        plt.tight_layout()

        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(
                output_path,
                bbox_inches="tight",
                dpi=300,
                facecolor="white",
                edgecolor="none",
                format="pdf",
            )
            plt.close()
        else:
            plt.show()
