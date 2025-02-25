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
            "CDS": "#2196F3",  # blue
            "UTR": "#FFA500",  # orange
            "start_codon": "#FF0000",  # red
            "stop_codon": "#800080",  # purple
            "truncation": "#FF1493",  # deep pink
            "alternative_start": "#FFD700",  # yellow
        }

        # Mutation impact colors
        self.mutation_colors = {
            "missense variant": "#FF9900",  # orange
            "synonymous variant": "#00CC00",  # green
            "nonsense variant": "#CC0000",  # dark red
            "inframe deletion": "#FF3333",  # light red
            "inframe insertion": "#FF6666",  # lighter red
            "frameshift variant": "#FF0000",  # red
            "splice variant": "#9933FF",  # purple
            "start lost variant": "#990000",  # dark red
            "5 prime utr variant": "#FFB366",  # light orange
            "3 prime utr variant": "#FFB366",  # light orange
            "other variant": "#999999",  # gray
            "unknown": "#000000",  # black
        }

        # Source-specific markers
        self.source_markers = {
            "gnomAD": "o",  # circle
            "ClinVar": "^",  # triangle
            "Aggregator": "s",  # square
        }

        # Track dimensions
        self.track_height = 0.15
        self.codon_height = 0.2

        # Y-positions for mutation tracks
        self.mutation_track_positions = {
            "gnomAD": 0.7,
            "ClinVar": 0.8,
            "Aggregator": 0.9,
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

        if features.empty:
            raise ValueError(f"No features found for transcript {transcript_id}")

        # Get transcript bounds and create figure
        transcript_info = features[features["feature_type"] == "transcript"].iloc[0]
        transcript_start = transcript_info["start"]
        transcript_end = transcript_info["end"]
        span = transcript_end - transcript_start

        fig, ax = plt.subplots(figsize=(15, 4))

        # Draw base transcript line
        ax.hlines(
            y=0.4,
            xmin=transcript_start,
            xmax=transcript_end,
            color="black",
            linewidth=1.5,
        )

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

            elif feature["feature_type"] in ["5UTR", "3UTR"]:
                rect = Rectangle(
                    (feature["start"], 0.325),
                    width,
                    self.track_height,
                    facecolor=self.feature_colors["UTR"],
                )
                ax.add_patch(rect)

            elif feature["feature_type"] in ["start_codon", "stop_codon"]:
                rect = Rectangle(
                    (feature["start"], 0.3),
                    width,
                    self.codon_height,
                    facecolor=self.feature_colors[feature["feature_type"]],
                )
                ax.add_patch(rect)

        # Plot alternative start sites
        if alt_features is not None and not alt_features.empty:
            for _, alt_feature in alt_features.iterrows():
                alt_start = alt_feature["start"]
                alt_end = alt_feature["end"]

                # Draw truncation bracket
                bracket_height = 0.05
                bracket_y = 0.2

                ax.hlines(
                    y=bracket_y,
                    xmin=alt_start - 0.5,
                    xmax=alt_end + 0.5,
                    color=self.feature_colors["truncation"],
                    linewidth=1.5,
                    label="Truncation",
                )

                ax.vlines(
                    x=[alt_start, alt_end - 0.5],
                    ymin=bracket_y,
                    ymax=bracket_y + bracket_height,
                    color=self.feature_colors["truncation"],
                    linewidth=0.5,
                )

                # Alternative start marker
                black_bar_y = 0.4
                alt_start_height = self.codon_height
                alt_start_ymin = black_bar_y - (alt_start_height / 2)
                alt_start_ymax = black_bar_y + (alt_start_height / 2)

                ax.vlines(
                    x=alt_end,
                    ymin=alt_start_ymin,
                    ymax=alt_start_ymax,
                    color=self.feature_colors["alternative_start"],
                    linewidth=2,
                    label="Alternative start",
                )

                if "start_codon" in alt_feature:
                    plt.text(
                        alt_end - span * 0.01,
                        alt_start_ymax + 0.05,
                        alt_feature["start_codon"],
                        fontsize=8,
                        rotation=45,
                    )

        # Plot mutations if provided
        if mutations_df is not None and not mutations_df.empty:
            # Create legend elements
            source_legend_elements = []
            unique_impacts = mutations_df["impact"].unique()
            mutation_legend_elements = self._create_mutation_legend_groups(
                unique_impacts
            )

            # Set jitter parameters
            jitter_amount = 0.02  # Adjust this to control spread
            np.random.seed(42)  # For reproducibility

            # Plot mutations by source with jitter
            for source, base_y_pos in self.mutation_track_positions.items():
                source_data = mutations_df[mutations_df["source"] == source]

                if not source_data.empty:
                    # Add source to legend
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
                    )  # Smaller legend markers

                    # Generate jittered y-positions for this source
                    jittered_y = base_y_pos + np.random.uniform(
                        -jitter_amount, jitter_amount, len(source_data)
                    )

                    # Plot mutations with jitter
                    for (_, mutation), y_pos in zip(source_data.iterrows(), jittered_y):
                        color = self._get_mutation_color(mutation["impact"])
                        ax.scatter(
                            mutation["position"],
                            y_pos,
                            marker=self.source_markers[source],
                            c=color,
                            s=8,  # Smaller dots
                            alpha=0.8,
                        )  # Slight transparency

        # Create legends
        feature_legend_elements = [
            plt.Rectangle(
                (0, 0), 1, 1, facecolor=self.feature_colors["CDS"], label="CDS"
            ),
            plt.Rectangle(
                (0, 0), 1, 1, facecolor=self.feature_colors["UTR"], label="UTR"
            ),
            plt.Rectangle(
                (0, 0),
                1,
                1,
                facecolor=self.feature_colors["start_codon"],
                label="Start codon",
            ),
            plt.Rectangle(
                (0, 0),
                1,
                1,
                facecolor=self.feature_colors["stop_codon"],
                label="Stop codon",
            ),
        ]

        if alt_features is not None and not alt_features.empty:
            feature_legend_elements.extend(
                [
                    Line2D(
                        [0],
                        [0],
                        color=self.feature_colors["truncation"],
                        label="Truncation",
                        linewidth=1.5,
                    ),
                    Line2D(
                        [0],
                        [0],
                        color=self.feature_colors["alternative_start"],
                        label="Alternative start",
                        linewidth=2,
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

        # Calculate zoom region based on alt_features
        zoom_start = alt_features["start"].min() - padding
        zoom_end = alt_features["end"].max() + padding

        # Create figure
        fig, ax = plt.subplots(figsize=(15, 4))

        # Draw base transcript line
        ax.hlines(y=0.4, xmin=zoom_start, xmax=zoom_end, color="black", linewidth=1.5)

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

            elif feature["feature_type"] in ["5UTR", "3UTR"]:
                rect = Rectangle(
                    (feature["start"], 0.325),
                    width,
                    self.track_height,
                    facecolor=self.feature_colors["UTR"],
                )
                ax.add_patch(rect)

            elif feature["feature_type"] in ["start_codon", "stop_codon"]:
                rect = Rectangle(
                    (feature["start"], 0.3),
                    width,
                    self.codon_height,
                    facecolor=self.feature_colors[feature["feature_type"]],
                )
                ax.add_patch(rect)

        # Plot alternative start sites
        for _, alt_feature in alt_features.iterrows():
            alt_start = alt_feature["start"]
            alt_end = alt_feature["end"]

            # Draw truncation bracket
            bracket_height = 0.05
            bracket_y = 0.2

            ax.hlines(
                y=bracket_y,
                xmin=alt_start,
                xmax=alt_end,
                color=self.feature_colors["truncation"],
                linewidth=1.5,
                label="Truncation",
            )

            ax.vlines(
                x=[alt_start, alt_end],
                ymin=bracket_y,
                ymax=bracket_y + bracket_height,
                color=self.feature_colors["truncation"],
                linewidth=1,
            )

            # Alternative start marker
            black_bar_y = 0.4
            alt_start_height = self.codon_height
            alt_start_ymin = black_bar_y - (alt_start_height / 2)
            alt_start_ymax = black_bar_y + (alt_start_height / 2)

            ax.vlines(
                x=alt_end,
                ymin=alt_start_ymin,
                ymax=alt_start_ymax,
                color=self.feature_colors["alternative_start"],
                linewidth=2,
                label="Alternative start",
            )

            if "start_codon" in alt_feature:
                plt.text(
                    alt_end - 5,
                    alt_start_ymax + 0.05,
                    alt_feature["start_codon"],
                    fontsize=8,
                    rotation=45,
                )

        # Plot mutations if provided (only those in zoom region)
        if mutations_df is not None and not mutations_df.empty:
            zoomed_mutations = mutations_df[
                (mutations_df["position"] >= zoom_start)
                & (mutations_df["position"] <= zoom_end)
            ].copy()

            if not zoomed_mutations.empty:
                source_legend_elements = []
                unique_impacts = zoomed_mutations["impact"].unique()
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
                            color = self._get_mutation_color(mutation["impact"])
                            ax.scatter(
                                mutation["position"],
                                y_pos,
                                marker=self.source_markers[source],
                                c=color,
                                s=8,
                                alpha=0.8,
                            )

        # Create and add legends
        feature_legend_elements = [
            plt.Rectangle(
                (0, 0), 1, 1, facecolor=self.feature_colors["CDS"], label="CDS"
            ),
            plt.Rectangle(
                (0, 0), 1, 1, facecolor=self.feature_colors["UTR"], label="UTR"
            ),
            plt.Rectangle(
                (0, 0),
                1,
                1,
                facecolor=self.feature_colors["start_codon"],
                label="Start codon",
            ),
            plt.Rectangle(
                (0, 0),
                1,
                1,
                facecolor=self.feature_colors["stop_codon"],
                label="Stop codon",
            ),
            Line2D(
                [0],
                [0],
                color=self.feature_colors["truncation"],
                label="Truncation",
                linewidth=1.5,
            ),
            Line2D(
                [0],
                [0],
                color=self.feature_colors["alternative_start"],
                label="Alternative start",
                linewidth=2,
            ),
        ]

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
            )
            plt.close()
        else:
            plt.show()
