"""Alternative isoform data handling module.

This module provides functionality for loading and processing alternative
isoform data from BED format files, with support for visualizing
alternative start sites.
"""

import pandas as pd
from typing import Optional, Dict, List


class AlternativeIsoform:
    """Handles alternative isoform data from BED format files.

    The name field in the BED file is expected to contain gene information in the format:
    ENSG00000260916.6_CCPG1_AUG_TruncationToAnno
    """

    def __init__(self):
        """Initialize an empty isoform handler."""
        self.isoforms = pd.DataFrame()

    def load_bed(self, file_path: str) -> None:
        """Load data from BED format file.

        Args:
            file_path: Path to the BED file
        """
        # Read BED format with standard columns
        self.isoforms = pd.read_csv(
            file_path,
            sep="\t",
            names=["chrom", "start", "end", "name", "score", "strand"],
        )

        # Parse the name field which contains gene information
        def parse_name(name: str) -> Dict[str, str]:
            parts = name.split("_")
            if len(parts) >= 4:
                return {
                    "gene_id": parts[0],
                    "gene_name": parts[1],
                    "start_codon": parts[2],
                    "isoform_type": parts[3],
                }
            # Handle malformed names gracefully
            return {
                "gene_id": name,
                "gene_name": name,
                "start_codon": "UNK",
                "isoform_type": "UNK",
            }

        # Parse name field and add as new columns
        info_df = pd.DataFrame(self.isoforms["name"].apply(parse_name).tolist())
        self.isoforms = pd.concat([self.isoforms, info_df], axis=1)

        # Add columns needed for visualization
        self.isoforms["feature_type"] = "alternative_start"
        self.isoforms["source"] = "truncation"

        # Convert position columns to integers
        self.isoforms["start"] = self.isoforms["start"].astype(int)
        self.isoforms["end"] = self.isoforms["end"].astype(int)

    def get_visualization_features(self, gene_name: str) -> pd.DataFrame:
        """Get features formatted for the GenomeVisualizer.

        Args:
            gene_name: Name of the gene to get features for

        Returns:
            DataFrame: Features formatted for visualization

        Raises:
            ValueError: If no data has been loaded
        """
        if self.isoforms.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        features = self.isoforms[self.isoforms["gene_name"] == gene_name].copy()

        if features.empty:
            return pd.DataFrame()

        # Format features for visualization to match GenomeVisualizer expected format
        viz_features = pd.DataFrame(
            {
                "chromosome": features["chrom"],
                "source": "truncation",
                "feature_type": "alternative_start",
                "start": features["start"],
                "end": features["end"],
                "score": features["score"],
                "strand": features["strand"],
                "frame": ".",
                "gene_id": features["gene_id"],
                "transcript_id": features["gene_id"] + "_alt",
                "gene_name": features["gene_name"],
                "start_codon": features["start_codon"],
            }
        )

        return viz_features

    def get_gene_list(self) -> List[str]:
        """Get list of all genes in the dataset.

        Returns:
            List of gene names

        Raises:
            ValueError: If no data has been loaded
        """
        if self.isoforms.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        return sorted(self.isoforms["gene_name"].unique().tolist())

    def get_stats(self) -> Dict:
        """Get basic statistics about the loaded isoforms.

        Returns:
            Dictionary containing various statistics

        Raises:
            ValueError: If no data has been loaded
        """
        if self.isoforms.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        stats = {
            "total_sites": len(self.isoforms),
            "unique_genes": len(self.get_gene_list()),
            "chromosomes": sorted(self.isoforms["chrom"].unique().tolist()),
            "start_codons": sorted(self.isoforms["start_codon"].unique().tolist()),
        }

        return stats

    def get_alternative_starts(self, gene_name: Optional[str] = None) -> pd.DataFrame:
        """Get alternative start sites, optionally filtered by gene name.

        Args:
            gene_name: Gene name to filter by

        Returns:
            DataFrame: Alternative start sites with positions and metadata

        Raises:
            ValueError: If no data has been loaded
        """
        if self.isoforms.empty:
            raise ValueError("No data loaded. Please load data first with load_bed().")

        if gene_name is not None:
            return self.isoforms[self.isoforms["gene_name"] == gene_name]
        return self.isoforms.copy()
