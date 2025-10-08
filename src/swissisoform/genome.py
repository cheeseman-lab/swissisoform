"""Genome sequence and annotation handling module.

This module provides the GenomeHandler class for working with genome sequences
and annotation files. It supports loading and querying FASTA genome sequences
and GTF annotation files.
"""

from Bio import SeqIO
import pandas as pd
import logging
from typing import Dict, Optional, Tuple

logger = logging.getLogger(__name__)


class GenomeHandler:
    """Handles access to genome sequence data and feature annotations.

    Attributes:
        genome (Dict[str, SeqRecord]): Dictionary of sequence records keyed by chromosome.
        annotations (pd.DataFrame): DataFrame containing genomic feature annotations.
        gtf_path (Optional[str]): Path to the genome annotation GTF file.
    """

    def __init__(self, genome_path: str, gtf_path: Optional[str] = None):
        """Initialize the GenomeHandler with genome sequence and annotations.

        Args:
            genome_path (str): Path to the genome FASTA file.
            gtf_path (Optional[str]): Path to the genome annotation GTF file.

        Returns:
            None
        """
        self.genome = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
        self.gtf_path = gtf_path
        if gtf_path:
            self.load_annotations(gtf_path)

    def load_annotations(self, gtf_path: str) -> None:
        """Load and parse GTF annotations into a pandas DataFrame.

        Args:
            gtf_path (str): Path to the genome annotation GTF file.

        Returns:
            None
        """
        features_list = []

        with open(gtf_path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")

                if len(fields) != 9:
                    continue

                # Parse attributes - GENCODE format is different from NCBI RefSeq
                attrs = {}
                attribute_string = fields[8]
                attr_pairs = [x.strip() for x in attribute_string.split(";")]

                for attr in attr_pairs:
                    if not attr:
                        continue

                    # GENCODE format: key "value"
                    match = attr.strip().split(" ", 1)
                    if len(match) == 2:
                        key = match[0].strip()
                        # Remove quotes from value
                        value = match[1].strip().strip('"')
                        attrs[key] = value

                feature_dict = {
                    "chromosome": fields[0],
                    "source": fields[1],
                    "feature_type": fields[2],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "score": fields[5],
                    "strand": fields[6],
                    "frame": fields[7],
                    **attrs,
                }

                features_list.append(feature_dict)

        # Convert to DataFrame
        df = pd.DataFrame(features_list)

        # Reorder columns
        base_cols = [
            "chromosome",
            "source",
            "feature_type",
            "start",
            "end",
            "score",
            "strand",
            "frame",
        ]
        other_cols = [col for col in df.columns if col not in base_cols]
        self.annotations = df[base_cols + other_cols]

    def find_gene_features(self, gene_name: str) -> pd.DataFrame:
        """Find all features associated with a gene name.

        Args:
            gene_name (str): Name of the gene.

        Returns:
            pd.DataFrame: DataFrame containing features associated with the gene.

        Raises:
            ValueError: If no GTF file was loaded.
        """
        if not hasattr(self, "annotations"):
            raise ValueError("No GTF file loaded")

        return self.annotations[self.annotations["gene_name"] == gene_name]

    def get_sequence(self, chrom: str, start: int, end: int, strand: str = "+") -> str:
        """Get sequence for a genomic region.

        Args:
            chrom (str): Chromosome name.
            start (int): Start position (1-based).
            end (int): End position (inclusive).
            strand (str): '+' for forward strand, '-' for reverse.

        Returns:
            str: Nucleotide sequence for the specified region.

        Raises:
            ValueError: If chromosome not found in genome.
        """
        # Try to find chromosome in the genome, checking different prefix formats
        if chrom in self.genome:
            seq_id = chrom
        elif chrom.startswith("chr") and chrom[3:] in self.genome:
            seq_id = chrom[3:]  # Try without "chr" prefix
        elif not chrom.startswith("chr") and f"chr{chrom}" in self.genome:
            seq_id = f"chr{chrom}"  # Try with "chr" prefix
        else:
            raise ValueError(f"Chromosome {chrom} not found in genome")

        seq = self.genome[seq_id].seq[start - 1 : end]
        return seq if strand == "+" else seq.reverse_complement()

    def get_transcript_features(self, transcript_id: str) -> pd.DataFrame:
        """Get all features associated with a transcript ID.

        Args:
            transcript_id (str): Transcript ID.

        Returns:
            pd.DataFrame: DataFrame containing features associated with the transcript.

        Raises:
            ValueError: If no GTF file was loaded.
        """
        if not hasattr(self, "annotations"):
            raise ValueError("No GTF file loaded")

        return self.annotations[self.annotations["transcript_id"] == transcript_id]

    def get_gene_stats(self, gene_name: str) -> Optional[Dict]:
        """Get basic statistics for a gene.

        Args:
            gene_name (str): Name of the gene.

        Returns:
            Optional[Dict]: Dictionary of gene statistics or None if gene not found.
        """
        features = self.find_gene_features(gene_name)
        if features.empty:
            return None

        stats = {
            "transcripts": len(features[features["feature_type"] == "transcript"]),
            "exons": len(features[features["feature_type"] == "exon"]),
            "cds": len(features[features["feature_type"] == "CDS"]),
            "utrs": len(features[features["feature_type"].isin(["5UTR", "3UTR"])]),
            "chromosomes": features["chromosome"].unique().tolist(),
            "total_span": features["end"].max() - features["start"].min() + 1,
        }

        return stats

    def get_transcript_ids(
        self, gene_name: str, standard_chroms_only: bool = True
    ) -> pd.DataFrame:
        """Get all transcript IDs for a given gene name.

        Args:
            gene_name (str): Name of the gene.
            standard_chroms_only (bool): If True, only return transcripts on standard chromosomes
                (chr1-22, chrX, chrY, chrM).

        Returns:
            pd.DataFrame: DataFrame containing transcript information.

        Raises:
            ValueError: If no GTF file was loaded.
        """
        if not hasattr(self, "annotations"):
            raise ValueError("No GTF file loaded")

        # Get all transcripts first
        transcripts = self.annotations[
            (self.annotations["gene_name"] == gene_name)
            & (self.annotations["feature_type"] == "transcript")
        ]

        # Basic transcript info
        transcript_info = transcripts[
            ["transcript_id", "chromosome", "start", "end", "strand"]
        ].drop_duplicates()

        if standard_chroms_only:
            # Define standard chromosomes for both formats (with and without "chr" prefix)
            standard_chromosomes = set()
            for i in range(1, 23):
                standard_chromosomes.add(f"chr{i}")
                standard_chromosomes.add(str(i))
            standard_chromosomes.update(["chrX", "chrY", "chrM", "X", "Y", "MT", "M"])

            # Filter and count
            original_count = len(transcript_info)
            transcript_info = transcript_info[
                transcript_info["chromosome"].isin(standard_chromosomes)
            ]
            filtered_count = len(transcript_info)

            # Print filtering info if any transcripts were filtered
            if filtered_count < original_count:
                logger.info(
                    f"Note: Filtered out {original_count - filtered_count} transcript(s) on non-standard chromosomes"
                )
                logger.info(
                    "Non-standard chromosomes:",
                    set(transcripts["chromosome"].unique()) - standard_chromosomes,
                )

        return transcript_info

    def get_transcript_sequence(self, transcript_id: str) -> Dict:
        """Get the full genomic sequence for a transcript.

        Args:
            transcript_id (str): Transcript ID.

        Returns:
            Dict: Dictionary containing transcript sequence and associated metadata.

        Raises:
            ValueError: If no GTF file was loaded.
            IndexError: If transcript not found.
        """
        if not hasattr(self, "annotations"):
            raise ValueError("No GTF file loaded")

        # Get transcript features
        transcript = self.annotations[
            (self.annotations["transcript_id"] == transcript_id)
            & (self.annotations["feature_type"] == "transcript")
        ]

        if transcript.empty:
            raise IndexError(f"Transcript {transcript_id} not found")

        transcript = transcript.iloc[0]

        # Get the sequence
        sequence = self.get_sequence(
            transcript["chromosome"],
            transcript["start"],
            transcript["end"],
            transcript["strand"],
        )

        return {
            "transcript_id": transcript_id,
            "chromosome": transcript["chromosome"],
            "start": transcript["start"],
            "end": transcript["end"],
            "strand": transcript["strand"],
            "sequence": str(sequence),
        }

    def get_transcript_features_with_sequence(self, transcript_id: str) -> Dict:
        """Get all features and sequence for a transcript.

        Args:
            transcript_id (str): Transcript ID.

        Returns:
            Dict: Dictionary containing transcript features and sequence information.
        """
        # Get transcript features
        features = self.get_transcript_features(transcript_id)
        sequence_info = self.get_transcript_sequence(transcript_id)

        return {"features": features, "sequence": sequence_info}
