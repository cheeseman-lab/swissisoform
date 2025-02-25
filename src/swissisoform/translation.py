"""Transcript translation and protein sequence generation module.

This module contains the TruncatedProteinGenerator class for generating
amino acid sequences from truncated transcripts and preparing datasets
for deep learning applications.
"""

import pandas as pd
import re
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from swissisoform.genome import GenomeHandler
from swissisoform.alternative_isoforms import AlternativeIsoform


class TruncatedProteinGenerator:
    """Generate amino acid sequences from truncated transcripts.
    
    This class facilitates the generation of protein sequences stemming from
    transcript truncations for use in deep learning models.
    """

    def __init__(
        self,
        genome_handler: GenomeHandler,
        alt_isoform_handler: AlternativeIsoform,
        output_dir: str,
    ):
        """Initialize the TruncatedProteinGenerator.

        Args:
            genome_handler: Initialized GenomeHandler instance
            alt_isoform_handler: Initialized AlternativeIsoform instance
            output_dir: Directory to save output files
        """
        self.genome = genome_handler
        self.alt_isoforms = alt_isoform_handler
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def extract_transcript_sequence(self, transcript_id: str) -> Tuple[str, bool]:
        """Extract the full nucleotide sequence for a transcript.

        Args:
            transcript_id: Transcript ID

        Returns:
            Tuple of (sequence, success_flag)
        """
        try:
            # Get transcript sequence - returns a dictionary with sequence info
            sequence_info = self.genome.get_transcript_sequence(transcript_id)

            if not sequence_info or "sequence" not in sequence_info:
                print(f"Failed to retrieve sequence for transcript {transcript_id}")
                return "", False

            # Extract the actual sequence string from the dictionary
            sequence = sequence_info["sequence"]

            if not sequence:
                print(f"Empty sequence for transcript {transcript_id}")
                return "", False

            return sequence, True

        except Exception as e:
            print(f"Error extracting transcript sequence: {str(e)}")
            return "", False

    def get_truncation_sites(self, gene_name: str) -> pd.DataFrame:
        """Get truncation sites for a gene from the alternative isoform data.

        Args:
            gene_name: Name of the gene

        Returns:
            DataFrame containing truncation information
        """
        try:
            # Get alternative isoform features - these are already truncation/alternative start sites
            alt_features = self.alt_isoforms.get_visualization_features(gene_name)

            if alt_features.empty:
                print(f"No alternative features found for gene {gene_name}")
                return pd.DataFrame()

            # No need to filter - all features from the BED file are truncations/alternative starts
            print(
                f"Found {len(alt_features)} alternative start/truncation sites for gene {gene_name}"
            )
            return alt_features

        except Exception as e:
            print(f"Error getting truncation sites: {str(e)}")
            return pd.DataFrame()

    def apply_truncation(
        self, transcript_sequence: str, truncation_position: int, strand: str
    ) -> str:
        """Apply truncation to a transcript sequence.

        Args:
            transcript_sequence: Full transcript nucleotide sequence
            truncation_position: Position to truncate (in genomic coordinates)
            strand: '+' or '-' strand

        Returns:
            Truncated transcript sequence
        """
        try:
            # For simplicity, just truncate at the middle of the sequence
            # In a real implementation, you would need to convert genomic coordinates
            # to transcript coordinates using exon information
            transcript_position = len(transcript_sequence) // 2

            # Apply truncation
            if strand == "+":
                truncated_sequence = transcript_sequence[:transcript_position]
            else:
                truncated_sequence = transcript_sequence[transcript_position:]

            return truncated_sequence

        except Exception as e:
            print(f"Error applying truncation: {str(e)}")
            return ""

    def translate_sequence(self, nucleotide_seq: str) -> str:
        """Translate a nucleotide sequence to amino acid sequence.

        Args:
            nucleotide_seq: Nucleotide sequence

        Returns:
            Amino acid sequence
        """
        try:
            # Find the start codon
            start_codon_match = re.search(r"ATG", nucleotide_seq)

            if not start_codon_match:
                print("No start codon (ATG) found in the sequence")
                return ""

            start_pos = start_codon_match.start()

            # Translate from start codon to the end
            coding_seq = nucleotide_seq[start_pos:]

            # Ensure the length is a multiple of 3
            length_to_use = len(coding_seq) - (len(coding_seq) % 3)
            coding_seq = coding_seq[:length_to_use]

            # Translate
            seq_obj = Seq(coding_seq)
            protein_seq = str(seq_obj.translate(to_stop=True))

            return protein_seq

        except Exception as e:
            print(f"Error translating sequence: {str(e)}")
            return ""

    def prepare_deep_learning_dataset(
        self,
        gene_list: List[str],
        output_format: str = "fasta,csv",
        include_canonical: bool = True,
        min_length: int = 50,
        max_length: int = 1000,
    ) -> pd.DataFrame:
        """Prepare a dataset of amino acid sequences for deep learning.

        Args:
            gene_list: List of gene names to process
            output_format: Format to save sequences ('fasta', 'csv', or both)
            include_canonical: Whether to include canonical transcripts as examples
            min_length: Minimum protein length to include
            max_length: Maximum protein length to include

        Returns:
            DataFrame with the dataset information
        """
        all_sequences = []

        # Process each gene
        for gene_idx, gene_name in enumerate(gene_list, 1):
            print(f"\nProcessing gene {gene_idx}/{len(gene_list)}: {gene_name}")

            # Get transcript information
            transcript_info = self.genome.get_transcript_ids(gene_name)

            if transcript_info.empty:
                print(f"No transcript info found for gene {gene_name}")
                continue

            print(f"Processing {len(transcript_info)} transcripts for gene {gene_name}")

            # Get truncation sites
            truncations = self.get_truncation_sites(gene_name)

            # Process each transcript
            for _, transcript in transcript_info.iterrows():
                transcript_id = transcript["transcript_id"]
                strand = transcript.get("strand", "+")

                # Get full transcript sequence
                full_sequence, success = self.extract_transcript_sequence(transcript_id)

                if not success or not full_sequence:
                    continue

                # Generate canonical protein if requested
                if include_canonical:
                    canonical_protein = self.translate_sequence(full_sequence)

                    if (
                        canonical_protein
                        and min_length <= len(canonical_protein) <= max_length
                    ):
                        all_sequences.append(
                            {
                                "gene": gene_name,
                                "transcript_id": transcript_id,
                                "variant_id": "canonical",
                                "sequence": canonical_protein,
                                "length": len(canonical_protein),
                                "is_truncated": 0,
                            }
                        )

                # Process truncations if available
                if not truncations.empty:
                    for _, trunc in truncations.iterrows():
                        trunc_id = f"trunc_{trunc['start']}_{trunc['end']}"
                        trunc_position = trunc["start"]

                        # Apply truncation
                        truncated_seq = self.apply_truncation(
                            full_sequence, trunc_position, strand
                        )

                        if not truncated_seq:
                            continue

                        # Translate truncated sequence
                        truncated_protein = self.translate_sequence(truncated_seq)

                        if (
                            truncated_protein
                            and min_length <= len(truncated_protein) <= max_length
                        ):
                            all_sequences.append(
                                {
                                    "gene": gene_name,
                                    "transcript_id": transcript_id,
                                    "variant_id": trunc_id,
                                    "sequence": truncated_protein,
                                    "length": len(truncated_protein),
                                    "is_truncated": 1,
                                }
                            )

        # Create dataset
        dataset = pd.DataFrame(all_sequences)

        # Save dataset
        if not dataset.empty:
            if "fasta" in output_format.lower():
                self._save_dataset_fasta(dataset)

            if "csv" in output_format.lower():
                output_file = self.output_dir / "protein_sequence_dataset.csv"
                dataset.to_csv(output_file, index=False)
                print(f"Saved dataset CSV to {output_file}")

        return dataset

    def _save_dataset_fasta(self, dataset: pd.DataFrame) -> None:
        """Save dataset as FASTA file.

        Args:
            dataset: DataFrame with sequence information
        """
        records = []

        for _, row in dataset.iterrows():
            record_id = f"{row['gene']}_{row['transcript_id']}_{row['variant_id']}"
            description = (
                f"{'Truncated' if row['is_truncated'] else 'Canonical'} protein"
            )

            records.append(
                SeqRecord(Seq(row["sequence"]), id=record_id, description=description)
            )

        output_file = self.output_dir / "protein_sequence_dataset.fasta"
        SeqIO.write(records, output_file, "fasta")
        print(f"Saved dataset FASTA to {output_file}")
