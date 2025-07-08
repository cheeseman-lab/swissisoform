"""Transcript translation and protein sequence generation module.

This module contains the TruncatedProteinGenerator class for generating
amino acid sequences from truncated transcripts and preparing datasets
for protein sequence analysis.
"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional

import pandas as pd
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from swissisoform.genome import GenomeHandler
from swissisoform.alternative_isoforms import AlternativeIsoform


class TruncatedProteinGenerator:
    """Generate amino acid sequences from truncated transcripts.

    This class facilitates the generation of protein sequences stemming from
    transcript truncations for downstream analysis and modeling.
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

    def get_truncation_sites(self, gene_name: str) -> pd.DataFrame:
        """Get truncation sites for a gene from the alternative isoform data.

        Args:
            gene_name: Name of the gene

        Returns:
            DataFrame containing truncation information
        """
        try:
            # Get alternative isoform features
            alt_features = self.alt_isoforms.get_visualization_features(gene_name)

            if alt_features.empty:
                print(f"No alternative features found for gene {gene_name}")
                return pd.DataFrame()

            print(
                f"Found {len(alt_features)} alternative start/truncation sites for gene {gene_name}"
            )
            return alt_features

        except Exception as e:
            print(f"Error getting truncation sites: {str(e)}")
            return pd.DataFrame()

    def apply_truncation(
        self, transcript_id: str, truncation_feature: pd.Series
    ) -> str:
        """Apply truncation to a transcript by removing the truncated segment and joining the sequence before and after.

        Args:
            transcript_id: Transcript ID to truncate
            truncation_feature: Feature containing truncation information

        Returns:
            Truncated transcript sequence
        """
        try:
            # Get transcript data
            transcript_data = self.genome.get_transcript_features_with_sequence(
                transcript_id
            )
            if not transcript_data or "sequence" not in transcript_data.get(
                "sequence", {}
            ):
                print(f"Failed to retrieve sequence for {transcript_id}")
                return ""

            # Extract necessary information
            full_sequence = transcript_data["sequence"]["sequence"]
            transcript_start = transcript_data["sequence"]["start"]
            strand = transcript_data["sequence"]["strand"]

            # Get truncation information
            trunc_start = truncation_feature["start"]
            trunc_end = truncation_feature["end"]

            # Convert to transcript coordinates
            rel_trunc_start = trunc_start - transcript_start
            rel_trunc_end = (
                trunc_end - transcript_start + 1
            )  # +1 because end is inclusive

            if strand == "-":
                # For negative strand, we need to reverse the coordinates
                rel_trunc_start, rel_trunc_end = (
                    len(full_sequence) - rel_trunc_end,
                    len(full_sequence) - rel_trunc_start,
                )

            # Ensure the truncation coordinates are within the transcript bounds
            rel_trunc_start = max(0, min(rel_trunc_start, len(full_sequence)))
            rel_trunc_end = max(0, min(rel_trunc_end, len(full_sequence)))

            # Extract segments before and after truncation
            before_segment = full_sequence[:rel_trunc_start]
            after_segment = full_sequence[rel_trunc_end:]

            # Join segments to create truncated sequence
            truncated_sequence = before_segment + after_segment

            return truncated_sequence

        except Exception as e:
            print(f"Error applying truncation: {str(e)}")
            return ""

    def build_truncated_cds(
        self, transcript_id: str, truncation_feature: pd.Series
    ) -> str:
        """Build a properly truncated CDS by analyzing how the truncation affects each CDS region.

        This provides better reading frame preservation than simply truncating the transcript.

        Args:
            transcript_id: Transcript ID
            truncation_feature: Feature containing truncation information

        Returns:
            Truncated CDS sequence
        """
        try:
            # Get transcript data
            transcript_data = self.genome.get_transcript_features_with_sequence(
                transcript_id
            )
            if not transcript_data or "sequence" not in transcript_data.get(
                "sequence", {}
            ):
                return ""

            # Extract necessary information
            features = transcript_data["features"]
            full_sequence = transcript_data["sequence"]["sequence"]
            transcript_start = transcript_data["sequence"]["start"]
            strand = transcript_data["sequence"]["strand"]

            # Get CDS regions
            cds_regions = features[features["feature_type"] == "CDS"].copy()
            if cds_regions.empty:
                return ""

            cds_regions.sort_values("start", inplace=True)

            # For mock truncation (to get canonical sequence), return full CDS
            if truncation_feature["start"] < 0 or truncation_feature["end"] < 0:
                canonical_cds = ""
                for _, cds in cds_regions.iterrows():
                    rel_start = cds["start"] - transcript_start
                    rel_end = cds["end"] - transcript_start + 1

                    if strand == "-":
                        rel_start = len(full_sequence) - rel_end
                        rel_end = len(full_sequence) - (cds["start"] - transcript_start)

                    cds_seq = full_sequence[rel_start:rel_end]
                    canonical_cds += cds_seq

                # Ensure CDS length is divisible by 3
                remainder = len(canonical_cds) % 3
                if remainder > 0:
                    canonical_cds = canonical_cds[: len(canonical_cds) - remainder]

                return canonical_cds

            # Get truncation information
            trunc_start = truncation_feature["start"]
            trunc_end = truncation_feature["end"]

            # Convert to transcript coordinates
            rel_trunc_start = trunc_start - transcript_start
            rel_trunc_end = trunc_end - transcript_start + 1

            if strand == "-":
                rel_trunc_start, rel_trunc_end = (
                    len(full_sequence) - rel_trunc_end,
                    len(full_sequence) - rel_trunc_start,
                )

            # Build truncated CDS by examining each CDS region
            truncated_cds = ""

            for _, cds in cds_regions.iterrows():
                # Calculate relative CDS positions
                rel_cds_start = cds["start"] - transcript_start
                rel_cds_end = cds["end"] - transcript_start + 1

                if strand == "-":
                    rel_cds_start = len(full_sequence) - rel_cds_end
                    rel_cds_end = len(full_sequence) - (cds["start"] - transcript_start)

                # CDS region is entirely before the truncation
                if rel_cds_end <= rel_trunc_start:
                    cds_seq = full_sequence[rel_cds_start:rel_cds_end]
                    truncated_cds += cds_seq

                # CDS region is entirely after the truncation
                elif rel_cds_start >= rel_trunc_end:
                    cds_seq = full_sequence[rel_cds_start:rel_cds_end]
                    truncated_cds += cds_seq

                # CDS region overlaps with the truncation
                else:
                    # Keep the part before truncation if any
                    if rel_cds_start < rel_trunc_start:
                        before_part = full_sequence[rel_cds_start:rel_trunc_start]
                        truncated_cds += before_part

                    # Keep the part after truncation if any
                    if rel_cds_end > rel_trunc_end:
                        after_part = full_sequence[rel_trunc_end:rel_cds_end]
                        truncated_cds += after_part

            # Ensure CDS length is divisible by 3
            remainder = len(truncated_cds) % 3
            if remainder > 0:
                truncated_cds = truncated_cds[: len(truncated_cds) - remainder]

            return truncated_cds

        except Exception as e:
            print(f"Error building truncated CDS: {str(e)}")
            return ""

    def translate_sequence(self, nucleotide_seq: str) -> str:
        """Translate a nucleotide sequence to amino acid sequence.

        Args:
            nucleotide_seq: Nucleotide sequence

        Returns:
            Amino acid sequence
        """
        try:
            if not nucleotide_seq or len(nucleotide_seq) < 3:
                return ""

            # Ensure sequence length is divisible by 3
            remainder = len(nucleotide_seq) % 3
            if remainder > 0:
                nucleotide_seq = nucleotide_seq[: len(nucleotide_seq) - remainder]

            # Translate directly (without requiring start codon)
            seq_obj = Seq(nucleotide_seq)
            protein_seq = str(seq_obj.translate())

            return protein_seq

        except Exception as e:
            print(f"Error translating sequence: {str(e)}")
            return ""

    def generate_protein_variants(
        self,
        gene_name: str,
        output_format: str = "fasta",
        exclude_canonical: bool = False,
    ) -> Dict[str, Dict[str, str]]:
        """Generate amino acid sequences for all truncation variants of a gene.

        Args:
            gene_name: Name of the gene
            output_format: Format to save sequences ('fasta', 'csv', or None for no save)
            exclude_canonical: Whether to exclude canonical transcripts

        Returns:
            Dictionary with transcript IDs as keys, containing both canonical and
            truncated protein sequences for each transcript
        """
        result = {}

        try:
            # Get transcript information
            transcript_info = self.genome.get_transcript_ids(gene_name)

            if transcript_info.empty:
                print(f"No transcript info found for gene {gene_name}")
                return result

            print(f"Processing {len(transcript_info)} transcripts for gene {gene_name}")

            # Get truncation sites
            truncations = self.get_truncation_sites(gene_name)

            # Process each transcript
            for _, transcript in transcript_info.iterrows():
                transcript_id = transcript["transcript_id"]

                print(f"  Processing transcript {transcript_id}")

                # Generate canonical protein first
                canonical_cds = self.build_truncated_cds(
                    transcript_id,
                    pd.Series(
                        {
                            "start": -1,  # Set to impossible values to get full CDS
                            "end": -1,
                        }
                    ),
                )

                if not canonical_cds:
                    print(f"  Failed to get CDS for {transcript_id}, skipping")
                    continue

                canonical_protein = self.translate_sequence(canonical_cds)

                if not canonical_protein:
                    print(
                        f"  Failed to translate canonical sequence for {transcript_id}, skipping"
                    )
                    continue

                # Store canonical sequence
                result[transcript_id] = {"canonical": canonical_protein}

                # Process truncations if available
                if not truncations.empty:
                    for _, trunc in truncations.iterrows():
                        trunc_id = trunc.get(
                            "feature_id", f"trunc_{trunc['start']}_{trunc['end']}"
                        )

                        # Build truncated CDS
                        truncated_cds = self.build_truncated_cds(transcript_id, trunc)

                        if not truncated_cds:
                            print(
                                f"  Failed to build truncated CDS for {transcript_id} - {trunc_id}, skipping"
                            )
                            continue

                        # Translate truncated sequence
                        truncated_protein = self.translate_sequence(truncated_cds)

                        if not truncated_protein:
                            print(
                                f"  Failed to translate truncated sequence for {transcript_id} - {trunc_id}, skipping"
                            )
                            continue

                        # Store truncated sequence
                        result[transcript_id][trunc_id] = truncated_protein
                        print(
                            f"  Generated truncated protein for {transcript_id} - {trunc_id}"
                        )

            # Save results if requested
            if output_format and result:
                self._save_sequences(
                    gene_name, result, output_format, exclude_canonical
                )

            return result

        except Exception as e:
            print(f"Error generating sequences for gene {gene_name}: {str(e)}")
            return result

    def create_protein_sequence_dataset_pairs(
        self,
        gene_list: List[str],
        output_format: str = "fasta,csv",
        min_length: int = 50,
        max_length: int = 1000,
    ) -> pd.DataFrame:
        """Create a dataset of protein sequences from paired canonical and truncated transcripts.

        This function processes multiple genes and creates a dataset containing only
        canonical transcripts and their truncated versions where the truncation
        affects the transcript (overlaps with coding region).

        Args:
            gene_list: List of gene names to process
            output_format: Format to save sequences ('fasta', 'csv', or both)
            min_length: Minimum protein length to include
            max_length: Maximum protein length to include

        Returns:
            DataFrame with the dataset information
        """
        all_sequences = []
        total_pairs = 0
        skipped_genes = []

        # Process each gene
        for gene_idx, gene_name in enumerate(gene_list, 1):
            print(f"\nProcessing gene {gene_idx}/{len(gene_list)}: {gene_name}")

            try:
                # Get transcript information
                transcript_info = self.genome.get_transcript_ids(gene_name)
                if transcript_info.empty:
                    print(
                        f"  ├─ No transcript info found for gene {gene_name}, skipping"
                    )
                    skipped_genes.append(gene_name)
                    continue

                # Get truncation sites
                truncations = self.alt_isoforms.get_visualization_features(gene_name)
                if truncations.empty:
                    print(
                        f"  ├─ No truncation features found for gene {gene_name}, skipping"
                    )
                    skipped_genes.append(gene_name)
                    continue

                print(
                    f"  ├─ Found {len(transcript_info)} transcripts and {len(truncations)} truncation sites"
                )

                # For each transcript, identify overlapping truncations and generate pairs
                gene_pairs = 0
                for _, transcript in transcript_info.iterrows():
                    transcript_id = transcript["transcript_id"]
                    transcript_start = transcript["start"]
                    transcript_end = transcript["end"]
                    transcript_chrom = transcript["chromosome"]

                    # Get CDS regions to determine coding regions
                    cds_features = self.genome.annotations[
                        (self.genome.annotations["transcript_id"] == transcript_id)
                        & (self.genome.annotations["feature_type"] == "CDS")
                    ]

                    if cds_features.empty:
                        print(
                            f"  │  ├─ No CDS features for transcript {transcript_id}, skipping"
                        )
                        continue

                    # Get CDS bounds
                    cds_start = cds_features["start"].min()
                    cds_end = cds_features["end"].max()

                    # Generate canonical protein first
                    canonical_cds = self.build_truncated_cds(
                        transcript_id,
                        pd.Series(
                            {"start": -1, "end": -1}
                        ),  # Special values to get full CDS
                    )

                    if not canonical_cds:
                        print(
                            f"  │  ├─ Failed to get CDS for {transcript_id}, skipping"
                        )
                        continue

                    canonical_protein = self.translate_sequence(canonical_cds)

                    if not canonical_protein:
                        print(
                            f"  │  ├─ Failed to translate canonical sequence for {transcript_id}, skipping"
                        )
                        continue

                    if not (min_length <= len(canonical_protein) <= max_length):
                        print(
                            f"  │  ├─ Canonical protein length {len(canonical_protein)} outside range {min_length}-{max_length}, skipping"
                        )
                        continue

                    # Store canonical sequence
                    canonical_entry = {
                        "gene": gene_name,
                        "transcript_id": transcript_id,
                        "variant_id": "canonical",
                        "sequence": canonical_protein,
                        "length": len(canonical_protein),
                        "is_truncated": 0,
                    }

                    # Check for truncations that overlap with the CDS
                    affecting_truncations = []
                    for idx, trunc in truncations.iterrows():
                        trunc_start = trunc["start"]
                        trunc_end = trunc["end"]
                        trunc_chrom = trunc["chromosome"]

                        # Skip truncations on different chromosomes
                        if trunc_chrom != transcript_chrom:
                            continue

                        # Check for overlap with CDS
                        if not (trunc_end < cds_start or trunc_start > cds_end):
                            trunc_id = f"trunc_{idx}"
                            if "start_codon" in trunc and not pd.isna(
                                trunc["start_codon"]
                            ):
                                trunc_id = f"trunc_{trunc['start_codon']}_{trunc_start}_{trunc_end}"

                            affecting_truncations.append((idx, trunc_id, trunc))

                    if not affecting_truncations:
                        print(
                            f"  │  ├─ No truncations affect CDS of transcript {transcript_id}, skipping"
                        )
                        continue

                    # Process truncations that affect this transcript
                    transcript_pairs = 0
                    for trunc_idx, trunc_id, trunc in affecting_truncations:
                        # Build truncated CDS
                        truncated_cds = self.build_truncated_cds(transcript_id, trunc)

                        if not truncated_cds:
                            print(
                                f"  │  │  ├─ Failed to build truncated CDS for {transcript_id} - {trunc_id}, skipping"
                            )
                            continue

                        if truncated_cds == canonical_cds:
                            print(
                                f"  │  │  ├─ Truncation {trunc_id} does not affect protein sequence of {transcript_id}, skipping"
                            )
                            continue

                        # Translate truncated sequence
                        truncated_protein = self.translate_sequence(truncated_cds)

                        if not truncated_protein:
                            print(
                                f"  │  │  ├─ Failed to translate truncated sequence for {transcript_id} - {trunc_id}, skipping"
                            )
                            continue

                        if not (min_length <= len(truncated_protein) <= max_length):
                            print(
                                f"  │  │  ├─ Truncated protein length {len(truncated_protein)} outside range {min_length}-{max_length}, skipping"
                            )
                            continue

                        # Check if the sequences are actually different
                        if truncated_protein == canonical_protein:
                            print(
                                f"  │  │  ├─ Truncation {trunc_id} results in identical protein, skipping"
                            )
                            continue

                        # Add the canonical sequence for this pair
                        all_sequences.append(canonical_entry)

                        # Add the truncated sequence
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

                        print(
                            f"  │  │  ├─ Generated pair: {transcript_id} + {trunc_id}"
                        )
                        transcript_pairs += 1
                        gene_pairs += 1
                        total_pairs += 1

                    if transcript_pairs > 0:
                        print(
                            f"  │  ├─ Generated {transcript_pairs} pairs for transcript {transcript_id}"
                        )
                    else:
                        # Remove the canonical entry if no truncated versions were created
                        # (This ensures we only have true pairs)
                        all_sequences = [
                            seq
                            for seq in all_sequences
                            if seq["transcript_id"] != transcript_id
                        ]

                if gene_pairs > 0:
                    print(
                        f"  └─ Generated {gene_pairs} total pairs for gene {gene_name}"
                    )
                else:
                    print(f"  └─ No valid pairs generated for gene {gene_name}")
                    skipped_genes.append(gene_name)

            except Exception as e:
                print(f"  └─ Error processing gene {gene_name}: {str(e)}")
                skipped_genes.append(gene_name)

        # Create dataset
        dataset = pd.DataFrame(all_sequences)

        # Save dataset
        if not dataset.empty:
            if "fasta" in output_format.lower():
                self._save_dataset_fasta(dataset)

            if "csv" in output_format.lower():
                output_file = self.output_dir / "protein_sequence_dataset_pairs.csv"
                dataset.to_csv(output_file, index=False)
                print(f"Saved dataset CSV to {output_file}")

        print(f"\nPairs generation summary:")
        print(f"  ├─ Total transcript-truncation pairs: {total_pairs}")
        print(f"  ├─ Total sequences: {len(all_sequences)}")
        print(
            f"  ├─ Genes with pairs: {len(gene_list) - len(skipped_genes)}/{len(gene_list)}"
        )

        if skipped_genes:
            print(f"  └─ Skipped genes: {len(skipped_genes)}")
            if len(skipped_genes) <= 10:
                for gene in skipped_genes:
                    print(f"      - {gene}")
            else:
                for gene in skipped_genes[:10]:
                    print(f"      - {gene}")
                print(f"      - ... and {len(skipped_genes) - 10} more")

        return dataset

    def create_protein_sequence_dataset(
        self,
        gene_list: List[str],
        output_format: str = "fasta,csv",
        include_canonical: bool = True,
        min_length: int = 50,
        max_length: int = 1000,
    ) -> pd.DataFrame:
        """Create a dataset of protein sequences from canonical and truncated transcripts.

        This function processes multiple genes and creates a standardized dataset
        of protein sequences suitable for various analyses.

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

            # Generate sequences
            sequences = self.generate_protein_variants(
                gene_name=gene_name,
                output_format=None,  # Don't save individual files
            )

            # Flatten sequence dictionary for dataset
            for transcript_id, variants in sequences.items():
                # Add canonical sequence if requested
                if include_canonical and "canonical" in variants:
                    canonical_seq = variants["canonical"]

                    if min_length <= len(canonical_seq) <= max_length:
                        all_sequences.append(
                            {
                                "gene": gene_name,
                                "transcript_id": transcript_id,
                                "variant_id": "canonical",
                                "sequence": canonical_seq,
                                "length": len(canonical_seq),
                                "is_truncated": 0,
                            }
                        )

                # Add truncated sequences
                for variant_id, seq in variants.items():
                    if variant_id == "canonical":
                        continue

                    if min_length <= len(seq) <= max_length:
                        all_sequences.append(
                            {
                                "gene": gene_name,
                                "transcript_id": transcript_id,
                                "variant_id": variant_id,
                                "sequence": seq,
                                "length": len(seq),
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

    def _save_sequences(
        self,
        gene_name: str,
        sequences: Dict[str, Dict[str, str]],
        output_format: str = "fasta",
        exclude_canonical: bool = False,
    ) -> None:
        """Save generated sequences to files.

        Args:
            gene_name: Name of the gene
            sequences: Dictionary of sequences keyed by transcript ID
            output_format: Format to save sequences ('fasta' or 'csv')
            exclude_canonical: Whether to exclude canonical transcripts
        """
        gene_dir = self.output_dir / gene_name
        gene_dir.mkdir(exist_ok=True)

        if output_format.lower() == "fasta":
            # Save as FASTA
            records = []

            for transcript_id, variants in sequences.items():
                # Add canonical sequence if requested
                if not exclude_canonical and "canonical" in variants:
                    records.append(
                        SeqRecord(
                            Seq(variants["canonical"]),
                            id=f"{transcript_id}_canonical",
                            description=f"{gene_name} canonical protein",
                        )
                    )

                # Add truncated sequences
                for variant_id, seq in variants.items():
                    if variant_id == "canonical":
                        continue

                    records.append(
                        SeqRecord(
                            Seq(seq),
                            id=f"{transcript_id}_{variant_id}",
                            description=f"{gene_name} truncated protein variant {variant_id}",
                        )
                    )

            # Write to file
            output_file = gene_dir / f"{gene_name}_protein_sequences.fasta"
            SeqIO.write(records, output_file, "fasta")
            print(f"Saved protein sequences to {output_file}")

        elif output_format.lower() == "csv":
            # Save as CSV
            rows = []

            for transcript_id, variants in sequences.items():
                # Add canonical sequence if requested
                if not exclude_canonical and "canonical" in variants:
                    rows.append(
                        {
                            "gene": gene_name,
                            "transcript_id": transcript_id,
                            "variant_id": "canonical",
                            "sequence": variants["canonical"],
                        }
                    )

                # Add truncated sequences
                for variant_id, seq in variants.items():
                    if variant_id == "canonical":
                        continue

                    rows.append(
                        {
                            "gene": gene_name,
                            "transcript_id": transcript_id,
                            "variant_id": variant_id,
                            "sequence": seq,
                        }
                    )

            # Write to file
            output_file = gene_dir / f"{gene_name}_protein_sequences.csv"
            pd.DataFrame(rows).to_csv(output_file, index=False)
            print(f"Saved protein sequences to {output_file}")

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
