"""Transcript translation and protein sequence generation module.

This module contains the TruncatedProteinGenerator class for generating
amino acid sequences from truncated transcripts and preparing datasets
for protein sequence analysis.
"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

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

    def get_preferred_transcript(
        self, gene_name: str, preferred_transcripts: Optional[Set[str]] = None
    ) -> Optional[str]:
        """Get the preferred transcript for a gene.

        Args:
            gene_name: Name of the gene
            preferred_transcripts: Set of preferred transcript IDs

        Returns:
            String transcript ID or None if not found
        """
        transcript_info = self.genome.get_transcript_ids(gene_name)
        if transcript_info.empty:
            return None

        # Filter by preferred transcripts if provided
        if preferred_transcripts:
            filtered = transcript_info[
                transcript_info["transcript_id"].isin(preferred_transcripts)
            ]
            if not filtered.empty:
                transcript_info = filtered

        # Return the first transcript
        return transcript_info.iloc[0]["transcript_id"]

    def extract_canonical_protein(self, transcript_id: str) -> Optional[Dict]:
        """Extract canonical (full) protein sequence using start/stop codon annotations.

        Args:
            transcript_id: Transcript ID to process

        Returns:
            Dict with coding_sequence, protein, and metadata or None if failed
        """
        # Get features and transcript data
        features = self.genome.get_transcript_features(transcript_id)
        transcript_data = self.genome.get_transcript_features_with_sequence(transcript_id)

        if not transcript_data or "sequence" not in transcript_data:
            return None

        transcript_start = transcript_data["sequence"]["start"]
        transcript_end = transcript_data["sequence"]["end"]
        strand = transcript_data["sequence"]["strand"]
        chromosome = transcript_data["sequence"]["chromosome"]

        # Find start_codon and stop_codon annotations
        start_codons = features[features["feature_type"] == "start_codon"]
        stop_codons = features[features["feature_type"] == "stop_codon"]

        if start_codons.empty:
            return None

        start_codon_start = start_codons.iloc[0]["start"]
        start_codon_end = start_codons.iloc[0]["end"]

        stop_codon_start = None
        stop_codon_end = None
        if not stop_codons.empty:
            stop_codon_start = stop_codons.iloc[0]["start"]
            stop_codon_end = stop_codons.iloc[0]["end"]

        # Extract the genomic range
        if strand == "+":
            extract_start = start_codon_start
            extract_end = stop_codon_end if stop_codon_end else None
        else:
            extract_start = stop_codon_start if stop_codon_start else None
            extract_end = start_codon_end

        # Get CDS regions and find overlapping ones
        cds_regions = features[features["feature_type"] == "CDS"].copy()
        overlapping_cds = []

        for _, cds in cds_regions.iterrows():
            cds_start = cds["start"]
            cds_end = cds["end"]

            # Check if this CDS overlaps with our extraction range
            if extract_end and (cds_start > extract_end or cds_end < extract_start):
                continue
            if extract_start and cds_end < extract_start:
                continue
            if extract_end and cds_start > extract_end:
                continue

            # Calculate effective boundaries within this CDS
            effective_start = max(cds_start, extract_start) if extract_start else cds_start
            effective_end = min(cds_end, extract_end) if extract_end else cds_end

            overlapping_cds.append({
                'start': effective_start,
                'end': effective_end,
                'length': effective_end - effective_start + 1
            })

        # Sort CDS regions properly for extraction
        if strand == "+":
            overlapping_cds.sort(key=lambda x: x['start'])
        else:
            overlapping_cds.sort(key=lambda x: x['start'], reverse=True)

        # Extract sequence from each CDS region
        coding_sequence = ""
        for cds in overlapping_cds:
            cds_seq = self.genome.get_sequence(
                chromosome,
                cds['start'],
                cds['end'],
                strand
            )
            coding_sequence += str(cds_seq)

        # Translate
        if len(coding_sequence) >= 3:
            protein = str(Seq(coding_sequence).translate())
        else:
            protein = ""

        return {
            'coding_sequence': coding_sequence,
            'protein': protein,
            'strand': strand,
            'transcript_id': transcript_id
        }

    def extract_truncated_protein(
        self, transcript_id: str, truncation_feature: pd.Series
    ) -> Optional[Dict]:
        """Extract protein sequence starting from the base after truncation.

        Args:
            transcript_id: Transcript ID to process
            truncation_feature: Series with 'start' and 'end' positions of truncation region

        Returns:
            Dict with coding_sequence, protein, and metadata or None if failed
        """
        # Get features and transcript data
        features = self.genome.get_transcript_features(transcript_id)
        transcript_data = self.genome.get_transcript_features_with_sequence(transcript_id)

        if not transcript_data or "sequence" not in transcript_data:
            return None

        transcript_start = transcript_data["sequence"]["start"]
        transcript_end = transcript_data["sequence"]["end"]
        strand = transcript_data["sequence"]["strand"]
        chromosome = transcript_data["sequence"]["chromosome"]

        # Get truncation positions
        truncation_start = truncation_feature["start"]
        truncation_end = truncation_feature["end"]

        # Find the original stop codon (our endpoint)
        stop_codons = features[features["feature_type"] == "stop_codon"]

        stop_codon_start = None
        stop_codon_end = None
        if not stop_codons.empty:
            stop_codon_start = stop_codons.iloc[0]["start"]
            stop_codon_end = stop_codons.iloc[0]["end"]

        # Start at the NEXT base after truncation
        if strand == "+":
            alt_start_pos = truncation_end + 1
            extract_end = stop_codon_end if stop_codon_end else None
        else:
            alt_start_pos = truncation_start
            extract_end = stop_codon_start if stop_codon_start else None

        # Get CDS regions that overlap with our new range
        cds_regions = features[features["feature_type"] == "CDS"].copy()
        overlapping_cds = []

        for _, cds in cds_regions.iterrows():
            cds_start = cds["start"]
            cds_end = cds["end"]

            if strand == "+":
                if cds_end < alt_start_pos:
                    continue
                if extract_end and cds_start > extract_end:
                    continue

                effective_start = max(cds_start, alt_start_pos)
                effective_end = min(cds_end, extract_end) if extract_end else cds_end
            else:
                if cds_start > alt_start_pos:
                    continue
                if extract_end and cds_end < extract_end:
                    continue

                effective_start = max(cds_start, extract_end) if extract_end else cds_start
                effective_end = min(cds_end, alt_start_pos)

            if effective_end >= effective_start:
                overlapping_cds.append({
                    'start': effective_start,
                    'end': effective_end,
                    'length': effective_end - effective_start + 1
                })

        # Sort CDS regions properly for extraction
        if strand == "+":
            overlapping_cds.sort(key=lambda x: x['start'])
        else:
            overlapping_cds.sort(key=lambda x: x['start'], reverse=True)

        # Extract sequence from each CDS region
        coding_sequence = ""
        for cds in overlapping_cds:
            cds_seq = self.genome.get_sequence(
                chromosome,
                cds['start'],
                cds['end'],
                strand
            )
            coding_sequence += str(cds_seq)

        # Translate the sequence directly
        if len(coding_sequence) >= 3:
            # Ensure length is divisible by 3
            remainder = len(coding_sequence) % 3
            if remainder > 0:
                coding_sequence = coding_sequence[:-remainder]

            protein = str(Seq(coding_sequence).translate())
        else:
            protein = ""

        return {
            'coding_sequence': coding_sequence,
            'protein': protein,
            'strand': strand,
            'transcript_id': transcript_id,
            'truncation_start': truncation_start,
            'truncation_end': truncation_end,
            'alternative_start_pos': alt_start_pos,
            'total_cds_regions': len(overlapping_cds)
        }

    def extract_gene_proteins(
        self, 
        gene_name: str, 
        preferred_transcripts: Optional[Set[str]] = None
    ) -> Optional[Dict]:
        """Extract both canonical and truncated proteins for a gene.

        Args:
            gene_name: Name of the gene
            preferred_transcripts: Set of preferred transcript IDs

        Returns:
            Dict with canonical and truncated results or None if failed
        """
        # Get preferred transcript
        transcript_id = self.get_preferred_transcript(gene_name, preferred_transcripts)
        if not transcript_id:
            return None

        # Get truncation features
        truncation_features = self.alt_isoforms.get_visualization_features(gene_name)
        if truncation_features.empty:
            return None

        # Extract canonical protein
        canonical_result = self.extract_canonical_protein(transcript_id)
        if not canonical_result:
            return None

        # Extract truncated protein (using first truncation)
        truncation = truncation_features.iloc[0]
        truncated_result = self.extract_truncated_protein(transcript_id, truncation)
        if not truncated_result:
            return None

        return {
            'gene_name': gene_name,
            'transcript_id': transcript_id,
            'canonical': canonical_result,
            'truncated': truncated_result,
            'truncation': truncation
        }

    def generate_protein_variants(
        self,
        gene_name: str,
        preferred_transcripts: Optional[Set[str]] = None,
        output_format: str = "fasta",
        exclude_canonical: bool = False,
    ) -> Dict[str, Dict[str, str]]:
        """Generate amino acid sequences for all truncation variants of a gene.

        Args:
            gene_name: Name of the gene
            preferred_transcripts: Set of preferred transcript IDs
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
                return result

            # Filter by preferred transcripts if provided
            if preferred_transcripts:
                filtered = transcript_info[
                    transcript_info["transcript_id"].isin(preferred_transcripts)
                ]
                if not filtered.empty:
                    transcript_info = filtered

            # Get truncation sites
            truncations = self.alt_isoforms.get_visualization_features(gene_name)

            # Process each transcript
            for _, transcript in transcript_info.iterrows():
                transcript_id = transcript["transcript_id"]

                # Generate canonical protein
                canonical_result = self.extract_canonical_protein(transcript_id)
                if not canonical_result:
                    continue

                canonical_protein = canonical_result['protein']
                if not canonical_protein:
                    continue

                # Store canonical sequence
                result[transcript_id] = {"canonical": canonical_protein}

                # Process truncations if available
                if not truncations.empty:
                    for _, trunc in truncations.iterrows():
                        trunc_id = f"trunc_{trunc['start']}_{trunc['end']}"
                        if "start_codon" in trunc and not pd.isna(trunc["start_codon"]):
                            trunc_id = f"trunc_{trunc['start_codon']}_{trunc['start']}_{trunc['end']}"

                        # Extract truncated protein
                        truncated_result = self.extract_truncated_protein(transcript_id, trunc)
                        if not truncated_result:
                            continue

                        truncated_protein = truncated_result['protein']
                        if not truncated_protein:
                            continue

                        # Store truncated sequence
                        result[transcript_id][trunc_id] = truncated_protein

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
        preferred_transcripts: Optional[Set[str]] = None,
        output_format: str = "fasta,csv",
        min_length: int = 50,
        max_length: int = 1000,
    ) -> pd.DataFrame:
        """Create a dataset of protein sequences from paired canonical and truncated transcripts.

        Args:
            gene_list: List of gene names to process
            preferred_transcripts: Set of preferred transcript IDs
            output_format: Format to save sequences ('fasta', 'csv', or both)
            min_length: Minimum protein length to include
            max_length: Maximum protein length to include

        Returns:
            DataFrame with the dataset information
        """
        all_sequences = []
        total_pairs = 0
        successful_genes = 0
        skipped_genes = []

        # Process each gene
        for gene_idx, gene_name in enumerate(gene_list, 1):
            try:
                # Extract both canonical and truncated proteins
                gene_result = self.extract_gene_proteins(gene_name, preferred_transcripts)
                
                if not gene_result:
                    skipped_genes.append(gene_name)
                    continue

                canonical_protein = gene_result['canonical']['protein']
                truncated_protein = gene_result['truncated']['protein']
                transcript_id = gene_result['transcript_id']

                # Check length constraints
                if not (min_length <= len(canonical_protein) <= max_length):
                    skipped_genes.append(gene_name)
                    continue

                if not (min_length <= len(truncated_protein) <= max_length):
                    skipped_genes.append(gene_name)
                    continue

                # Check if the sequences are actually different
                if truncated_protein == canonical_protein:
                    skipped_genes.append(gene_name)
                    continue

                # Add canonical sequence
                all_sequences.append({
                    "gene": gene_name,
                    "transcript_id": transcript_id,
                    "variant_id": "canonical",
                    "sequence": canonical_protein,
                    "length": len(canonical_protein),
                    "is_truncated": 0,
                })

                # Add truncated sequence
                trunc_id = f"trunc_{gene_result['truncation']['start']}_{gene_result['truncation']['end']}"
                all_sequences.append({
                    "gene": gene_name,
                    "transcript_id": transcript_id,
                    "variant_id": trunc_id,
                    "sequence": truncated_protein,
                    "length": len(truncated_protein),
                    "is_truncated": 1,
                })

                total_pairs += 1
                successful_genes += 1

            except Exception as e:
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

        print(f"Generated {total_pairs} transcript-truncation pairs from {successful_genes}/{len(gene_list)} genes")
        if skipped_genes:
            print(f"Skipped {len(skipped_genes)} genes due to missing data or constraints")

        return dataset

    def create_protein_sequence_dataset(
        self,
        gene_list: List[str],
        preferred_transcripts: Optional[Set[str]] = None,
        output_format: str = "fasta,csv",
        include_canonical: bool = True,
        min_length: int = 50,
        max_length: int = 1000,
    ) -> pd.DataFrame:
        """Create a dataset of protein sequences from canonical and truncated transcripts.

        Args:
            gene_list: List of gene names to process
            preferred_transcripts: Set of preferred transcript IDs
            output_format: Format to save sequences ('fasta', 'csv', or both)
            include_canonical: Whether to include canonical transcripts as examples
            min_length: Minimum protein length to include
            max_length: Maximum protein length to include

        Returns:
            DataFrame with the dataset information
        """
        all_sequences = []
        successful_genes = 0

        # Process each gene
        for gene_idx, gene_name in enumerate(gene_list, 1):
            # Generate sequences
            sequences = self.generate_protein_variants(
                gene_name=gene_name,
                preferred_transcripts=preferred_transcripts,
                output_format=None,  # Don't save individual files
            )

            if not sequences:
                continue

            successful_genes += 1

            # Flatten sequence dictionary for dataset
            for transcript_id, variants in sequences.items():
                # Add canonical sequence if requested
                if include_canonical and "canonical" in variants:
                    canonical_seq = variants["canonical"]

                    if min_length <= len(canonical_seq) <= max_length:
                        all_sequences.append({
                            "gene": gene_name,
                            "transcript_id": transcript_id,
                            "variant_id": "canonical",
                            "sequence": canonical_seq,
                            "length": len(canonical_seq),
                            "is_truncated": 0,
                        })

                # Add truncated sequences
                for variant_id, seq in variants.items():
                    if variant_id == "canonical":
                        continue

                    if min_length <= len(seq) <= max_length:
                        all_sequences.append({
                            "gene": gene_name,
                            "transcript_id": transcript_id,
                            "variant_id": variant_id,
                            "sequence": seq,
                            "length": len(seq),
                            "is_truncated": 1,
                        })

        # Create dataset
        dataset = pd.DataFrame(all_sequences)

        # Save dataset
        if not dataset.empty:
            if "fasta" in output_format.lower():
                self._save_dataset_fasta(dataset)

            if "csv" in output_format.lower():
                output_file = self.output_dir / "protein_sequence_dataset.csv"
                dataset.to_csv(output_file, index=False)

        print(f"Generated protein dataset from {successful_genes}/{len(gene_list)} genes")
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

        elif output_format.lower() == "csv":
            # Save as CSV
            rows = []

            for transcript_id, variants in sequences.items():
                # Add canonical sequence if requested
                if not exclude_canonical and "canonical" in variants:
                    rows.append({
                        "gene": gene_name,
                        "transcript_id": transcript_id,
                        "variant_id": "canonical",
                        "sequence": variants["canonical"],
                    })

                # Add truncated sequences
                for variant_id, seq in variants.items():
                    if variant_id == "canonical":
                        continue

                    rows.append({
                        "gene": gene_name,
                        "transcript_id": transcript_id,
                        "variant_id": variant_id,
                        "sequence": seq,
                    })

            # Write to file
            output_file = gene_dir / f"{gene_name}_protein_sequences.csv"
            pd.DataFrame(rows).to_csv(output_file, index=False)

    def _save_dataset_fasta(self, dataset: pd.DataFrame) -> None:
        """Save dataset as FASTA file.

        Args:
            dataset: DataFrame with sequence information
        """
        records = []

        for _, row in dataset.iterrows():
            record_id = f"{row['gene']}_{row['transcript_id']}_{row['variant_id']}"
            description = f"{'Truncated' if row['is_truncated'] else 'Canonical'} protein"

            records.append(
                SeqRecord(Seq(row["sequence"]), id=record_id, description=description)
            )

        output_file = self.output_dir / "protein_sequence_dataset.fasta"
        SeqIO.write(records, output_file, "fasta")