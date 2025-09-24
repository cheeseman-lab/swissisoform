"""Enhanced transcript translation and protein sequence generation module.

This module contains the AlternativeProteinGenerator class for generating
amino acid sequences from alternative start sites (both truncations and extensions),
with optional mutation integration capabilities.
"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import asyncio
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from swissisoform.genome import GenomeHandler
from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.mutations import MutationHandler


class AlternativeProteinGenerator:
    """Generate amino acid sequences from alternative start sites.

    This class facilitates the generation of protein sequences stemming from
    alternative start sites (both truncations and extensions) for downstream
    analysis and modeling, with optional mutation integration capabilities.

    Attributes:
        genome (GenomeHandler): Genome sequence and annotation handler.
        alt_isoforms (AlternativeIsoform): Alternative isoform handler.
        output_dir (Path): Directory to save output files.
        mutation_handler (Optional[MutationHandler]): Handler for mutation integration.
        debug (bool): Enable debug mode for detailed output.
    """


    class ValidationCache:
        """Simple in-memory cache for storing variant validation results.

        This class stores consequences of transcript-specific variants along with
        coding sequences and genomic-to-coding position mappings for reuse across
        validation calls.
        """

        def __init__(self):
            """Initialize empty caches for results, coding sequences, and position maps."""
            self.results = {}  # key -> consequence
            self.coding_sequences = {}  # transcript_id -> coding sequence
            self.position_maps = {}  # transcript_id -> {genomic_pos: coding_pos}

        def get_key(self, transcript_id: str, pos: int, ref: str, alt: str) -> str:
            """Generate a unique cache key for a variant.

            Args:
                transcript_id: Identifier of the transcript.
                pos: Genomic position of the variant.
                ref: Reference allele.
                alt: Alternate allele.

            Returns:
                A string key encoding the transcript, position, and allele change.
            """
            return f"{transcript_id}:{pos}:{ref}>{alt}"

        def get_cached_result(self, transcript_id: str, pos: int, ref: str, alt: str):
            """Retrieve a cached consequence for a variant, if available.

            Args:
                transcript_id: Identifier of the transcript.
                pos: Genomic position of the variant.
                ref: Reference allele.
                alt: Alternate allele.

            Returns:
                The cached consequence string, or None if not cached.
            """
            return self.results.get(self.get_key(transcript_id, pos, ref, alt))

        def cache_result(
            self, transcript_id: str, pos: int, ref: str, alt: str, consequence: str
        ):
            """Store the consequence of a variant in the cache.

            Args:
                transcript_id: Identifier of the transcript.
                pos: Genomic position of the variant.
                ref: Reference allele.
                alt: Alternate allele.
                consequence: Predicted consequence to cache.
            """
            self.results[self.get_key(transcript_id, pos, ref, alt)] = consequence

    def __init__(
        self,
        genome_handler: GenomeHandler,
        alt_isoform_handler: AlternativeIsoform,
        output_dir: str,
        mutation_handler: Optional[MutationHandler] = None,
        debug: bool = False,
    ):
        """Initialize the AlternativeProteinGenerator.

        Args:
            genome_handler (GenomeHandler): Initialized GenomeHandler instance.
            alt_isoform_handler (AlternativeIsoform): Initialized AlternativeIsoform instance.
            output_dir (str): Directory to save output files.
            mutation_handler (Optional[MutationHandler]): Optional MutationHandler for mutation integration.
            debug (bool): Enable debug mode for detailed output.

        Returns:
            None
        """
        self.genome = genome_handler
        self.alt_isoforms = alt_isoform_handler
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.mutation_handler = mutation_handler
        self.debug = debug
        self.validation_cache = self.ValidationCache()

    def _debug_print(self, message: str):
        """Print debug message if debug mode is enabled.

        Args:
            message (str): Message to print.

        Returns:
            None
        """
        if self.debug:
            print(f"DEBUG: {message}")

    def extract_canonical_protein(self, transcript_id: str) -> Optional[Dict]:
        """Extract canonical (full) protein sequence using start/stop codon annotations.

        Args:
            transcript_id (str): Transcript ID to process.

        Returns:
            Optional[Dict]: Dictionary with coding_sequence, protein, and metadata, or None if failed.
        """
        # Get features and transcript data
        features = self.genome.get_transcript_features(transcript_id)
        transcript_data = self.genome.get_transcript_features_with_sequence(
            transcript_id
        )

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
            effective_start = (
                max(cds_start, extract_start) if extract_start else cds_start
            )
            effective_end = min(cds_end, extract_end) if extract_end else cds_end

            overlapping_cds.append(
                {
                    "start": effective_start,
                    "end": effective_end,
                    "length": effective_end - effective_start + 1,
                }
            )

        # Sort CDS regions properly for extraction
        if strand == "+":
            overlapping_cds.sort(key=lambda x: x["start"])
        else:
            overlapping_cds.sort(key=lambda x: x["start"], reverse=True)

        # Extract sequence from each CDS region
        coding_sequence = ""
        for cds in overlapping_cds:
            cds_seq = self.genome.get_sequence(
                chromosome, cds["start"], cds["end"], strand
            )
            coding_sequence += str(cds_seq)

        # Translate
        if len(coding_sequence) >= 3:
            protein = str(Seq(coding_sequence).translate())
        else:
            protein = ""

        self._debug_print(
            f"Canonical sequence length: {len(coding_sequence)} bp, {len(protein)} AA"
        )

        return {
            "coding_sequence": coding_sequence,
            "protein": protein,
            "strand": strand,
            "transcript_id": transcript_id,
        }

    def extract_alternative_protein(
        self, transcript_id: str, alternative_feature: pd.Series
    ) -> Optional[Dict]:
        """Extract protein sequence from alternative start site.

        Handles both truncations and extensions based on region_type.

        Args:
            transcript_id (str): Transcript ID to process.
            alternative_feature (pd.Series): Series from get_translation_features() containing region info.

        Returns:
            Optional[Dict]: Dictionary with coding_sequence, protein, and metadata, or None if failed.
        """
        region_type = alternative_feature.get("region_type", "unknown")

        self._debug_print(f"Processing {region_type} for transcript {transcript_id}")

        if region_type == "extension":
            return self._extract_extension_protein(transcript_id, alternative_feature)
        elif region_type == "truncation":
            return self._extract_truncation_protein(transcript_id, alternative_feature)
        else:
            self._debug_print(f"Unknown region type: {region_type}")
            return None

    def _extract_extension_protein(
        self, transcript_id: str, extension_feature: pd.Series
    ) -> Optional[Dict]:
        """Extract full extended protein sequence with improved region detection."""
        # Get features and transcript data
        features = self.genome.get_transcript_features(transcript_id)
        transcript_data = self.genome.get_transcript_features_with_sequence(
            transcript_id
        )

        if not transcript_data or "sequence" not in transcript_data:
            return None

        strand = transcript_data["sequence"]["strand"]
        chromosome = transcript_data["sequence"]["chromosome"]

        # Get extension positions
        extension_start = int(extension_feature["start"])
        extension_end = int(extension_feature["end"])

        # Find canonical start and stop codons
        start_codons = features[features["feature_type"] == "start_codon"]
        stop_codons = features[features["feature_type"] == "stop_codon"]

        if start_codons.empty:
            self._debug_print("No start codon found for transcript")
            return None
        if stop_codons.empty:
            self._debug_print("No stop codon found for transcript")
            return None

        canonical_start_codon_start = int(start_codons.iloc[0]["start"])
        canonical_start_codon_end = int(start_codons.iloc[0]["end"])
        stop_codon_start = int(stop_codons.iloc[0]["start"])
        stop_codon_end = int(stop_codons.iloc[0]["end"])

        self._debug_print(
            f"Extension region: {extension_start}-{extension_end} ({strand} strand)"
        )
        self._debug_print(
            f"Canonical start codon: {canonical_start_codon_start}-{canonical_start_codon_end}"
        )

        # PART 1: Extract extension region - improved logic
        extension_sequence = ""

        # Get all exonic regions (exons and UTRs) within the extension span
        exonic_features = features[
            (
                features["feature_type"].isin(
                    ["exon", "five_prime_utr", "5_prime_utr", "UTR"]
                )
            )
            | (features["feature_type"].str.contains("UTR", na=False))
        ].copy()

        if exonic_features.empty:
            # Fallback to using transcript boundaries
            transcript_start = int(transcript_data["sequence"]["start"])
            transcript_end = int(transcript_data["sequence"]["end"])

            if strand == "+":
                extract_start = max(extension_start, transcript_start)
                extract_end = min(extension_end, canonical_start_codon_start - 1)
            else:
                extract_start = max(extension_start, canonical_start_codon_end + 1)
                extract_end = min(extension_end, transcript_end)

            if extract_start <= extract_end:
                try:
                    seq = self.genome.get_sequence(
                        chromosome, extract_start, extract_end, strand
                    )
                    extension_sequence = str(seq)
                    self._debug_print(
                        f"Fallback extraction: {extract_start}-{extract_end}, {len(extension_sequence)} bp"
                    )
                except Exception as e:
                    self._debug_print(f"Error in fallback extraction: {e}")
        else:
            # Use exonic features for extension extraction
            extension_regions = []

            for _, feature in exonic_features.iterrows():
                feature_start = int(feature["start"])
                feature_end = int(feature["end"])

                # Find overlap with extension region, excluding canonical start codon
                if strand == "+":
                    overlap_start = max(feature_start, extension_start)
                    overlap_end = min(
                        feature_end, min(extension_end, canonical_start_codon_start - 1)
                    )
                else:
                    overlap_start = max(
                        feature_start,
                        max(extension_start, canonical_start_codon_end + 1),
                    )
                    overlap_end = min(feature_end, extension_end)

                if overlap_start <= overlap_end:
                    extension_regions.append((overlap_start, overlap_end))

            # Remove duplicates and sort
            extension_regions = list(set(extension_regions))
            if strand == "+":
                extension_regions.sort(key=lambda x: x[0])
            else:
                extension_regions.sort(key=lambda x: x[0], reverse=True)

            # Extract sequences from each region
            for start, end in extension_regions:
                try:
                    seq = self.genome.get_sequence(chromosome, start, end, strand)
                    seq_str = str(seq)
                    extension_sequence += seq_str

                    self._debug_print(
                        f"Extension region {start}-{end}: {len(seq_str)} bp"
                    )
                except Exception as e:
                    self._debug_print(
                        f"Error extracting extension region {start}-{end}: {e}"
                    )

        self._debug_print(f"Extension sequence length: {len(extension_sequence)} bp")

        # PART 2: Extract CDS regions from canonical start to stop codon
        if strand == "+":
            cds_extract_start = canonical_start_codon_start
            cds_extract_end = stop_codon_end
        else:
            cds_extract_start = stop_codon_start
            cds_extract_end = canonical_start_codon_end

        # Get CDS regions that overlap with canonical range
        cds_regions = features[features["feature_type"] == "CDS"].copy()
        overlapping_cds = []

        for _, cds in cds_regions.iterrows():
            cds_start = int(cds["start"])
            cds_end = int(cds["end"])

            # Check if this CDS overlaps with our extraction range
            if cds_extract_end and (
                cds_start > cds_extract_end or cds_end < cds_extract_start
            ):
                continue
            if cds_extract_start and cds_end < cds_extract_start:
                continue
            if cds_extract_end and cds_start > cds_extract_end:
                continue

            # Calculate effective boundaries within this CDS
            effective_start = (
                max(cds_start, cds_extract_start) if cds_extract_start else cds_start
            )
            effective_end = (
                min(cds_end, cds_extract_end) if cds_extract_end else cds_end
            )

            if effective_end >= effective_start:
                overlapping_cds.append(
                    {
                        "start": effective_start,
                        "end": effective_end,
                        "length": effective_end - effective_start + 1,
                    }
                )

        # Sort CDS regions properly for extraction
        if strand == "+":
            overlapping_cds.sort(key=lambda x: x["start"])
        else:
            overlapping_cds.sort(key=lambda x: x["start"], reverse=True)

        # Extract sequence from each CDS region
        cds_sequence = ""
        for cds in overlapping_cds:
            try:
                cds_seq = self.genome.get_sequence(
                    chromosome, cds["start"], cds["end"], strand
                )
                cds_sequence += str(cds_seq)
            except Exception as e:
                self._debug_print(
                    f"Error extracting CDS {cds['start']}-{cds['end']}: {e}"
                )

        self._debug_print(f"CDS sequence length: {len(cds_sequence)} bp")
        self._debug_print(f"Total CDS regions used: {len(overlapping_cds)}")

        # PART 3: Combine extension + CDS sequences
        full_coding_sequence = extension_sequence + cds_sequence
        self._debug_print(f"Combined sequence length: {len(full_coding_sequence)} bp")

        # Early validation - if we have no extension sequence, this might not be a valid extension
        if len(extension_sequence) == 0:
            self._debug_print(
                "⚠️  WARNING: No extension sequence extracted - this may not be a valid extension"
            )
            # Continue anyway, as this might be an annotation issue

        # Translate the combined sequence
        protein = ""
        if len(full_coding_sequence) >= 3:
            # Ensure length is divisible by 3 for clean translation
            remainder = len(full_coding_sequence) % 3
            if remainder > 0:
                full_coding_sequence = full_coding_sequence[:-remainder]
                self._debug_print(f"Trimmed {remainder} bp for clean translation")

            try:
                protein = str(Seq(full_coding_sequence).translate())
                self._debug_print(f"Extended protein length: {len(protein)} AA")

                # Check for premature stops
                premature_stops = self._count_premature_stops(protein)
                if premature_stops > 0:
                    self._debug_print(
                        f"❌ EXCLUDING: {premature_stops} premature stop codon(s) detected"
                    )
                    return None  # Exit early - don't create this protein
            except Exception as e:
                self._debug_print(f"Translation error: {e}")
                return None

        else:
            self._debug_print("Combined sequence too short for translation")
            return None

        return {
            "coding_sequence": full_coding_sequence,
            "protein": protein,
            "strand": strand,
            "transcript_id": transcript_id,
            "extension_start": extension_start,
            "extension_end": extension_end,
            "alternative_start_pos": extension_start
            if strand == "+"
            else extension_end,
            "extension_sequence_length": len(extension_sequence),
            "cds_sequence_length": len(cds_sequence),
            "total_cds_regions": len(overlapping_cds),
            "extraction_method": "improved_exonic_plus_cds",
            "region_type": "extension",
        }

    def _extract_truncation_protein(
        self, transcript_id: str, truncation_feature: pd.Series
    ) -> Optional[Dict]:
        """Extract protein sequence starting from alternative start position.

        Args:
            transcript_id (str): Transcript ID to process.
            truncation_feature (pd.Series): Series with truncation region information.

        Returns:
            Optional[Dict]: Dictionary with coding_sequence, protein, and metadata, or None if failed.
        """
        # Get features and transcript data
        features = self.genome.get_transcript_features(transcript_id)
        transcript_data = self.genome.get_transcript_features_with_sequence(
            transcript_id
        )

        if not transcript_data or "sequence" not in transcript_data:
            return None

        strand = transcript_data["sequence"]["strand"]
        chromosome = transcript_data["sequence"]["chromosome"]

        # Get alternative start position
        alt_start_pos = truncation_feature.get("alternative_start_pos")
        if alt_start_pos is None:
            self._debug_print("No alternative start position in truncation feature")
            return None

        alt_start_pos = int(alt_start_pos)
        self._debug_print(
            f"Alternative start position: {alt_start_pos} ({strand} strand)"
        )

        # Get CDS regions
        cds_regions = features[features["feature_type"] == "CDS"].copy()
        if cds_regions.empty:
            self._debug_print("No CDS regions found for transcript")
            return None

        # Canonical CDS start for removed region calculation
        if strand == "+":
            canonical_start = cds_regions["start"].min()
            removed_start, removed_end = canonical_start, alt_start_pos - 1
            removed_len = max(0, removed_end - removed_start + 1)
        else:
            canonical_start = cds_regions["end"].max()
            removed_start, removed_end = alt_start_pos + 1, canonical_start
            removed_len = max(0, removed_end - removed_start + 1)

        self._debug_print(
            f"Truncation region (removed): {removed_start}-{removed_end} ({strand} strand)"
        )
        self._debug_print(f"Truncation sequence length: {removed_len} bp")

        # Stop codon for truncation end
        stop_codons = features[features["feature_type"] == "stop_codon"]
        if stop_codons.empty:
            self._debug_print("No stop codon found for transcript")
            return None

        stop_codon_start = stop_codons.iloc[0]["start"]
        stop_codon_end = stop_codons.iloc[0]["end"]

        if strand == "+":
            extract_start = alt_start_pos
            extract_end = stop_codon_end
        else:
            extract_start = stop_codon_start
            extract_end = alt_start_pos

        # Gather CDS overlaps (remaining region to translate)
        overlapping_cds = []
        for _, cds in cds_regions.iterrows():
            cds_start, cds_end = cds["start"], cds["end"]
            if extract_end < cds_start or extract_start > cds_end:
                continue

            effective_start = max(cds_start, extract_start)
            effective_end = min(cds_end, extract_end)

            if effective_end >= effective_start:
                overlapping_cds.append({"start": effective_start, "end": effective_end})

        if not overlapping_cds:
            self._debug_print("No overlapping CDS regions found for remaining sequence")
            return None

        # Order CDS properly
        if strand == "+":
            overlapping_cds.sort(key=lambda x: x["start"])
        else:
            overlapping_cds.sort(key=lambda x: x["start"], reverse=True)

        # ADD THIS DEBUG SECTION (similar to extensions):
        self._debug_print(f"CDS sequence construction for truncation:")
        self._debug_print(f"Total overlapping CDS regions: {len(overlapping_cds)}")

        # Build coding sequence
        coding_sequence = ""
        for i, cds in enumerate(overlapping_cds):
            seq = self.genome.get_sequence(chromosome, cds["start"], cds["end"], strand)
            seq_str = str(seq)
            coding_sequence += seq_str

            # Print each CDS region info (like extensions do)
            self._debug_print(
                f"CDS region {i + 1}: {cds['start']}-{cds['end']}, length {len(seq_str)} bp"
            )

        self._debug_print(
            f"Final truncated CDS sequence length: {len(coding_sequence)} bp"
        )
        self._debug_print(f"Total CDS regions used: {len(overlapping_cds)}")

        # Translate
        protein = ""
        if len(coding_sequence) >= 3:
            remainder = len(coding_sequence) % 3
            if remainder > 0:
                coding_sequence = coding_sequence[:-remainder]

            protein = str(Seq(coding_sequence).translate())
            self._debug_print(f"Truncated protein length: {len(protein)} AA")

            # Check for premature stops and exit
            premature_stops = self._count_premature_stops(protein)
            if premature_stops > 0:
                self._debug_print(
                    f"❌ EXCLUDING: {premature_stops} premature stop codon(s) detected"
                )
                return None  # Exit early - don't create this protein

        else:
            self._debug_print("Coding sequence too short for translation")

        return {
            "coding_sequence": coding_sequence,
            "protein": protein,
            "strand": strand,
            "transcript_id": transcript_id,
            "alternative_start_pos": alt_start_pos,
            "removed_region": (removed_start, removed_end),
            "removed_length": removed_len,
            "remaining_cds_regions": len(overlapping_cds),
            "extraction_method": "alternative_start_plus_cds",
            "region_type": "truncation",
        }

    def extract_gene_proteins(
        self, gene_name: str, top_n_per_type_per_transcript: Optional[int] = 1
    ) -> Optional[List[Dict]]:
        """Extract canonical and alternative proteins for all features of a gene.

        Args:
            gene_name (str): Name of the gene.
            top_n_per_type_per_transcript (Optional[int]): Number of top features per type per transcript.

        Returns:
            Optional[List[Dict]]: List of dicts with canonical and alternative results for each feature, or None if failed.
        """
        self._debug_print(f"Processing gene: {gene_name}")

        # Get alternative features
        alt_features = self.alt_isoforms.get_translation_features(
            gene_name, top_n_per_type_per_transcript
        )
        if alt_features.empty:
            self._debug_print(f"No alternative features found for {gene_name}")
            return None

        all_pairs = []
        canonical_cache = {}  # Cache canonical proteins by transcript_id

        # Process each alternative feature
        for idx, feature in alt_features.iterrows():
            transcript_id = feature["transcript_id"]
            feature_name = feature.get("name", f"feature_{idx}")

            self._debug_print(
                f"Processing feature {feature_name} for transcript {transcript_id}"
            )

            # Get canonical protein (use cache if available)
            if transcript_id not in canonical_cache:
                canonical_result = self.extract_canonical_protein(transcript_id)
                canonical_cache[transcript_id] = canonical_result
            else:
                canonical_result = canonical_cache[transcript_id]

            if not canonical_result:
                self._debug_print(
                    f"Could not extract canonical protein for {transcript_id}"
                )
                continue

            # Extract alternative protein for this specific feature
            alternative_result = self.extract_alternative_protein(
                transcript_id, feature
            )
            if not alternative_result:
                self._debug_print(
                    f"Could not extract alternative protein for feature {feature_name}"
                )
                continue

            # ===== ADD PREMATURE STOP CHECK FOR ALTERNATIVE PROTEINS =====
            alternative_protein = alternative_result["protein"]
            premature_stops = self._count_premature_stops(alternative_protein)
            region_type = feature.get("region_type", "unknown")

            if premature_stops > 0:
                region_type = feature.get("region_type", "unknown")
                self._debug_print(
                    f"❌ BLOCKING {region_type.upper()}: {premature_stops} premature stop codon(s) detected"
                )
                self._debug_print(
                    f"❌ Gene: {gene_name}, Transcript: {transcript_id}, Feature: {feature_name}"
                )
                continue  # Skip this feature entirely

            self._debug_print(
                f"✅ {region_type.upper()} protein clean: {len(alternative_protein)} AA"
            )
            # ===== END PREMATURE STOP CHECK =====

            all_pairs.append(
                {
                    "gene_name": gene_name,
                    "transcript_id": transcript_id,
                    "feature_id": feature_name,
                    "region_type": feature.get("region_type", "unknown"),
                    "canonical": canonical_result,
                    "alternative": alternative_result,
                    "feature_info": feature.to_dict(),
                }
            )

        return all_pairs if all_pairs else None

    # ===== MUTATION INTEGRATION METHODS =====

    def _is_single_bp_substitution(self, hgvsc: str) -> bool:
        """Check if HGVS represents a single BP substitution we want to process.

        Args:
            hgvsc (str): HGVS coding notation (e.g., c.76C>T, c.-25G>A).

        Returns:
            bool: True if this is a single BP substitution pattern.
        """
        if not hgvsc or pd.isna(hgvsc) or str(hgvsc) == "nan":
            return False

        hgvsc = str(hgvsc).strip()

        # Match single BP substitution patterns: c.76C>T, c.-25G>A, etc.
        # Accepts: c.123A>T, c.-45G>C
        # Rejects: c.76delA, c.76_77insG, c.76_78delATC, c.76A>TG
        pattern = r"^c\.-?\d+[ATCG]>[ATCG]$"

        return bool(re.match(pattern, hgvsc, re.IGNORECASE))

    def _validate_single_bp_variant(self, mutation_row: pd.Series) -> bool:
        """Validate if this mutation represents a single BP change we want to process.

        Args:
            mutation_row (pd.Series): Row from mutations DataFrame.

        Returns:
            bool: True if this variant should be processed.
        """
        impact = mutation_row.get("impact", "")
        hgvsc = mutation_row.get("hgvsc", "")

        # Only process specific impact types
        valid_impacts = ["missense variant", "5 prime UTR variant"]
        if impact not in valid_impacts:
            return False

        # Validate single BP substitution pattern
        if not self._is_single_bp_substitution(hgvsc):
            return False

        # Additional validation: ensure we have required columns
        required_cols = ["position", "reference", "alternate"]
        for col in required_cols:
            if col not in mutation_row or pd.isna(mutation_row[col]):
                return False

        return True

    def _deduplicate_mutations(self, mutations: pd.DataFrame) -> pd.DataFrame:
        """Remove duplicate variants based on genomic position and alleles.

        Args:
            mutations (pd.DataFrame): DataFrame with mutations.

        Returns:
            pd.DataFrame: Deduplicated DataFrame.
        """
        if mutations.empty:
            return mutations

        # Define uniqueness criteria - genomic coordinates
        unique_cols = ["position", "reference", "alternate"]

        # Keep first occurrence
        deduplicated = mutations.drop_duplicates(subset=unique_cols, keep="first")

        if len(deduplicated) < len(mutations):
            dropped_count = len(mutations) - len(deduplicated)
            self._debug_print(f"Removed {dropped_count} duplicate variants")

        return deduplicated

    def _apply_mutation_to_sequence(
        self,
        transcript_id: str,
        genomic_pos: int,
        ref_allele: str,
        alt_allele: str,
        sequence_type: str = "canonical",
        extension_feature: Optional[pd.Series] = None,
    ) -> Optional[Dict]:
        """Apply mutation to canonical or extension sequence with improved multi-base handling."""
        if self.debug:
            print(
                f"            ├─ Applying mutation {ref_allele}>{alt_allele} at position {genomic_pos} to {sequence_type} sequence"
            )

        # Get the base sequence and metadata
        if sequence_type == "extension":
            if extension_feature is None:
                if self.debug:
                    print(
                        "            └─ ❌ Extension feature required for extension mutations"
                    )
                return None
            base_result = self.extract_alternative_protein(
                transcript_id, extension_feature
            )
            if not base_result:
                if self.debug:
                    print("            └─ ❌ Could not extract base extension protein")
                return None
        else:
            base_result = self.extract_canonical_protein(transcript_id)
            if not base_result:
                if self.debug:
                    print("            └─ ❌ Could not extract base canonical protein")
                return None

        original_coding_sequence = base_result["coding_sequence"]
        original_protein = base_result["protein"]
        strand = base_result["strand"]

        if self.debug:
            print(
                f"            ├─ Base sequence: {len(original_coding_sequence)} bp CDS, {len(original_protein)} AA"
            )
            print(f"            ├─ Transcript strand: {strand}")

        # Get transcript data for validation and mapping
        transcript_data = self.genome.get_transcript_features_with_sequence(
            transcript_id
        )
        if not transcript_data:
            if self.debug:
                print("            └─ ❌ No transcript data available")
            return None

        chromosome = transcript_data["sequence"]["chromosome"]

        # Validate genomic reference allele
        if ref_allele and len(ref_allele) > 0:
            try:
                actual_genomic_forward = self.genome.get_sequence(
                    chromosome, genomic_pos, genomic_pos + len(ref_allele) - 1, "+"
                )
                actual_genomic_forward = str(actual_genomic_forward).upper()

                if self.debug:
                    print(
                        f"            ├─ Genomic sequence (+ strand): {actual_genomic_forward}"
                    )
                    print(f"            ├─ Expected reference: {ref_allele}")

                if actual_genomic_forward != ref_allele.upper():
                    if self.debug:
                        print(
                            f"            └─ ❌ Reference allele mismatch: expected {ref_allele}, found {actual_genomic_forward}"
                        )
                    return None

            except Exception as e:
                if self.debug:
                    print(
                        f"            └─ ❌ Could not get sequence at position {genomic_pos}: {e}"
                    )
                return None

        # Convert alleles to transcript orientation
        complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}

        if strand == "+":
            transcript_ref = ref_allele.upper()
            transcript_alt = alt_allele.upper()
        else:
            # FIX: Proper reverse complement for multi-base sequences
            transcript_ref = (
                "".join([complement_map.get(b, b) for b in ref_allele.upper()[::-1]])
                if ref_allele
                else ""
            )
            transcript_alt = (
                "".join([complement_map.get(b, b) for b in alt_allele.upper()[::-1]])
                if alt_allele
                else ""
            )

        if self.debug:
            print(
                f"            ├─ Transcript alleles: {transcript_ref}>{transcript_alt}"
            )

        # Map genomic position to coding position
        mutation_coding_pos = self._map_genomic_to_coding_position(
            transcript_id,
            genomic_pos,
            original_coding_sequence,
            sequence_type,
            extension_feature,
        )

        if mutation_coding_pos is None:
            if self.debug:
                print(
                    f"            └─ ❌ Could not map genomic position {genomic_pos} to coding sequence"
                )
            return None

        if self.debug:
            print(f"            ├─ Mapped to coding position: {mutation_coding_pos}")

        # FIX: For negative strand multi-base mutations, we need to be more careful about validation
        ref_len = len(transcript_ref)
        if mutation_coding_pos + ref_len > len(original_coding_sequence):
            if self.debug:
                print(
                    f"            └─ ❌ Coding position {mutation_coding_pos}+{ref_len} beyond sequence length {len(original_coding_sequence)}"
                )
            return None

        # Validate reference in coding sequence
        if transcript_ref:
            current_bases = original_coding_sequence[
                mutation_coding_pos : mutation_coding_pos + ref_len
            ].upper()

            if current_bases != transcript_ref:
                if self.debug:
                    print(
                        f"            ├─ ❌ Transcript reference mismatch at coding pos {mutation_coding_pos}"
                    )
                    print(
                        f"            ├─     Expected: {transcript_ref} (from genomic {ref_allele})"
                    )
                    print(f"            ├─     Found: {current_bases}")
                    print(f"            ├─     Strand: {strand}")

                    # FIX: For multi-base negative strand, try adjusting position by 1
                    if (
                        strand == "-"
                        and len(ref_allele) > 1
                        and mutation_coding_pos > 0
                    ):
                        # Try position - 1
                        alt_pos = mutation_coding_pos - 1
                        if alt_pos + ref_len <= len(original_coding_sequence):
                            alt_current_bases = original_coding_sequence[
                                alt_pos : alt_pos + ref_len
                            ].upper()
                            if alt_current_bases == transcript_ref:
                                if self.debug:
                                    print(
                                        f"            ├─ ✅ Found match at adjusted position {alt_pos}"
                                    )
                                mutation_coding_pos = alt_pos
                                current_bases = alt_current_bases
                            else:
                                if self.debug:
                                    print(
                                        f"            ├─     Tried pos-1 ({alt_pos}): {alt_current_bases}"
                                    )

                    # If still no match, try other positions within a small window
                    if current_bases != transcript_ref and len(ref_allele) > 1:
                        for offset in [-2, -1, 1, 2]:
                            test_pos = mutation_coding_pos + offset
                            if 0 <= test_pos and test_pos + ref_len <= len(
                                original_coding_sequence
                            ):
                                test_bases = original_coding_sequence[
                                    test_pos : test_pos + ref_len
                                ].upper()
                                if test_bases == transcript_ref:
                                    if self.debug:
                                        print(
                                            f"            ├─ ✅ Found match at offset {offset} (position {test_pos})"
                                        )
                                    mutation_coding_pos = test_pos
                                    current_bases = test_bases
                                    break
                                elif self.debug:
                                    print(
                                        f"            ├─     Tried offset {offset} (pos {test_pos}): {test_bases}"
                                    )

                    # If we still don't have a match, fail
                    if current_bases != transcript_ref:
                        if self.debug:
                            print(
                                f"            └─ ❌ Could not find matching sequence after position adjustments"
                            )
                        return None

            if self.debug:
                print(
                    f"            ├─ Coding sequence reference validated: {current_bases} == {transcript_ref}"
                )

        # Apply mutation in transcript orientation
        mutated_coding_sequence = (
            original_coding_sequence[:mutation_coding_pos]
            + transcript_alt
            + original_coding_sequence[mutation_coding_pos + ref_len :]
        )

        if self.debug:
            print(
                f"            ├─ Applied mutation at coding pos {mutation_coding_pos}: {current_bases}>{transcript_alt}"
            )

        # Translate mutated sequence
        if len(mutated_coding_sequence) >= 3:
            remainder = len(mutated_coding_sequence) % 3
            if remainder > 0:
                mutated_coding_sequence = mutated_coding_sequence[:-remainder]
            mutated_protein = str(Seq(mutated_coding_sequence).translate())
        else:
            mutated_protein = ""

        # Check if protein actually changed
        if mutated_protein == original_protein:
            if self.debug:
                print(f"            └─ Synonymous mutation - protein unchanged")
            return None

        # Calculate amino acid difference
        aa_change = self._calculate_aa_difference(original_protein, mutated_protein)

        if self.debug:
            print(f"            ├─ Protein changed! AA change: {aa_change}")
            print(f"            └─ ✅ Mutation application successful")

        return {
            "coding_sequence": mutated_coding_sequence,
            "protein": mutated_protein,
            "strand": strand,
            "transcript_id": transcript_id,
            "sequence_type": sequence_type,
            "mutation_position": genomic_pos,
            "mutation_ref": ref_allele,
            "mutation_alt": alt_allele,
            "coding_position": mutation_coding_pos,
            "aa_change": aa_change,
            "protein_changed": True,
        }

    def _map_genomic_to_coding_position(
        self,
        transcript_id: str,
        genomic_pos: int,
        coding_sequence: str,
        sequence_type: str,
        extension_feature: Optional[pd.Series] = None,
    ) -> Optional[int]:
        """Improved genomic to coding position mapping."""
        # Use context-aware cache key
        if extension_feature is not None and sequence_type == "extension":
            cache_key = f"{transcript_id}_extension"
        else:
            cache_key = f"{transcript_id}_canonical"

        # Build cache if needed
        if cache_key not in self.validation_cache.position_maps:
            if sequence_type == "extension" and extension_feature is not None:
                pos_map = self._build_extension_position_map(
                    transcript_id, extension_feature
                )
            else:
                pos_map = self._build_canonical_position_map(transcript_id)

            self.validation_cache.position_maps[cache_key] = pos_map
        else:
            pos_map = self.validation_cache.position_maps[cache_key]

        # Get coding position
        coding_pos = pos_map.get(genomic_pos)

        if self.debug and coding_pos is not None:
            # Add context information
            region_type = "unknown"
            if sequence_type == "extension" and extension_feature is not None:
                extension_start = int(extension_feature["start"])
                extension_end = int(extension_feature["end"])

                features = self.genome.get_transcript_features(transcript_id)
                start_codons = features[features["feature_type"] == "start_codon"]
                if not start_codons.empty:
                    # Fix: Get transcript_data properly
                    transcript_data = self.genome.get_transcript_features_with_sequence(
                        transcript_id
                    )
                    canonical_start = (
                        start_codons.iloc[0]["start"]
                        if transcript_data
                        and transcript_data.get("sequence", {}).get("strand") == "+"
                        else start_codons.iloc[0]["end"]
                    )

                    if (
                        extension_start <= genomic_pos <= extension_end
                        and genomic_pos != canonical_start
                    ):
                        region_type = "extension"
                    else:
                        region_type = "CDS"
            else:
                region_type = "CDS"

            print(
                f"            ├─ Position {genomic_pos} mapped to coding pos {coding_pos} ({region_type})"
            )

            # Show sequence context
            if coding_pos < len(coding_sequence):
                context_start = max(0, coding_pos - 5)
                context_end = min(len(coding_sequence), coding_pos + 6)
                context = coding_sequence[context_start:context_end]
                marker_pos = coding_pos - context_start
                context_with_marker = (
                    context[:marker_pos]
                    + f"[{context[marker_pos]}]"
                    + context[marker_pos + 1 :]
                )
                print(f"            ├─ Coding sequence context: {context_with_marker}")

        return coding_pos

    def _map_genomic_to_cds_position(
        self,
        transcript_id: str,
        genomic_pos: int,
        coding_sequence: str,
        sequence_type: str,
        extension_feature: Optional[pd.Series] = None,
    ) -> Optional[int]:
        """Map genomic position to CDS position, accounting for any extension offset."""
        features = self.genome.get_transcript_features(transcript_id)
        transcript_data = self.genome.get_transcript_features_with_sequence(
            transcript_id
        )
        strand = transcript_data["sequence"]["strand"]

        start_codons = features[features["feature_type"] == "start_codon"]
        canonical_start = (
            start_codons.iloc[0]["start"]
            if strand == "+"
            else start_codons.iloc[0]["end"]
        )

        cds_regions = features[features["feature_type"] == "CDS"].copy()
        if cds_regions.empty:
            return None

        # Sort CDS regions the same way
        if strand == "+":
            cds_regions = cds_regions.sort_values("start")
        else:
            cds_regions = cds_regions.sort_values("start", ascending=False)

        # Calculate extension offset using UTR regions (not genomic span)
        extension_offset = 0
        if sequence_type == "extension" and extension_feature is not None:
            utr_features = features[features["feature_type"].str.contains("UTR")]

            for _, utr in utr_features.iterrows():
                utr_start = utr["start"]
                utr_end = utr["end"]

                if strand == "+":
                    overlap_start = max(utr_start, extension_feature["start"])
                    overlap_end = min(utr_end, canonical_start - 1)
                else:
                    overlap_start = max(utr_start, canonical_start + 1)
                    overlap_end = min(utr_end, extension_feature["end"])

                if overlap_start <= overlap_end:
                    extension_offset += overlap_end - overlap_start + 1

        # Find position within CDS regions using EXACT same logic as position map building
        coding_pos = 0

        for _, cds in cds_regions.iterrows():
            cds_start = int(cds["start"])
            cds_end = int(cds["end"])

            # Check if mutation falls within this CDS
            if cds_start <= genomic_pos <= cds_end:
                if strand == "+":
                    effective_start = max(cds_start, canonical_start)
                    if effective_start <= genomic_pos:
                        offset_in_cds = genomic_pos - effective_start
                        final_pos = extension_offset + coding_pos + offset_in_cds

                        if self.debug:
                            print(
                                f"            ├─ Mapped to coding position: {final_pos}"
                            )
                            if final_pos < len(coding_sequence):
                                context = coding_sequence[
                                    max(0, final_pos - 5) : final_pos + 6
                                ]
                                print(f"            ├─ Sequence context: {context}")

                        return final_pos
                else:
                    effective_end = min(cds_end, canonical_start)
                    if genomic_pos <= effective_end:
                        # FIX: Match the exact logic from _build_canonical_position_map
                        # The position map assigns: effective_end → coding_pos, effective_end-1 → coding_pos+1, etc.
                        offset_in_cds = effective_end - genomic_pos - 1
                        final_pos = extension_offset + coding_pos + offset_in_cds

                        if self.debug:
                            print(
                                f"            ├─ Mapped to coding position: {final_pos}"
                            )
                            if final_pos < len(coding_sequence):
                                context = coding_sequence[
                                    max(0, final_pos - 5) : final_pos + 6
                                ]
                                print(f"            ├─ Sequence context: {context}")

                        return final_pos

            # Add this CDS length using EXACT same logic as position map building
            if strand == "+":
                effective_start = max(cds_start, canonical_start)
                effective_end = cds_end
                if effective_start <= effective_end:
                    coding_pos += effective_end - effective_start + 1
            else:
                effective_start = cds_start
                effective_end = min(cds_end, canonical_start)
                if effective_start <= effective_end:
                    coding_pos += effective_end - effective_start + 1

        return None

    async def _get_mutations_in_region(
        self,
        gene_name: str,
        start: int,
        end: int,
        sources: List[str] = None,
        impact_types: List[str] = None,
    ) -> pd.DataFrame:
        """Get filtered, deduplicated single BP mutations within a specific genomic region.

        Args:
            gene_name (str): Gene name.
            start (int): Region start position.
            end (int): Region end position.
            sources (List[str], optional): List of mutation sources to include.
            impact_types (List[str], optional): List of impact types to include.

        Returns:
            pd.DataFrame: Filtered and deduplicated DataFrame of mutations.
        """
        self._debug_print(f"Fetching mutations for {gene_name} in region {start}-{end}")

        if not self.mutation_handler:
            self._debug_print("No mutation handler available")
            return pd.DataFrame()

        # Set defaults
        if sources is None:
            sources = ["clinvar"]
        if impact_types is None:
            impact_types = ["missense variant"]

        # Create region features for mutation handler
        region_features = pd.DataFrame(
            [
                {
                    "start": start,
                    "end": end,
                    "chromosome": "chr1",  # Will be corrected by mutation handler
                    "name": f"region_{start}_{end}",
                }
            ]
        )

        # Get mutations in this region from specified sources
        mutations = await self.mutation_handler.get_visualization_ready_mutations(
            gene_name=gene_name, alt_features=region_features, sources=sources
        )

        if mutations is None or mutations.empty:
            self._debug_print(f"No mutations found in region")
            return pd.DataFrame()

        # Filter to user-specified impact types
        mutations = mutations[mutations["impact"].isin(impact_types)].copy()

        if mutations.empty:
            self._debug_print(
                f"No mutations with specified impact types: {impact_types}"
            )
            return pd.DataFrame()

        # Filter to only single BP changes using validation
        # valid_mask = mutations.apply(self._validate_single_bp_variant, axis=1)
        # mutations = mutations[valid_mask].copy()

        if mutations.empty:
            self._debug_print(f"No valid single BP variants found")
            return pd.DataFrame()

        # Deduplicate by genomic coordinates
        mutations = self._deduplicate_mutations(mutations)

        self._debug_print(f"Found {len(mutations)} variants")
        return mutations

    def _calculate_aa_difference(
        self, original_protein: str, mutated_protein: str
    ) -> Optional[str]:
        """Calculate amino acid difference between two protein sequences.

        Args:
            original_protein (str): Original protein sequence.
            mutated_protein (str): Mutated protein sequence.

        Returns:
            Optional[str]: String in format "A123V" (original AA, position, new AA) or None if no single change.
        """
        if len(original_protein) != len(mutated_protein):
            # Handle length differences (indels, frameshifts, etc.)
            return f"len_{len(original_protein)}>{len(mutated_protein)}"

        differences = []
        for i, (orig_aa, mut_aa) in enumerate(zip(original_protein, mutated_protein)):
            if orig_aa != mut_aa:
                # Use 1-based positioning like standard notation
                differences.append(f"{orig_aa}{i + 1}{mut_aa}")

        if len(differences) == 0:
            return None  # No difference
        elif len(differences) == 1:
            return differences[0]  # Single AA change
        else:
            # Multiple changes - could truncate or show first few
            if len(differences) <= 3:
                return "_".join(differences)
            else:
                return f"{differences[0]}_and_{len(differences) - 1}_more"

    async def extract_gene_proteins_with_mutations(
        self,
        gene_name: str,
        include_mutations: bool = True,
        sources: List[str] = None,
        impact_types: List[str] = None,
        pre_validated_variants: Optional[Dict[str, Set[str]]] = None,  # NEW
        skip_validation: bool = False,  # NEW
    ) -> Optional[List[Dict]]:
        """Extract proteins for all alternative features with optional mutation integration.

        Args:
            gene_name (str): Name of the gene.
            include_mutations (bool): Whether to include mutation variants.
            sources (List[str], optional): List of mutation sources to include.
            impact_types (List[str], optional): Mutation impact types to include.
            pre_validated_variants (Optional[Dict[str, Set[str]]]): Pre-validated variant IDs by gene.
            skip_validation (bool): Skip validation loop (use pre-validated results).

        Returns:
            Optional[List[Dict]]: List of enhanced protein pair dictionaries, or None if failed.
        """
        self._debug_print(
            f"Processing gene {gene_name} with mutations={include_mutations}"
        )

        # Set defaults
        if sources is None:
            sources = ["clinvar"]
        if impact_types is None:
            impact_types = ["missense variant"]

        # Get base pairs using the simplified logic
        base_pairs = self.extract_gene_proteins(gene_name)
        if not base_pairs:
            self._debug_print(f"No base pairs found for {gene_name}")
            return None

        enhanced_pairs = []
        for pair in base_pairs:
            enhanced_pair = {
                "gene_name": gene_name,
                "transcript_id": pair["transcript_id"],
                "feature_id": pair["feature_id"],
                "region_type": pair["region_type"],
                "canonical": pair["canonical"],
                "alternative_base": pair["alternative"],
                "alternative_mutations": [],
            }

            if include_mutations and self.mutation_handler:
                # Get mutations in this alternative region
                feature_info = pair["feature_info"]
                mutations = await self._get_mutations_in_region(
                    gene_name=gene_name,
                    start=feature_info["start"],
                    end=feature_info["end"],
                    sources=sources,
                    impact_types=impact_types
                    if pre_validated_variants is None
                    else None,
                )

                # NEW: Filter to only pre-validated variants if provided
                if (
                    not mutations.empty
                    and pre_validated_variants
                    and gene_name in pre_validated_variants
                ):
                    pre_validated_ids = pre_validated_variants[gene_name]
                    mutations = mutations[
                        mutations["variant_id"].isin(pre_validated_ids)
                    ].copy()
                    self._debug_print(
                        f"Filtered to {len(mutations)} pre-validated mutations"
                    )

                if not mutations.empty:
                    self._debug_print(
                        f"Found {len(mutations)} mutations in alternative region"
                    )

                    # Process each mutation
                    for _, mutation in mutations.iterrows():
                        try:
                            variant_id = mutation.get("variant_id", "unknown")
                            print(f"            ├─ Processing variant: {variant_id}")
                            # Extract mutation data
                            genomic_pos = int(mutation["position"])
                            ref_allele = str(mutation["reference"]).upper()
                            alt_allele = str(mutation["alternate"]).upper()
                            impact = mutation.get("impact", "")
                            hgvsc = str(mutation.get("hgvsc", ""))
                            hgvsp = str(mutation.get("hgvsp", ""))

                            # Determine sequence type based on impact and region
                            sequence_type = "canonical"
                            extension_feature = None

                            if pair["region_type"] == "extension":
                                sequence_type = "extension"
                                extension_feature = pd.Series(feature_info)

                            # Apply mutation using unified method
                            mutation_result_data = self._apply_mutation_to_sequence(
                                pair["transcript_id"],
                                genomic_pos,
                                ref_allele,
                                alt_allele,
                                sequence_type,
                                extension_feature,
                            )

                            if mutation_result_data:
                                # The unified method already filtered out synonymous mutations
                                # and calculated aa_change

                                variant_type = f"{sequence_type}_mutated"

                                mutation_result = {
                                    "coding_sequence": mutation_result_data[
                                        "coding_sequence"
                                    ],
                                    "protein": mutation_result_data["protein"],
                                    "transcript_id": pair["transcript_id"],
                                    "variant_type": variant_type,
                                    "mutation": {
                                        "position": genomic_pos,
                                        "reference": ref_allele,
                                        "alternate": alt_allele,
                                        "hgvsc": hgvsc,
                                        "hgvsp": hgvsp,
                                        "impact": impact,
                                        # NEW: Use original impact if skipping validation
                                        "impact_validated": impact
                                        if skip_validation
                                        else mutation_result_data.get(
                                            "impact_validated", impact
                                        ),
                                        "variant_id": mutation.get("variant_id", ""),
                                        "source": mutation.get("source", ""),
                                        "calculated_aa_change": mutation_result_data.get(
                                            "aa_change", ""
                                        ),
                                        # NEW: Add validation status
                                        "validation_skipped": skip_validation,
                                    },
                                }

                                enhanced_pair["alternative_mutations"].append(
                                    mutation_result
                                )
                                mode = (
                                    "pre-validated" if skip_validation else "validated"
                                )
                                self._debug_print(
                                    f"Successfully created {sequence_type} mutation variant ({mode}) with AA change: {mutation_result_data.get('aa_change', 'N/A')}"
                                )
                            else:
                                self._debug_print(
                                    f"Mutation at {genomic_pos} did not result in protein change (synonymous or failed)"
                                )

                        except Exception as e:
                            self._debug_print(f"Error creating mutation variant: {e}")
                            continue
                else:
                    self._debug_print("No valid mutations found in this region")

            enhanced_pairs.append(enhanced_pair)

        return enhanced_pairs

    # ===== DATASET GENERATION METHODS =====

    def generate_protein_variants(
        self,
        gene_name: str,
        output_format: str = "fasta",
        exclude_canonical: bool = False,
    ) -> Dict[str, Dict[str, str]]:
        """Generate amino acid sequences for all alternative variants of a gene.

        Args:
            gene_name (str): Name of the gene.
            output_format (str): Format to save sequences ('fasta', 'csv', or None for no save).
            exclude_canonical (bool): Whether to exclude canonical transcripts.

        Returns:
            Dict[str, Dict[str, str]]: Dictionary with transcript IDs as keys, containing both canonical and alternative protein sequences for each transcript.
        """
        result = {}

        try:
            # Get alternative features - transcript IDs are embedded
            alt_features = self.alt_isoforms.get_translation_features(gene_name)
            if alt_features.empty:
                return result

            canonical_cache = {}

            # Process each feature
            for _, feature in alt_features.iterrows():
                transcript_id = feature["transcript_id"]
                feature_name = feature.get("name", f"feature_{feature.name}")
                region_type = feature.get("region_type", "unknown")

                # Get canonical protein (use cache)
                if transcript_id not in canonical_cache:
                    canonical_result = self.extract_canonical_protein(transcript_id)
                    canonical_cache[transcript_id] = canonical_result
                else:
                    canonical_result = canonical_cache[transcript_id]

                if not canonical_result:
                    continue

                canonical_protein = canonical_result["protein"]
                if not canonical_protein:
                    continue

                # Initialize transcript entry if not exists
                if transcript_id not in result:
                    result[transcript_id] = {"canonical": canonical_protein}

                # Extract alternative protein
                alternative_result = self.extract_alternative_protein(
                    transcript_id, feature
                )
                if not alternative_result:
                    continue

                alternative_protein = alternative_result["protein"]
                if not alternative_protein:
                    continue

                # Store alternative sequence with descriptive ID
                alt_id = f"{region_type}_{feature_name}"
                result[transcript_id][alt_id] = alternative_protein

            # Save results if requested
            if output_format and result:
                self._save_sequences(
                    gene_name, result, output_format, exclude_canonical
                )

            return result

        except Exception as e:
            print(f"Error generating sequences for gene {gene_name}: {str(e)}")
            return result

    async def create_protein_sequence_dataset_with_mutations(
        self,
        gene_list: List[str],
        include_mutations: bool = True,
        sources: List[str] = None,
        impact_types: List[str] = None,
        output_format: str = "fasta,csv",
        min_length: int = 10,
        max_length: int = 100000,
        pre_validated_variants: Optional[Dict[str, Set[str]]] = None,
        skip_validation: bool = False,
    ) -> pd.DataFrame:
        """Create a dataset with correct comparison sets for each region type."""
        all_sequences = []
        successful_genes = 0
        skipped_genes = []

        # Set defaults
        if sources is None:
            sources = ["clinvar"]
        if impact_types is None:
            impact_types = ["missense variant"]

        # Process each gene
        for gene_idx, gene_name in enumerate(gene_list, 1):
            try:
                self._debug_print(
                    f"Processing gene {gene_idx}/{len(gene_list)}: {gene_name}"
                )

                if include_mutations and self.mutation_handler:
                    enhanced_pairs = await self.extract_gene_proteins_with_mutations(
                        gene_name,
                        include_mutations=include_mutations,
                        sources=sources,
                        impact_types=impact_types,
                        pre_validated_variants=pre_validated_variants,
                        skip_validation=skip_validation,
                    )
                    if not enhanced_pairs:
                        skipped_genes.append(gene_name)
                        continue

                    for pair in enhanced_pairs:
                        canonical_protein = pair["canonical"]["protein"]
                        alternative_protein = pair["alternative_base"]["protein"]
                        transcript_id = pair["transcript_id"]
                        feature_id = pair["feature_id"]
                        region_type = pair["region_type"]

                        # Check length constraints
                        if not (min_length <= len(canonical_protein) <= max_length):
                            continue
                        if not (min_length <= len(alternative_protein) <= max_length):
                            continue
                        if alternative_protein == canonical_protein:
                            continue

                        # ===== COMPARISON SET 1: CANONICAL (always included) =====
                        all_sequences.append(
                            {
                                "gene": gene_name,
                                "transcript_id": transcript_id,
                                "variant_id": "canonical",
                                "sequence": canonical_protein,
                                "length": len(canonical_protein),
                                "variant_type": "canonical",
                                "region_type": region_type,
                                "comparison_set": f"{region_type}_analysis",  # Add comparison set identifier
                                # Mutation fields
                                "mutation_position": None,
                                "mutation_change": None,
                                "aa_change": None,
                                "hgvsc": None,
                                "hgvsp": None,
                                "mutation_impact": None,
                                "mutation_impact_validated": None,
                                "in_alt_start_site": None,
                                "mutation_source": None,
                                "clinvar_variant_id": None,
                            }
                        )

                        # ===== COMPARISON SET 2: ALTERNATIVE (truncated or extended) =====
                        all_sequences.append(
                            {
                                "gene": gene_name,
                                "transcript_id": transcript_id,
                                "variant_id": feature_id,
                                "sequence": alternative_protein,
                                "length": len(alternative_protein),
                                "variant_type": region_type,  # "truncation" or "extension"
                                "region_type": region_type,
                                "comparison_set": f"{region_type}_analysis",
                                # Mutation fields
                                "mutation_position": None,
                                "mutation_change": None,
                                "aa_change": None,
                                "hgvsc": None,
                                "hgvsp": None,
                                "mutation_impact": None,
                                "mutation_impact_validated": None,
                                "in_alt_start_site": None,
                                "mutation_source": None,
                                "clinvar_variant_id": None,
                            }
                        )

                        # ===== COMPARISON SET 3: MUTATIONS TO APPROPRIATE BASE =====
                        for mut_idx, mut_variant in enumerate(
                            pair["alternative_mutations"]
                        ):
                            mutated_protein = mut_variant["protein"]
                            mut_info = mut_variant["mutation"]

                            # Check length constraints
                            if not (min_length <= len(mutated_protein) <= max_length):
                                continue

                            # Determine correct mutation variant type and base comparison
                            if region_type == "truncation":
                                # TRUNCATION SET: canonical, truncated, canonical+mutations
                                variant_type = "canonical_mutated"
                                base_protein = canonical_protein
                                variant_description = "canonical protein with mutations in truncated region"
                            elif region_type == "extension":
                                # EXTENSION SET: canonical, extended, extended+mutations
                                variant_type = "extension_mutated"
                                base_protein = alternative_protein
                                variant_description = "extended protein with mutations in extension region"
                            else:
                                variant_type = f"{region_type}_mutated"
                                base_protein = canonical_protein
                                variant_description = (
                                    f"{region_type} protein with mutations"
                                )

                            # Skip if identical to base
                            if mutated_protein == base_protein:
                                continue

                            # Calculate AA difference from correct base
                            aa_difference = self._calculate_aa_difference(
                                base_protein, mutated_protein
                            )

                            # Create variant ID
                            base_variant_id = f"{variant_type}_{mut_info['position']}_{mut_info['reference']}>{mut_info['alternate']}"
                            if aa_difference:
                                variant_id = f"{base_variant_id}_{aa_difference}"
                            else:
                                variant_id = base_variant_id

                            all_sequences.append(
                                {
                                    "gene": gene_name,
                                    "transcript_id": transcript_id,
                                    "variant_id": variant_id,
                                    "sequence": mutated_protein,
                                    "length": len(mutated_protein),
                                    "variant_type": variant_type,  # "canonical_mutated" or "extension_mutated"
                                    "region_type": region_type,
                                    "comparison_set": f"{region_type}_analysis",
                                    "variant_description": variant_description,
                                    # Full mutation information
                                    "mutation_position": mut_info["position"],
                                    "mutation_change": f"{mut_info['reference']}>{mut_info['alternate']}",
                                    "aa_change": aa_difference,
                                    "hgvsc": mut_info.get("hgvsc", ""),
                                    "hgvsp": mut_info.get("hgvsp", ""),
                                    "mutation_impact": mut_info.get("impact", ""),
                                    "mutation_impact_validated": mut_info.get(
                                        "impact_validated", ""
                                    ),
                                    "in_alt_start_site": mut_info.get(
                                        "in_alt_start_site", False
                                    ),
                                    "validation_note": mut_info.get(
                                        "validation_note", ""
                                    ),
                                    "applied_to": mut_info.get("applied_to", ""),
                                    "mutation_source": mut_info.get("source", ""),
                                    "clinvar_variant_id": mut_info.get(
                                        "variant_id", ""
                                    ),
                                }
                            )

                else:
                    # Standard extraction without mutations - same 3-way comparison structure
                    gene_pairs = self.extract_gene_proteins(gene_name)
                    if not gene_pairs:
                        skipped_genes.append(gene_name)
                        continue

                    for pair in gene_pairs:
                        canonical_protein = pair["canonical"]["protein"]
                        alternative_protein = pair["alternative"]["protein"]
                        transcript_id = pair["transcript_id"]
                        feature_id = pair["feature_id"]
                        region_type = pair["region_type"]

                        # Check length constraints
                        if not (min_length <= len(canonical_protein) <= max_length):
                            continue
                        if not (min_length <= len(alternative_protein) <= max_length):
                            continue
                        if alternative_protein == canonical_protein:
                            continue

                        # Add canonical sequence
                        all_sequences.append(
                            {
                                "gene": gene_name,
                                "transcript_id": transcript_id,
                                "variant_id": "canonical",
                                "sequence": canonical_protein,
                                "length": len(canonical_protein),
                                "variant_type": "canonical",
                                "region_type": region_type,
                                "comparison_set": f"{region_type}_analysis",
                            }
                        )

                        # Add alternative sequence
                        all_sequences.append(
                            {
                                "gene": gene_name,
                                "transcript_id": transcript_id,
                                "variant_id": feature_id,
                                "sequence": alternative_protein,
                                "length": len(alternative_protein),
                                "variant_type": region_type,
                                "region_type": region_type,
                                "comparison_set": f"{region_type}_analysis",
                            }
                        )

                successful_genes += 1

            except Exception as e:
                # import traceback
                self._debug_print(f"Error processing gene {gene_name}: {str(e)}")
                # traceback.print_exc()  # <— this will show the true source
                skipped_genes.append(gene_name)

        # Create dataset
        dataset = pd.DataFrame(all_sequences)

        # Save dataset if requested
        if not dataset.empty and output_format:
            if "fasta" in output_format.lower():
                self._save_dataset_fasta(dataset, include_mutations)
            if "csv" in output_format.lower():
                self._save_dataset_csv(dataset, include_mutations)

        # Print summary with comparison set breakdown
        self._print_dataset_summary_with_comparison_sets(
            dataset, successful_genes, len(gene_list), skipped_genes, include_mutations
        )

        return dataset

    def create_protein_sequence_dataset_pairs(
        self,
        gene_list: List[str],
        output_format: str = "fasta,csv",
        min_length: int = 50,
        max_length: int = 1000,
    ) -> pd.DataFrame:
        """Create a dataset of protein sequences from paired canonical and alternative transcripts.

        Args:
            gene_list (List[str]): List of gene names to process.
            output_format (str): Format to save sequences ('fasta', 'csv', or both).
            min_length (int): Minimum protein length to include.
            max_length (int): Maximum protein length to include.

        Returns:
            pd.DataFrame: DataFrame with the dataset information.
        """
        all_sequences = []
        total_pairs = 0
        successful_genes = 0
        skipped_genes = []

        # Process each gene
        for gene_idx, gene_name in enumerate(gene_list, 1):
            try:
                # Extract both canonical and alternative proteins
                gene_result = self.extract_gene_proteins(gene_name)

                if not gene_result:
                    skipped_genes.append(gene_name)
                    continue

                for pair in gene_result:
                    canonical_protein = pair["canonical"]["protein"]
                    alternative_protein = pair["alternative"]["protein"]
                    transcript_id = pair["transcript_id"]
                    feature_id = pair["feature_id"]
                    region_type = pair["region_type"]

                    # Check length constraints
                    if not (min_length <= len(canonical_protein) <= max_length):
                        continue
                    if not (min_length <= len(alternative_protein) <= max_length):
                        continue
                    if alternative_protein == canonical_protein:
                        continue

                    # Add canonical sequence
                    all_sequences.append(
                        {
                            "gene": gene_name,
                            "transcript_id": transcript_id,
                            "variant_id": "canonical",
                            "sequence": canonical_protein,
                            "length": len(canonical_protein),
                            "is_alternative": 0,
                            "region_type": "canonical",
                        }
                    )

                    # Add alternative sequence
                    all_sequences.append(
                        {
                            "gene": gene_name,
                            "transcript_id": transcript_id,
                            "variant_id": feature_id,
                            "sequence": alternative_protein,
                            "length": len(alternative_protein),
                            "is_alternative": 1,
                            "region_type": region_type,
                        }
                    )

                    total_pairs += 1

                successful_genes += 1

            except Exception as e:
                skipped_genes.append(gene_name)

        # Create dataset
        dataset = pd.DataFrame(all_sequences)

        # Save dataset
        if not dataset.empty:
            if "fasta" in output_format.lower():
                self._save_dataset_fasta(
                    dataset, include_mutations=False, pairs_only=True
                )
            if "csv" in output_format.lower():
                output_file = self.output_dir / "protein_sequences_pairs.csv"
                dataset.to_csv(output_file, index=False)

        print(
            f"Generated {total_pairs} canonical-alternative pairs from {successful_genes}/{len(gene_list)} genes"
        )
        if skipped_genes:
            print(
                f"Skipped {len(skipped_genes)} genes due to missing data or constraints"
            )

        return dataset

    # ===== VALIDATION METHODS =====

    async def predict_consequence_by_translation(
        self,
        transcript_id: str,
        genomic_pos: int,
        ref_allele: str,
        alt_allele: str,
        current_feature: Optional[pd.Series] = None,
    ) -> str:
        """Predict the functional consequence of a mutation by applying it to the transcript and analyzing the resulting protein change.

        Args:
            transcript_id (str): Transcript identifier to analyze.
            genomic_pos (int): Genomic position of the mutation.
            ref_allele (str): Reference allele at the mutation position.
            alt_allele (str): Alternate allele to introduce at the mutation position.
            current_feature (Optional[pd.Series]): Current alternative feature being analyzed.

        Returns:
            str: Predicted consequence type.
        """
        if self.debug:
            print(
                f"        ├─ Predicting consequence for {ref_allele}>{alt_allele} at {genomic_pos}"
            )
            if current_feature is not None:
                feature_type = current_feature.get("region_type", "unknown")
                feature_range = f"{current_feature.get('start', '?')}-{current_feature.get('end', '?')}"
                print(
                    f"        │  ├─ Context: {feature_type} feature at {feature_range}"
                )

        try:
            # Clean up alleles
            ref_clean = (
                str(ref_allele).strip()
                if ref_allele and str(ref_allele) != "nan"
                else ""
            )
            alt_clean = (
                str(alt_allele).strip()
                if alt_allele and str(alt_allele) != "nan"
                else ""
            )

            if not ref_clean or not alt_clean:
                if self.debug:
                    print(
                        f"        └─ Skipping mutation: empty alleles (ref='{ref_clean}', alt='{alt_clean}')"
                    )
                return "unknown"

            ref_len = len(ref_clean)
            alt_len = len(alt_clean)
            length_change = alt_len - ref_len

            if self.debug:
                print(
                    f"        │  ├─ Allele lengths: ref={ref_len} ({ref_clean}), alt={alt_len} ({alt_clean})"
                )
                print(f"        │  ├─ Length change: {length_change}")

            # Simple length-based classification for indels
            if length_change != 0:
                if length_change % 3 == 0:
                    # In-frame - length change divisible by 3
                    consequence = (
                        "inframe insertion" if length_change > 0 else "inframe deletion"
                    )
                    if self.debug:
                        print(
                            f"        │  └─ In-frame indel: {consequence} ({length_change} bp change)"
                        )
                    return consequence
                else:
                    # Frameshift - length change not divisible by 3
                    if self.debug:
                        print(
                            f"        │  └─ Frameshift variant ({length_change} bp change, not divisible by 3)"
                        )
                    return "frameshift variant"

            # Same length variants - use protein translation for detailed analysis
            if self.debug:
                print(f"        │  ├─ Same length variant - using protein translation")

            # Determine sequence type based on current feature context
            if (
                current_feature is not None
                and current_feature.get("region_type") == "extension"
            ):
                if self.debug:
                    print(f"        │  ├─ Using EXTENSION sequence")

                base_result = self.extract_alternative_protein(
                    transcript_id, current_feature
                )
                if not base_result:
                    if self.debug:
                        print(
                            f"        │  └─ ❌ Could not extract base extension protein"
                        )
                    return "unknown"

                mutated_result = self._apply_mutation_to_sequence(
                    transcript_id,
                    genomic_pos,
                    ref_clean,
                    alt_clean,
                    "extension",
                    current_feature,
                )

            else:
                if self.debug:
                    print(f"        │  ├─ Using CANONICAL sequence")

                base_result = self.extract_canonical_protein(transcript_id)
                if not base_result:
                    if self.debug:
                        print(f"        │  └─ ❌ Could not extract canonical protein")
                    return "unknown"

                mutated_result = self._apply_mutation_to_sequence(
                    transcript_id, genomic_pos, ref_clean, alt_clean, "canonical"
                )

            if self.debug:
                print(f"        │  ├─ Base protein: {len(base_result['protein'])} AA")
                if mutated_result:
                    print(
                        f"        │  ├─ Mutated protein: {len(mutated_result['protein'])} AA"
                    )

            # Classify the change
            if not mutated_result:
                if self.debug:
                    print(f"        │  └─ No protein change detected - synonymous")
                return "synonymous variant"

            consequence = self.classify_protein_change(
                base_result["protein"],
                mutated_result["protein"],
                base_result["coding_sequence"],
                mutated_result["coding_sequence"],
            )

            if self.debug:
                print(f"        │  └─ Protein-based classification: {consequence}")

            return consequence

        except Exception as e:
            if self.debug:
                print(f"        │  └─ ❌ Exception during prediction: {str(e)}")
                import traceback

                print(f"        │     └─ Traceback: {traceback.format_exc()}")

            return "unknown"

    def classify_protein_change(
        self, orig_protein: str, mut_protein: str, orig_cds: str, mut_cds: str
    ) -> str:
        """Classify the type of protein sequence change resulting from a mutation.

        Args:
            orig_protein (str): Original protein sequence before mutation.
            mut_protein (str): Mutated protein sequence after mutation.
            orig_cds (str): Original coding DNA sequence.
            mut_cds (str): Mutated coding DNA sequence.

        Returns:
            str: Consequence type (e.g., 'synonymous variant', 'missense variant', 'nonsense variant', 'frameshift variant', 'inframe deletion', 'inframe insertion').
        """
        if self.debug:
            print(f"          ├─ Classifying protein change:")
            print(f"          │  ├─ Original protein: {len(orig_protein)} AA")
            print(f"          │  ├─ Mutated protein: {len(mut_protein)} AA")
            print(f"          │  ├─ Original CDS: {len(orig_cds)} bp")
            print(f"          │  └─ Mutated CDS: {len(mut_cds)} bp")

        # No protein change
        if mut_protein == orig_protein:
            if self.debug:
                print(
                    f"          └─ Classification: synonymous variant (no protein change)"
                )
            return "synonymous variant"

        # Check for premature stop codon
        if "*" in mut_protein and "*" not in orig_protein[: len(mut_protein)]:
            # Find position of stop codon
            stop_pos = mut_protein.find("*")
            if self.debug:
                print(
                    f"          └─ Classification: nonsense variant (stop codon at position {stop_pos + 1})"
                )
            return "nonsense variant"

        # Check for length changes
        if len(mut_protein) != len(orig_protein):
            cds_length_change = len(mut_cds) - len(orig_cds)
            protein_length_change = len(mut_protein) - len(orig_protein)

            if self.debug:
                print(f"          │  ├─ Length change detected:")
                print(f"          │  │  ├─ CDS change: {cds_length_change} bp")
                print(f"          │  │  └─ Protein change: {protein_length_change} AA")

            if cds_length_change % 3 != 0:
                if self.debug:
                    print(
                        f"          └─ Classification: frameshift variant (CDS change not divisible by 3)"
                    )
                return "frameshift variant"
            else:
                if len(mut_protein) < len(orig_protein):
                    if self.debug:
                        print(
                            f"          └─ Classification: inframe deletion ({-protein_length_change} AA deleted)"
                        )
                    return "inframe deletion"
                else:
                    if self.debug:
                        print(
                            f"          └─ Classification: inframe insertion ({protein_length_change} AA inserted)"
                        )
                    return "inframe insertion"

        # Same length, different sequence = missense
        # Find the position of the change
        if self.debug:
            differences_found = 0
            for i, (orig_aa, mut_aa) in enumerate(zip(orig_protein, mut_protein)):
                if orig_aa != mut_aa:
                    differences_found += 1
                    if differences_found <= 3:  # Show first 3 differences
                        print(
                            f"          │  ├─ AA change at position {i + 1}: {orig_aa}>{mut_aa}"
                        )
                    elif differences_found == 4:
                        print(f"          │  ├─ ... (more differences)")
                        break

            print(f"          │  ├─ Total differences: {differences_found}")
            print(f"          └─ Classification: missense variant")

        return "missense variant"

    async def predict_consequence_fast(
        self,
        transcript_id: str,
        genomic_pos: int,
        ref_allele: str,
        alt_allele: str,
        current_feature: Optional[pd.Series] = None,
    ) -> str:
        """FAST: Predict consequence with caching and simple rules."""
        # Input validation - only exclude truly invalid cases
        if genomic_pos is None:
            if self.debug:
                print(f"        └─ ❌ Invalid genomic position: None")
            return "unknown"

        # Ensure genomic_pos is integer
        try:
            genomic_pos = int(genomic_pos)
        except (ValueError, TypeError):
            if self.debug:
                print(
                    f"        └─ ❌ Invalid genomic position type: {type(genomic_pos)}"
                )
            return "unknown"

        # Check cache first
        cached = self.validation_cache.get_cached_result(
            transcript_id, genomic_pos, ref_allele, alt_allele
        )
        if cached:
            if self.debug:
                print(f"        └─ ⚡ Cached result: {cached}")
            return cached

        try:
            # Handle allele conversion and validation - allow empty strings for indels
            ref_allele = (
                str(ref_allele).strip().upper()
                if ref_allele and str(ref_allele) != "nan"
                else ""
            )
            alt_allele = (
                str(alt_allele).strip().upper()
                if alt_allele and str(alt_allele) != "nan"
                else ""
            )

            # Only exclude if BOTH alleles are invalid/missing
            if (
                (not ref_allele and not alt_allele)
                or ref_allele == "NAN"
                or alt_allele == "NAN"
            ):
                if self.debug:
                    print(
                        f"        └─ ❌ Both alleles invalid: '{ref_allele}', '{alt_allele}'"
                    )
                return "unknown"

            # Quick classification for obvious cases
            ref_len = len(ref_allele)
            alt_len = len(alt_allele)
            length_diff = alt_len - ref_len

            if self.debug:
                print(
                    f"        ├─ Analyzing {ref_allele}>{alt_allele} (lengths: {ref_len}→{alt_len}, diff: {length_diff})"
                )

            # RULE 1: Length changes (including cases where one allele is empty)
            if length_diff != 0:
                if length_diff % 3 == 0:
                    consequence = (
                        "inframe insertion" if length_diff > 0 else "inframe deletion"
                    )
                else:
                    consequence = "frameshift variant"

                self.validation_cache.cache_result(
                    transcript_id, genomic_pos, ref_allele, alt_allele, consequence
                )
                if self.debug:
                    print(f"        ├─ ⚡ Length-based classification: {consequence}")
                    print(
                        f"        │  └─ Reason: {abs(length_diff)} bp change {'in-frame' if length_diff % 3 == 0 else 'causes frameshift'}"
                    )
                return consequence

            # RULE 2: Single BP substitutions - use codon analysis
            if ref_len == 1 and alt_len == 1:
                if self.debug:
                    print(
                        f"        ├─ ⚡ Single BP substitution - analyzing codon context..."
                    )

                consequence = self._analyze_single_bp_fast(
                    transcript_id, genomic_pos, ref_allele, alt_allele, current_feature
                )

                self.validation_cache.cache_result(
                    transcript_id, genomic_pos, ref_allele, alt_allele, consequence
                )
                if self.debug:
                    print(f"        └─ ⚡ Codon-based classification: {consequence}")
                return consequence

            # RULE 3: No change (both empty or identical)
            if ref_allele == alt_allele:
                consequence = "synonymous variant"
                self.validation_cache.cache_result(
                    transcript_id, genomic_pos, ref_allele, alt_allele, consequence
                )
                if self.debug:
                    print(f"        └─ ⚡ No change: {consequence}")
                return consequence

            # RULE 4: Complex cases - fallback to full translation
            if self.debug:
                print(
                    f"        ├─ Complex variant ({ref_len}→{alt_len} bp), using full translation..."
                )

            consequence = await self.predict_consequence_by_translation(
                transcript_id, genomic_pos, ref_allele, alt_allele, current_feature
            )

            self.validation_cache.cache_result(
                transcript_id, genomic_pos, ref_allele, alt_allele, consequence
            )
            if self.debug:
                print(f"        └─ Full translation result: {consequence}")
            return consequence

        except Exception as e:
            if self.debug:
                print(f"        └─ ❌ Fast validation failed: {str(e)}")
                import traceback

                print(f"        └─ Traceback: {traceback.format_exc()}")
            return "unknown"

    def _analyze_single_bp_fast(
        self,
        transcript_id: str,
        genomic_pos: int,
        ref_allele: str,
        alt_allele: str,
        current_feature: Optional[pd.Series] = None,
    ) -> str:
        """FAST: Analyze single BP changes using codon-level analysis with better error handling."""
        try:
            # Get transcript info for strand debugging
            transcript_data = self.genome.get_transcript_features_with_sequence(
                transcript_id
            )
            strand = transcript_data["sequence"]["strand"] if transcript_data else "?"

            if self.debug:
                print(f"        │  ├─ 🧬 MUTATION CLASSIFICATION DEBUG")
                print(f"        │  │  ├─ Transcript: {transcript_id} ({strand} strand)")
                print(
                    f"        │  │  ├─ Genomic mutation: {ref_allele}>{alt_allele} at position {genomic_pos}"
                )

            # Create context-aware cache key
            if current_feature is not None:
                feature_type = current_feature.get("region_type", "canonical")
                cache_key = f"{transcript_id}_{feature_type}"
            else:
                feature_type = "canonical"
                cache_key = f"{transcript_id}_canonical"

            if self.debug:
                print(f"        │  │  ├─ Cache key: {cache_key}")

            # Check if we need to build sequence cache
            if cache_key not in self.validation_cache.coding_sequences:
                if self.debug:
                    print(
                        f"        │  │  ├─ Building {feature_type} sequence cache for {transcript_id}..."
                    )

                # Build sequence cache with error handling
                try:
                    self._build_sequence_cache(transcript_id, current_feature)

                    # Store with context-aware key
                    if transcript_id in self.validation_cache.coding_sequences:
                        self.validation_cache.coding_sequences[cache_key] = (
                            self.validation_cache.coding_sequences[transcript_id]
                        )
                        self.validation_cache.position_maps[cache_key] = (
                            self.validation_cache.position_maps[transcript_id]
                        )
                    else:
                        if self.debug:
                            print(f"        │  │  ├─ ❌ Failed to build sequence cache")
                        return "unknown"

                except Exception as e:
                    if self.debug:
                        print(f"        │  │  ├─ ❌ Error building cache: {str(e)}")
                    return "unknown"
            else:
                if self.debug:
                    print(
                        f"        │  │  ├─ Using cached {feature_type} sequence for {transcript_id}"
                    )

            # Retrieve from context-aware cache
            coding_seq = self.validation_cache.coding_sequences.get(cache_key)
            pos_map = self.validation_cache.position_maps.get(cache_key)

            if not coding_seq or not pos_map:
                if self.debug:
                    print(f"        │  │  ├─ ❌ Empty cache data")
                return "unknown"

            if self.debug:
                print(f"        │  │  ├─ Coding sequence: {len(coding_seq)} bp")
                print(f"        │  │  ├─ Position map: {len(pos_map)} positions")
                if pos_map:
                    min_pos = min(pos_map.keys()) if pos_map.keys() else None
                    max_pos = max(pos_map.keys()) if pos_map.keys() else None
                    print(f"        │  │  ├─ Position range: {min_pos} to {max_pos}")

            # Map genomic position to coding position
            coding_pos = pos_map.get(genomic_pos)

            if self.debug:
                print(
                    f"        │  │  ├─ Position mapping: {genomic_pos} → {coding_pos}"
                )

            if coding_pos is None:
                # Enhanced debugging for mapping failures
                if self.debug:
                    print(f"        │  │  ├─ ❌ Position mapping failed")

                    # Check nearby positions
                    nearby = {
                        k: v for k, v in pos_map.items() if abs(k - genomic_pos) <= 20
                    }
                    if nearby:
                        print(f"        │  │  ├─ Nearby mapped positions:")
                        for k, v in sorted(nearby.items())[:5]:
                            print(f"        │  │  │  ├─ {k} → {v}")

                    # Check if in extension region
                    if (
                        current_feature is not None
                        and current_feature.get("region_type") == "extension"
                    ):
                        ext_start = current_feature.get("start")
                        ext_end = current_feature.get("end")
                        if ext_start is not None and ext_end is not None:
                            if int(ext_start) <= genomic_pos <= int(ext_end):
                                print(
                                    f"        │  │  ├─ Position IS in extension region {ext_start}-{ext_end}"
                                )
                                return "5 prime UTR variant"  # Classify as UTR variant

                return "intronic variant"

            # Validate coding position bounds
            if coding_pos >= len(coding_seq):
                if self.debug:
                    print(
                        f"        │  │  ├─ ❌ BOUNDS ERROR: coding_pos {coding_pos} >= sequence length {len(coding_seq)}"
                    )
                return "unknown"

            # Find the codon
            codon_start = (coding_pos // 3) * 3
            codon_offset = coding_pos % 3
            codon_number = (coding_pos // 3) + 1

            # Extract original codon from coding sequence
            if codon_start + 3 > len(coding_seq):
                if self.debug:
                    print(
                        f"        │  │  ├─ ❌ Incomplete codon at position {codon_start}"
                    )
                return "unknown"

            original_codon = coding_seq[codon_start : codon_start + 3]

            # Handle strand-specific allele conversion
            if strand == "-":
                complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
                transcript_ref = complement_map.get(ref_allele, ref_allele)
                transcript_alt = complement_map.get(alt_allele, alt_allele)
                if self.debug:
                    print(
                        f"        │  │  ├─ Negative strand conversion: {ref_allele}>{alt_allele} → {transcript_ref}>{transcript_alt}"
                    )
            else:
                transcript_ref = ref_allele
                transcript_alt = alt_allele

            # Verify reference matches expectation
            expected_base = original_codon[codon_offset]
            reference_match = expected_base.upper() == transcript_ref.upper()

            if self.debug:
                print(f"        │  │  ├─ Codon analysis:")
                print(f"        │  │  │  ├─ Codon {codon_number}: {original_codon}")
                print(f"        │  │  │  ├─ Position {codon_offset + 1} in codon")
                print(f"        │  │  │  ├─ Expected base: {expected_base}")
                print(f"        │  │  │  ├─ Mutation ref: {transcript_ref}")
                print(f"        │  │  │  └─ Reference match: {reference_match}")

            if not reference_match:
                if self.debug:
                    print(
                        f"        │  │  ├─ ⚠️  Reference mismatch - proceeding with analysis"
                    )

            # Apply mutation to codon
            mutated_codon = (
                original_codon[:codon_offset]
                + transcript_alt
                + original_codon[codon_offset + 1 :]
            )

            # Translate codons
            try:
                original_aa = str(Seq(original_codon).translate())
                mutated_aa = str(Seq(mutated_codon).translate())

                if self.debug:
                    print(f"        │  │  ├─ Translation:")
                    print(f"        │  │  │  ├─ {original_codon} → {original_aa}")
                    print(f"        │  │  │  └─ {mutated_codon} → {mutated_aa}")
            except Exception as e:
                if self.debug:
                    print(f"        │  │  ├─ Translation failed: {e}")
                return "unknown"

            # Classify the change
            if original_aa == mutated_aa:
                consequence = "synonymous variant"
                reason = "no amino acid change"
            elif mutated_aa == "*":
                consequence = "nonsense variant"
                reason = f"introduces stop codon at position {codon_number}"
            elif original_aa == "*":
                consequence = "stop lost variant"
                reason = f"removes stop codon at position {codon_number}"
            else:
                consequence = "missense variant"
                reason = f"amino acid change {original_aa}{codon_number}{mutated_aa}"

            if self.debug:
                print(f"        │  │  └─ FINAL CLASSIFICATION: {consequence}")
                print(f"        │  │     └─ Reason: {reason}")

            return consequence

        except Exception as e:
            if self.debug:
                print(f"        │  │  └─ ❌ Single BP analysis failed: {str(e)}")
                import traceback

                print(f"        │  │     └─ Traceback: {traceback.format_exc()}")
            return "unknown"

    def _build_sequence_cache(
        self,
        transcript_id: str,
        current_feature: Optional[pd.Series] = None,
    ):
        """Build cached coding sequence and position mapping with better error handling."""
        if self.debug:
            print(f"        ├─ Building cache for {transcript_id}")

        try:
            # Determine sequence type
            if (
                current_feature is not None
                and current_feature.get("region_type") == "extension"
            ):
                # For extensions, use alternative sequence
                result = self.extract_alternative_protein(
                    transcript_id, current_feature
                )
                if not result:
                    if self.debug:
                        print(
                            f"        ├─ ❌ Could not extract extension sequence - trying canonical fallback"
                        )
                    # Fallback to canonical if extension fails
                    result = self.extract_canonical_protein(transcript_id)
                    if not result:
                        raise ValueError(
                            f"Could not extract any sequence for {transcript_id}"
                        )

                coding_seq = result["coding_sequence"]

                # Build position map for extension (with fallback)
                try:
                    pos_map = self._build_extension_position_map(
                        transcript_id, current_feature
                    )
                    if not pos_map:
                        if self.debug:
                            print(
                                f"        ├─ ❌ Extension position map empty - using canonical fallback"
                            )
                        pos_map = self._build_canonical_position_map(transcript_id)
                except Exception as e:
                    if self.debug:
                        print(f"        ├─ ❌ Extension position mapping failed: {e}")
                        print(
                            f"        ├─ Using canonical position mapping as fallback"
                        )
                    pos_map = self._build_canonical_position_map(transcript_id)

            else:
                # For canonical sequences
                result = self.extract_canonical_protein(transcript_id)
                if not result:
                    raise ValueError(
                        f"Could not extract canonical sequence for {transcript_id}"
                    )
                coding_seq = result["coding_sequence"]

                # Build position map for canonical
                pos_map = self._build_canonical_position_map(transcript_id)

            # Validate cache contents
            if not coding_seq:
                raise ValueError(f"Empty coding sequence for {transcript_id}")
            if not pos_map:
                if self.debug:
                    print(
                        f"        ├─ ⚠️  WARNING: Empty position map for {transcript_id}"
                    )
                # Create a minimal position map if completely empty
                # This is a fallback that might not be perfect but prevents crashes
                pos_map = {}

            # Cache both
            self.validation_cache.coding_sequences[transcript_id] = coding_seq
            self.validation_cache.position_maps[transcript_id] = pos_map

            if self.debug:
                print(f"        ├─ Cached sequence: {len(coding_seq)} bp")
                print(f"        └─ Cached positions: {len(pos_map)} mapped")

        except Exception as e:
            if self.debug:
                print(f"        └─ ❌ Cache building failed for {transcript_id}: {e}")
            # Set empty cache to prevent repeated failures
            self.validation_cache.coding_sequences[transcript_id] = ""
            self.validation_cache.position_maps[transcript_id] = {}
            raise

    def _build_canonical_position_map(self, transcript_id: str) -> Dict[int, int]:
        """Build genomic position -> coding position map for canonical sequence with better error handling."""
        try:
            features = self.genome.get_transcript_features(transcript_id)
            transcript_data = self.genome.get_transcript_features_with_sequence(
                transcript_id
            )

            if not transcript_data:
                if self.debug:
                    print(f"        │  ├─ ❌ No transcript data for {transcript_id}")
                return {}

            strand = transcript_data["sequence"]["strand"]

            # Get CDS regions and start codon
            cds_regions = features[features["feature_type"] == "CDS"].copy()
            start_codons = features[features["feature_type"] == "start_codon"]

            if cds_regions.empty:
                if self.debug:
                    print(f"        │  ├─ ❌ No CDS regions for {transcript_id}")
                return {}

            if start_codons.empty:
                if self.debug:
                    print(
                        f"        │  ├─ ⚠️  No start codon annotation for {transcript_id}"
                    )
                # Use first CDS as start
                canonical_start = (
                    cds_regions["start"].min()
                    if strand == "+"
                    else cds_regions["end"].max()
                )
            else:
                canonical_start = (
                    start_codons.iloc[0]["start"]
                    if strand == "+"
                    else start_codons.iloc[0]["end"]
                )

            if self.debug:
                print(
                    f"            ├  ├─ Building position map for {transcript_id} ({strand} strand)"
                )
                print(f"            ├  │  ├─ Canonical start: {canonical_start}")
                print(f"            ├  │  ├─ CDS regions: {len(cds_regions)}")

            # Sort CDS regions
            if strand == "+":
                cds_regions = cds_regions.sort_values("start")
            else:
                cds_regions = cds_regions.sort_values("start", ascending=False)

            pos_map = {}
            coding_pos = 0

            for idx, (_, cds) in enumerate(cds_regions.iterrows()):
                cds_start = int(cds["start"])
                cds_end = int(cds["end"])

                if self.debug and idx < 3:  # Show first few CDS regions
                    print(f"        │  │  ├─ CDS {idx + 1}: {cds_start}-{cds_end}")

                # Map positions from canonical start onward
                if strand == "+":
                    effective_start = max(cds_start, canonical_start)
                    effective_end = cds_end
                    if effective_start <= effective_end:
                        for genomic_pos in range(effective_start, effective_end + 1):
                            pos_map[genomic_pos] = coding_pos
                            coding_pos += 1
                else:
                    effective_start = cds_start
                    effective_end = min(cds_end, canonical_start)
                    if effective_start <= effective_end:
                        for genomic_pos in range(
                            effective_end, effective_start - 1, -1
                        ):
                            pos_map[genomic_pos] = coding_pos
                            coding_pos += 1

            if self.debug:
                print(
                    f"        │  │  └─ Final position map: {len(pos_map)} positions mapped"
                )
                # Show a few example mappings
                if pos_map:
                    example_positions = list(pos_map.items())[:3]
                    for genomic_pos, coding_pos in example_positions:
                        print(f"        │  │     ├─ {genomic_pos} → {coding_pos}")

            return pos_map

        except Exception as e:
            if self.debug:
                print(f"        │  └─ ❌ Position mapping failed: {e}")
            return {}

    def _build_extension_position_map(
        self, transcript_id: str, extension_feature: pd.Series
    ) -> Dict[int, int]:
        """Build position map for extension sequences with improved logic."""
        features = self.genome.get_transcript_features(transcript_id)
        transcript_data = self.genome.get_transcript_features_with_sequence(
            transcript_id
        )

        if not transcript_data:
            return {}

        strand = transcript_data["sequence"]["strand"]

        # Get canonical start codon
        start_codons = features[features["feature_type"] == "start_codon"]
        if start_codons.empty:
            if self.debug:
                print(f"        │  ├─ ❌ No start codon found for {transcript_id}")
            return {}

        canonical_start_pos = (
            start_codons.iloc[0]["start"]
            if strand == "+"
            else start_codons.iloc[0]["end"]
        )

        # Get extension boundaries
        extension_start = int(extension_feature["start"])
        extension_end = int(extension_feature["end"])

        if self.debug:
            print(f"        │  ├─ Building EXTENSION position map:")
            print(f"        │  │  ├─ Extension: {extension_start}-{extension_end}")
            print(f"        │  │  ├─ Canonical start: {canonical_start_pos}")
            print(f"        │  │  └─ Strand: {strand}")

        pos_map = {}
        coding_pos = 0

        # STEP 1: Map extension region (5' UTR portions only)
        utr_features = features[features["feature_type"].str.contains("UTR", na=False)]
        extension_utrs = []

        for _, utr in utr_features.iterrows():
            utr_start = int(utr["start"])
            utr_end = int(utr["end"])

            # Find overlap between UTR and extension region, excluding canonical start
            if strand == "+":
                # For + strand: extension is upstream of canonical start
                overlap_start = max(utr_start, extension_start)
                overlap_end = min(utr_end, min(extension_end, canonical_start_pos - 1))
            else:
                # For - strand: extension is downstream of canonical start
                overlap_start = max(
                    utr_start, max(extension_start, canonical_start_pos + 1)
                )
                overlap_end = min(utr_end, extension_end)

            if overlap_start <= overlap_end:
                extension_utrs.append((overlap_start, overlap_end))

        # Sort UTR regions for correct order
        if strand == "+":
            extension_utrs.sort(key=lambda x: x[0])  # 5' to 3'
        else:
            extension_utrs.sort(key=lambda x: x[0], reverse=True)  # 3' to 5'

        # Map extension UTR positions
        for start, end in extension_utrs:
            if strand == "+":
                for genomic_pos in range(start, end + 1):
                    pos_map[genomic_pos] = coding_pos
                    coding_pos += 1
            else:
                for genomic_pos in range(end, start - 1, -1):
                    pos_map[genomic_pos] = coding_pos
                    coding_pos += 1

        extension_length = coding_pos
        if self.debug:
            print(f"        │  │  ├─ Extension sequence mapped: {extension_length} bp")

        # STEP 2: Map CDS regions starting from canonical start
        cds_regions = features[features["feature_type"] == "CDS"].copy()
        if cds_regions.empty:
            if self.debug:
                print(f"        │  │  └─ ❌ No CDS regions found")
            return pos_map

        # Sort CDS regions
        if strand == "+":
            cds_regions = cds_regions.sort_values("start")
        else:
            cds_regions = cds_regions.sort_values("start", ascending=False)

        for _, cds in cds_regions.iterrows():
            cds_start = int(cds["start"])
            cds_end = int(cds["end"])

            # Only include CDS regions from canonical start onward
            if strand == "+":
                effective_start = max(cds_start, canonical_start_pos)
                effective_end = cds_end
                if effective_start <= effective_end:
                    for genomic_pos in range(effective_start, effective_end + 1):
                        pos_map[genomic_pos] = coding_pos
                        coding_pos += 1
            else:
                effective_start = cds_start
                effective_end = min(cds_end, canonical_start_pos)
                if effective_start <= effective_end:
                    for genomic_pos in range(effective_end, effective_start - 1, -1):
                        pos_map[genomic_pos] = coding_pos
                        coding_pos += 1

        total_cds_length = coding_pos - extension_length

        if self.debug:
            print(f"        │  │  ├─ CDS sequence mapped: {total_cds_length} bp")
            print(f"        │  │  └─ Total positions mapped: {len(pos_map)}")

            # Show mapping examples
            if pos_map:
                examples = list(pos_map.items())[:5]
                print(f"        │  │     Examples:")
                for genomic_pos, coding_pos in examples:
                    region_type = (
                        "extension" if coding_pos < extension_length else "CDS"
                    )
                    print(
                        f"        │  │     ├─ {genomic_pos} → {coding_pos} ({region_type})"
                    )

        return pos_map

    # ===== HELPER METHODS =====

    def validate_extension_protein(
        self,
        canonical_protein: str,
        extended_protein: str,
        gene_name: str,
        transcript_id: str,
        verbose: bool = True,
    ) -> bool:
        """Validate that an extension protein makes biological sense.

        This method logs warnings but always returns True to allow sequence generation.
        """
        # Basic checks - these should still fail validation
        if not extended_protein or not canonical_protein:
            if verbose:
                print(f"    ❌ Empty protein sequence for {gene_name}:{transcript_id}")
            return False  # This is a real failure

        # All other checks become warnings but don't fail validation
        warnings_found = False

        # Extended should be longer than canonical
        if len(extended_protein) <= len(canonical_protein):
            if verbose:
                print(
                    f"    ⚠️  WARNING: Extension shorter than canonical for {gene_name}:{transcript_id}"
                )
                print(
                    f"        Canonical: {len(canonical_protein)} AA, Extended: {len(extended_protein)} AA"
                )
            warnings_found = True

        # Check for premature stop codons
        premature_stops = self._count_premature_stops(extended_protein)
        if premature_stops > 0:
            if verbose:
                print(
                    f"    ⚠️  WARNING: Extension has {premature_stops} premature stop codon(s) for {gene_name}:{transcript_id}"
                )
                stop_positions = [
                    i for i, aa in enumerate(extended_protein) if aa == "*"
                ]
                print(f"        Stop positions: {stop_positions}")
            warnings_found = True

        # Check for reasonable extension length
        extension_length = len(extended_protein) - len(canonical_protein)
        if extension_length > 200:
            if verbose:
                print(
                    f"    ⚠️  WARNING: Very long extension ({extension_length} AA) for {gene_name}:{transcript_id}"
                )
            warnings_found = True

        # Check sequence composition
        if not self._check_sequence_composition(
            extended_protein, gene_name, transcript_id, verbose, warn_only=True
        ):
            warnings_found = True

        # Always return True but log if warnings were found
        if warnings_found and verbose:
            print(
                f"    📝 Extension sequence generated despite warnings for {gene_name}:{transcript_id}"
            )

        return True  # Always proceed with sequence generation

    def validate_truncation_protein(
        self,
        canonical_protein: str,
        truncated_protein: str,
        gene_name: str,
        transcript_id: str,
        verbose: bool = True,
    ) -> bool:
        """Validate that a truncation protein makes biological sense.

        This method logs warnings but always returns True to allow sequence generation.
        """
        # Basic checks - these should still fail validation
        if not truncated_protein or not canonical_protein:
            if verbose:
                print(f"    ❌ Empty protein sequence for {gene_name}:{transcript_id}")
            return False  # This is a real failure

        # All other checks become warnings but don't fail validation
        warnings_found = False

        # Truncated should be shorter than canonical
        if len(truncated_protein) >= len(canonical_protein):
            if verbose:
                print(
                    f"    ⚠️  WARNING: Truncation longer than canonical for {gene_name}:{transcript_id}"
                )
                print(
                    f"        Canonical: {len(canonical_protein)} AA, Truncated: {len(truncated_protein)} AA"
                )
            warnings_found = True

        # Check for premature stop codons
        premature_stops = self._count_premature_stops(truncated_protein)
        if premature_stops > 0:
            if verbose:
                print(
                    f"    ⚠️  WARNING: Truncation has {premature_stops} premature stop codon(s) for {gene_name}:{transcript_id}"
                )
                stop_positions = [
                    i for i, aa in enumerate(truncated_protein) if aa == "*"
                ]
                print(f"        Stop positions: {stop_positions}")
            warnings_found = True

        # Check for very short proteins
        if len(truncated_protein) < 20:
            if verbose:
                print(
                    f"    ⚠️  WARNING: Very short truncation ({len(truncated_protein)} AA) for {gene_name}:{transcript_id}"
                )
            warnings_found = True

        # Check sequence composition
        if not self._check_sequence_composition(
            truncated_protein, gene_name, transcript_id, verbose, warn_only=True
        ):
            warnings_found = True

        # Check for extreme truncation
        truncation_length = len(canonical_protein) - len(truncated_protein)
        if truncation_length > len(canonical_protein) * 0.8:
            if verbose:
                print(
                    f"    ⚠️  WARNING: Extreme truncation ({truncation_length} AA removed, {truncation_length / len(canonical_protein) * 100:.1f}%) for {gene_name}:{transcript_id}"
                )
            warnings_found = True

        # Always return True but log if warnings were found
        if warnings_found and verbose:
            print(
                f"    📝 Truncation sequence generated despite warnings for {gene_name}:{transcript_id}"
            )

        return True  # Always proceed with sequence generation

    def _count_premature_stops(self, protein_sequence: str) -> int:
        """Count premature stop codons in protein sequence.

        A stop codon is considered premature if it's not at the very end.

        Args:
            protein_sequence (str): Protein sequence to check.

        Returns:
            int: Number of premature stop codons.
        """
        if not protein_sequence:
            return 0

        # Count all stop codons except the last character (which is expected to be a stop)
        premature_stops = protein_sequence[:-1].count("*")
        return premature_stops

    def _check_sequence_composition(
        self,
        protein_sequence: str,
        gene_name: str,
        transcript_id: str,
        verbose: bool = True,
        warn_only: bool = False,
    ) -> bool:
        """Check if protein sequence composition is reasonable.

        Args:
            protein_sequence (str): Protein sequence to check.
            gene_name (str): Gene name for logging.
            transcript_id (str): Transcript ID for logging.
            verbose (bool): Whether to print warnings.
            warn_only (bool): If True, only warn but don't fail validation.

        Returns:
            bool: True if composition is reasonable, False otherwise.
        """
        if not protein_sequence:
            return False

        warnings_found = False
        sequence_length = len(protein_sequence)

        # Count unusual amino acids - this should still fail
        unusual_count = protein_sequence.count("X")
        if unusual_count > 0:
            if verbose:
                print(
                    f"    ❌ {unusual_count} unknown amino acid(s) (X) in {gene_name}:{transcript_id}"
                )
            if not warn_only:
                return False
            warnings_found = True

        # Check for excessive proline
        proline_count = protein_sequence.count("P")
        proline_percent = (proline_count / sequence_length) * 100
        if proline_percent > 15:
            if verbose:
                print(
                    f"    ⚠️  WARNING: High proline content ({proline_percent:.1f}%) in {gene_name}:{transcript_id}"
                )
            warnings_found = True

        # Check for valid amino acid characters - this should still fail
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY*")
        invalid_aas = set(protein_sequence) - valid_aas
        if invalid_aas:
            if verbose:
                print(
                    f"    ❌ Invalid amino acid characters in {gene_name}:{transcript_id}: {invalid_aas}"
                )
            if not warn_only:
                return False
            warnings_found = True

        return not warnings_found if not warn_only else True

    def validate_protein_pair(
        self,
        canonical_protein: str,
        alternative_protein: str,
        region_type: str,
        gene_name: str,
        transcript_id: str,
        verbose: bool = True,
    ) -> bool:
        """Validate a canonical-alternative protein pair.

        Args:
            canonical_protein (str): Canonical protein sequence.
            alternative_protein (str): Alternative protein sequence.
            region_type (str): Type of alternative region ('extension' or 'truncation').
            gene_name (str): Gene name for logging.
            transcript_id (str): Transcript ID for logging.
            verbose (bool): Whether to print validation warnings.

        Returns:
            bool: True if the pair is valid, False otherwise.
        """
        # Check for identical sequences (no real alternative)
        if canonical_protein == alternative_protein:
            if verbose:
                print(
                    f"    ⚠️  Identical canonical and alternative sequences for {gene_name}:{transcript_id}"
                )
            return False

        # Validate based on region type
        if region_type == "extension":
            return self.validate_extension_protein(
                canonical_protein,
                alternative_protein,
                gene_name,
                transcript_id,
                verbose,
            )
        elif region_type == "truncation":
            return self.validate_truncation_protein(
                canonical_protein,
                alternative_protein,
                gene_name,
                transcript_id,
                verbose,
            )
        else:
            if verbose:
                print(
                    f"    ⚠️  Unknown region type '{region_type}' for {gene_name}:{transcript_id}"
                )
            # For unknown types, just do basic checks
            return (
                alternative_protein.startswith("M")
                and len(alternative_protein) >= 20
                and self._count_premature_stops(alternative_protein) == 0
            )

    def _save_sequences(
        self,
        gene_name: str,
        sequences: Dict[str, Dict[str, str]],
        output_format: str = "fasta",
        exclude_canonical: bool = False,
    ) -> None:
        """Save generated sequences to files.

        Args:
            gene_name (str): Name of the gene.
            sequences (Dict[str, Dict[str, str]]): Dictionary of sequences.
            output_format (str): Output format ('fasta' or 'csv').
            exclude_canonical (bool): Whether to exclude canonical sequences.

        Returns:
            None
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

                # Add alternative sequences
                for variant_id, seq in variants.items():
                    if variant_id == "canonical":
                        continue

                    records.append(
                        SeqRecord(
                            Seq(seq),
                            id=f"{transcript_id}_{variant_id}",
                            description=f"{gene_name} alternative protein variant {variant_id}",
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
                    rows.append(
                        {
                            "gene": gene_name,
                            "transcript_id": transcript_id,
                            "variant_id": "canonical",
                            "sequence": variants["canonical"],
                        }
                    )

                # Add alternative sequences
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

    def _save_dataset_fasta(
        self,
        dataset: pd.DataFrame,
        include_mutations: bool = False,
        pairs_only: bool = False,
    ) -> None:
        """Save dataset as FASTA file.

        Args:
            dataset (pd.DataFrame): DataFrame of sequences.
            include_mutations (bool): Whether to include mutation variants.
            pairs_only (bool): Whether to save only canonical-alternative pairs.

        Returns:
            None
        """
        records = []

        for _, row in dataset.iterrows():
            if include_mutations and "variant_type" in row:
                record_id = f"{row['gene']}_{row['transcript_id']}_{row['variant_id']}"
                description = f"{row['variant_type']} protein"

                if (
                    "canonical_mutated" in str(row.get("variant_type", ""))
                    and "mutation_change" in row
                ):
                    description += f" with mutation {row['mutation_change']}"
            else:
                record_id = f"{row['gene']}_{row['transcript_id']}_{row['variant_id']}"
                is_alt = row.get("is_alternative", 0)
                region_type = row.get("region_type", "unknown")
                description = f"{'Alternative' if is_alt else 'Canonical'} protein ({region_type})"

            records.append(
                SeqRecord(Seq(row["sequence"]), id=record_id, description=description)
            )

        if include_mutations:
            output_file = self.output_dir / "protein_sequences_with_mutations.fasta"
        elif pairs_only:
            output_file = self.output_dir / "protein_sequences_pairs.fasta"
        else:
            output_file = self.output_dir / "protein_sequences.fasta"

        SeqIO.write(records, output_file, "fasta")

    def _save_dataset_csv(
        self, dataset: pd.DataFrame, include_mutations: bool = False
    ) -> None:
        """Save dataset as CSV file.

        Args:
            dataset (pd.DataFrame): DataFrame of sequences.
            include_mutations (bool): Whether to include mutation variants.

        Returns:
            None
        """
        if include_mutations:
            output_file = self.output_dir / "protein_sequences_with_mutations.csv"
        else:
            output_file = self.output_dir / "protein_sequences.csv"

        dataset.to_csv(output_file, index=False)

    def _print_dataset_summary(
        self,
        dataset: pd.DataFrame,
        successful_genes: int,
        total_genes: int,
        skipped_genes: List[str],
        include_mutations: bool = False,
    ) -> None:
        """Print summary of dataset generation.

        Args:
            dataset (pd.DataFrame): DataFrame of sequences.
            successful_genes (int): Number of successfully processed genes.
            total_genes (int): Total number of genes.
            skipped_genes (List[str]): List of skipped genes.
            include_mutations (bool): Whether mutations were included.

        Returns:
            None
        """
        print(f"\nProtein Sequence Generation Summary:")
        print(f"  ├─ Genes processed successfully: {successful_genes}/{total_genes}")

        if not dataset.empty:
            if include_mutations and "variant_type" in dataset.columns:
                canonical_count = len(dataset[dataset["variant_type"] == "canonical"])
                truncation_count = len(dataset[dataset["variant_type"] == "truncation"])
                extension_count = len(dataset[dataset["variant_type"] == "extension"])
                mutated_count = len(
                    dataset[dataset["variant_type"] == "canonical_mutated"]
                )

                print(f"  ├─ Total sequences: {len(dataset)}")
                print(f"  ├─ Canonical sequences: {canonical_count}")
                print(f"  ├─ Truncation sequences: {truncation_count}")
                print(f"  ├─ Extension sequences: {extension_count}")
                print(f"  ├─ Mutated sequences: {mutated_count}")
            else:
                alternative_count = len(dataset[dataset.get("is_alternative", 0) == 1])
                canonical_count = len(dataset[dataset.get("is_alternative", 0) == 0])

                print(f"  ├─ Total sequences: {len(dataset)}")
                print(f"  ├─ Canonical sequences: {canonical_count}")
                print(f"  ├─ Alternative sequences: {alternative_count}")

                # Break down by region type
                if "region_type" in dataset.columns:
                    region_counts = dataset["region_type"].value_counts()
                    print(f"  ├─ By region type:")
                    for region_type, count in region_counts.items():
                        print(f"  │  ├─ {region_type}: {count}")

            print(f"  ├─ Average sequence length: {dataset['length'].mean():.1f}")
            print(
                f"  ├─ Sequence length range: {dataset['length'].min()}-{dataset['length'].max()}"
            )

            genes_with_data = dataset["gene"].nunique()
            print(f"  └─ Genes with valid sequences: {genes_with_data}/{total_genes}")
        else:
            print("  └─ No valid sequences generated")

        if skipped_genes:
            print(
                f"\nSkipped genes ({len(skipped_genes)}): {', '.join(skipped_genes[:10])}"
            )
            if len(skipped_genes) > 10:
                print(f"  ... and {len(skipped_genes) - 10} more")

    def _print_dataset_summary_with_comparison_sets(
        self,
        dataset: pd.DataFrame,
        successful_genes: int,
        total_genes: int,
        skipped_genes: List[str],
        include_mutations: bool = False,
    ) -> None:
        """Print summary with comparison set breakdown."""
        print(f"\nProtein Sequence Generation Summary:")
        print(f"  ├─ Genes processed successfully: {successful_genes}/{total_genes}")

        if not dataset.empty:
            print(f"  ├─ Total sequences: {len(dataset)}")

            # Show comparison sets
            if "comparison_set" in dataset.columns:
                comparison_counts = dataset["comparison_set"].value_counts()
                print(f"  ├─ Comparison sets:")
                for comp_set, count in comparison_counts.items():
                    print(f"  │  ├─ {comp_set}: {count} sequences")

            # Show variant types within each comparison set
            if (
                "variant_type" in dataset.columns
                and "comparison_set" in dataset.columns
            ):
                print(f"  ├─ Variant breakdown by comparison set:")
                for comp_set in dataset["comparison_set"].unique():
                    subset = dataset[dataset["comparison_set"] == comp_set]
                    variant_counts = subset["variant_type"].value_counts()
                    print(f"  │  ├─ {comp_set}:")
                    for variant, count in variant_counts.items():
                        print(f"  │  │  ├─ {variant}: {count}")

            print(f"  ├─ Average sequence length: {dataset['length'].mean():.1f}")
            print(
                f"  ├─ Sequence length range: {dataset['length'].min()}-{dataset['length'].max()}"
            )

            genes_with_data = dataset["gene"].nunique()
            print(f"  └─ Genes with valid sequences: {genes_with_data}/{total_genes}")
        else:
            print("  └─ No valid sequences generated")

        if skipped_genes:
            print(
                f"\nSkipped genes ({len(skipped_genes)}): {', '.join(skipped_genes[:10])}"
            )
            if len(skipped_genes) > 10:
                print(f"  ... and {len(skipped_genes) - 10} more")
