"""Enhanced transcript translation and protein sequence generation module.

This module contains the TruncatedProteinGenerator class for generating
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
    """

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
            genome_handler: Initialized GenomeHandler instance
            alt_isoform_handler: Initialized AlternativeIsoform instance
            output_dir: Directory to save output files
            mutation_handler: Optional MutationHandler for mutation integration
            debug: Enable debug mode for detailed output
        """
        self.genome = genome_handler
        self.alt_isoforms = alt_isoform_handler
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.mutation_handler = mutation_handler
        self.debug = debug

    def _debug_print(self, message: str):
        """Print debug message if debug mode is enabled."""
        if self.debug:
            print(f"DEBUG: {message}")

    def extract_canonical_protein(self, transcript_id: str) -> Optional[Dict]:
        """Extract canonical (full) protein sequence using start/stop codon annotations.

        Args:
            transcript_id: Transcript ID to process

        Returns:
            Dict with coding_sequence, protein, and metadata or None if failed
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
            transcript_id: Transcript ID to process
            alternative_feature: Series from get_visualization_features() containing region info

        Returns:
            Dict with coding_sequence, protein, and metadata or None if failed
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
        """Extract full extended protein sequence using hybrid approach.

        For extensions:
        1. Get raw genomic sequence from extension start to canonical start (5'UTR portion)
        2. Get CDS sequence from canonical start to stop codon
        3. Concatenate and translate

        This avoids translating introns/3'UTR while including the extension region.

        Args:
            transcript_id: Transcript ID to process
            extension_feature: Series with extension region information

        Returns:
            Dict with coding_sequence, protein, and metadata or None if failed
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

        # Get extension positions
        extension_start = extension_feature["start"]
        extension_end = extension_feature["end"]

        # Find canonical start and stop codons
        start_codons = features[features["feature_type"] == "start_codon"]
        stop_codons = features[features["feature_type"] == "stop_codon"]

        if start_codons.empty:
            self._debug_print("No start codon found for transcript")
            return None
        if stop_codons.empty:
            self._debug_print("No stop codon found for transcript")
            return None

        canonical_start_codon_start = start_codons.iloc[0]["start"]
        canonical_start_codon_end = start_codons.iloc[0]["end"]
        stop_codon_start = stop_codons.iloc[0]["start"]
        stop_codon_end = stop_codons.iloc[0]["end"]

        # PART 1: Extract extension region (5'UTR portion)
        # This is raw genomic sequence from extension start to canonical start
        if strand == "+":
            extension_region_start = extension_start
            extension_region_end = (
                canonical_start_codon_start - 1
            )  # Stop before canonical start
        else:
            extension_region_start = (
                canonical_start_codon_end  # Start at canonical start
            )
            extension_region_end = extension_end

        self._debug_print(
            f"Extension region: {extension_region_start}-{extension_region_end} ({strand} strand)"
        )

        try:
            extension_sequence = self.genome.get_sequence(
                chromosome, extension_region_start, extension_region_end, strand
            )
            extension_sequence = str(extension_sequence)
            self._debug_print(
                f"Extension sequence length: {len(extension_sequence)} bp"
            )
        except Exception as e:
            self._debug_print(f"Error extracting extension sequence: {e}")
            return None

        # PART 2: Extract CDS regions from canonical start to stop codon
        # This is the same logic as canonical protein extraction
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
            cds_start = cds["start"]
            cds_end = cds["end"]

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
            cds_seq = self.genome.get_sequence(
                chromosome, cds["start"], cds["end"], strand
            )
            cds_sequence += str(cds_seq)

        self._debug_print(f"CDS sequence length: {len(cds_sequence)} bp")
        self._debug_print(f"Total CDS regions used: {len(overlapping_cds)}")

        # PART 3: Combine extension + CDS sequences
        full_coding_sequence = extension_sequence + cds_sequence
        self._debug_print(f"Combined sequence length: {len(full_coding_sequence)} bp")

        # Translate the combined sequence
        if len(full_coding_sequence) >= 3:
            # Ensure length is divisible by 3 for clean translation
            remainder = len(full_coding_sequence) % 3
            if remainder > 0:
                full_coding_sequence = full_coding_sequence[:-remainder]

            protein = str(Seq(full_coding_sequence).translate())
            self._debug_print(f"Extended protein length: {len(protein)} AA")

            # Basic validation
            if not protein.startswith("M"):
                self._debug_print(
                    f"Warning: Extended protein doesn't start with M, starts with: {protein[:5]}"
                )
        else:
            self._debug_print("Combined sequence too short for translation")
            protein = ""

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
            "extraction_method": "hybrid_5utr_plus_cds",
            "region_type": "extension",
        }

    def _extract_truncation_protein(
        self, transcript_id: str, truncation_feature: pd.Series
    ) -> Optional[Dict]:
        """Extract protein sequence starting from the base after truncation.

        Args:
            transcript_id: Transcript ID to process
            truncation_feature: Series with truncation region information

        Returns:
            Dict with coding_sequence, protein, and metadata or None if failed
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

                effective_start = (
                    max(cds_start, extract_end) if extract_end else cds_start
                )
                effective_end = min(cds_end, alt_start_pos)

            # Handle the length check with special case for single nucleotide CDS
            if effective_end > effective_start:
                # Normal case - include the region
                overlapping_cds.append(
                    {
                        "start": effective_start,
                        "end": effective_end,
                        "length": effective_end - effective_start + 1,
                    }
                )
            elif effective_end == effective_start:
                # Boundary case - check if this is a legitimate single nucleotide CDS
                if cds_start == cds_end:
                    # Legitimate single nucleotide CDS - include it
                    overlapping_cds.append(
                        {
                            "start": effective_start,
                            "end": effective_end,
                            "length": 1,
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

        # Translate the sequence directly
        if len(coding_sequence) >= 3:
            # Ensure length is divisible by 3
            remainder = len(coding_sequence) % 3
            if remainder > 0:
                coding_sequence = coding_sequence[
                    remainder:
                ]  # Remove from beginning due to BED intersect issues

            protein = str(Seq(coding_sequence).translate())
        else:
            protein = ""

        return {
            "coding_sequence": coding_sequence,
            "protein": protein,
            "strand": strand,
            "transcript_id": transcript_id,
            "truncation_start": truncation_start,
            "truncation_end": truncation_end,
            "alternative_start_pos": alt_start_pos,
            "total_cds_regions": len(overlapping_cds),
            "region_type": "truncation",
        }

    def extract_gene_proteins(self, gene_name: str) -> Optional[List[Dict]]:
        """Extract canonical and alternative proteins for all features of a gene.

        Args:
            gene_name: Name of the gene

        Returns:
            List of dicts with canonical and alternative results for each feature
        """
        self._debug_print(f"Processing gene: {gene_name}")

        # Get alternative features - transcript IDs are already embedded
        alt_features = self.alt_isoforms.get_visualization_features(gene_name)
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

    def _parse_hgvs_to_alleles(self, hgvsc: str) -> Tuple[Optional[str], Optional[str]]:
        """Parse HGVS coding notation to extract reference and alternate alleles."""
        if not hgvsc or pd.isna(hgvsc) or str(hgvsc) == "nan":
            return None, None

        hgvsc = str(hgvsc).strip()

        # Handle simple substitutions: c.76C>T, c.52C>T, etc.
        if ">" in hgvsc:
            try:
                parts = hgvsc.split(">")
                if len(parts) == 2:
                    left_part = parts[0].strip()
                    alt_allele = parts[1].strip()

                    # Extract reference allele: c.76C -> C
                    match = re.search(r"[ATCG]$", left_part, re.IGNORECASE)
                    if match:
                        ref_allele = match.group(0).upper()
                        alt_allele = alt_allele.upper()
                        return ref_allele, alt_allele
            except Exception as e:
                self._debug_print(f"Error parsing substitution HGVS '{hgvsc}': {e}")
                return None, None

        # Handle deletions, insertions, etc. (for future expansion)
        elif "del" in hgvsc.lower():
            self._debug_print(f"Deletion mutations not yet supported: {hgvsc}")
            return None, None
        elif "ins" in hgvsc.lower():
            self._debug_print(f"Insertion mutations not yet supported: {hgvsc}")
            return None, None

        self._debug_print(f"Could not parse HGVS notation: {hgvsc}")
        return None, None

    def _apply_mutation_to_canonical_sequence(
        self, transcript_id: str, genomic_pos: int, ref_allele: str, alt_allele: str
    ) -> Optional[str]:
        """Apply mutation to canonical coding sequence using validated genomic position."""
        self._debug_print(
            f"Applying mutation {ref_allele}>{alt_allele} at position {genomic_pos}"
        )

        # Get the original canonical protein result
        canonical_result = self.extract_canonical_protein(transcript_id)
        if not canonical_result:
            self._debug_print(f"Could not extract canonical protein")
            return None

        original_coding_sequence = canonical_result["coding_sequence"]
        strand = canonical_result["strand"]

        # Get transcript features to map genomic position to coding sequence position
        features = self.genome.get_transcript_features(transcript_id)
        transcript_data = self.genome.get_transcript_features_with_sequence(
            transcript_id
        )

        if not transcript_data:
            return None

        chromosome = transcript_data["sequence"]["chromosome"]

        # Validate reference allele at genomic position
        actual_base = self.genome.get_sequence(
            chromosome, genomic_pos, genomic_pos, "+"
        )
        actual_base = str(actual_base).upper()

        # Convert HGVS alleles to genomic coordinates for validation
        if strand == "-":
            complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
            genomic_ref_allele = complement_map.get(ref_allele, ref_allele)
            genomic_alt_allele = complement_map.get(alt_allele, alt_allele)
        else:
            genomic_ref_allele = ref_allele
            genomic_alt_allele = alt_allele

        if actual_base != genomic_ref_allele.upper():
            self._debug_print(
                f"Reference mismatch at {genomic_pos}: expected {genomic_ref_allele}, found {actual_base}"
            )
            return None

        # Map genomic position to coding sequence position
        start_codons = features[features["feature_type"] == "start_codon"]
        if start_codons.empty:
            return None

        start_codon_start = start_codons.iloc[0]["start"]
        cds_regions = features[features["feature_type"] == "CDS"].copy()

        if cds_regions.empty:
            return None

        # Sort CDS regions
        if strand == "+":
            cds_regions = cds_regions.sort_values("start")
        else:
            cds_regions = cds_regions.sort_values("start", ascending=False)

        # Find which CDS contains our mutation and the position within the coding sequence
        coding_pos = 0
        found_position = False

        for _, cds in cds_regions.iterrows():
            if cds["start"] <= genomic_pos <= cds["end"]:
                # Found the CDS containing our mutation
                if strand == "+":
                    # For positive strand
                    if (
                        cds["start"] >= start_codon_start
                    ):  # Only count CDS after start codon
                        offset_in_cds = genomic_pos - cds["start"]
                        final_coding_pos = coding_pos + offset_in_cds
                        found_position = True
                        break
                else:
                    # For negative strand
                    offset_from_end = cds["end"] - genomic_pos
                    final_coding_pos = coding_pos + offset_from_end
                    found_position = True
                    break

            # Add length of this CDS to coding position counter
            if strand == "+":
                if cds["end"] >= start_codon_start:
                    if cds["start"] < start_codon_start:
                        coding_pos += cds["end"] - start_codon_start + 1
                    else:
                        coding_pos += cds["end"] - cds["start"] + 1
            else:
                coding_pos += cds["end"] - cds["start"] + 1

        if not found_position:
            self._debug_print(
                f"Could not map genomic position {genomic_pos} to coding sequence position"
            )
            return None

        self._debug_print(
            f"Mapped genomic position {genomic_pos} to coding position {final_coding_pos}"
        )

        # Apply the mutation to the coding sequence
        if final_coding_pos >= len(original_coding_sequence):
            self._debug_print(
                f"Coding position {final_coding_pos} is beyond sequence length {len(original_coding_sequence)}"
            )
            return None

        # Verify the reference matches what we expect in the coding sequence
        current_base = original_coding_sequence[final_coding_pos]

        if current_base.upper() != ref_allele.upper():
            self._debug_print(
                f"Reference mismatch in coding sequence at pos {final_coding_pos}: expected {ref_allele}, found {current_base}"
            )
            return None

        # Apply the mutation
        mutated_coding_sequence = (
            original_coding_sequence[:final_coding_pos]
            + alt_allele.upper()
            + original_coding_sequence[final_coding_pos + 1 :]
        )

        self._debug_print(
            f"Applied mutation {ref_allele}>{alt_allele} at coding position {final_coding_pos}"
        )

        return mutated_coding_sequence

    async def _get_mutations_in_region(
        self, gene_name: str, start: int, end: int, impact_types: List[str] = None
    ) -> pd.DataFrame:
        """Get mutations within a specific genomic region."""
        self._debug_print(f"Fetching mutations for {gene_name} in region {start}-{end}")

        if not self.mutation_handler:
            self._debug_print("No mutation handler available")
            return pd.DataFrame()

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

        # Get mutations in this region
        mutations = await self.mutation_handler.get_visualization_ready_mutations(
            gene_name=gene_name, alt_features=region_features, sources=["clinvar"]
        )

        if mutations is None or mutations.empty:
            self._debug_print(f"No mutations found in region")
            return pd.DataFrame()

        # Filter to only missense variants for now
        mutations = mutations[mutations["impact"] == "missense variant"].copy()
        # Filter by impact types if specified (but only missense variants will remain)
        if impact_types:
            mutations = mutations[mutations["impact"].isin(impact_types)]

        self._debug_print(
            f"Class is only equipped to handle missense variants, temporarily."
        )
        return mutations

    def _calculate_aa_difference(
        self, original_protein: str, mutated_protein: str
    ) -> Optional[str]:
        """Calculate amino acid difference between two protein sequences.

        Args:
            original_protein: Original protein sequence
            mutated_protein: Mutated protein sequence

        Returns:
            String in format "A123V" (original AA, position, new AA) or None if no single change
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
        impact_types: List[str] = None,
    ) -> Optional[List[Dict]]:
        """Extract proteins for all alternative features with optional mutation integration.

        Args:
            gene_name: Name of the gene
            include_mutations: Whether to include mutation variants
            impact_types: Mutation impact types to include

        Returns:
            List of enhanced protein pair dictionaries
        """
        self._debug_print(
            f"Processing gene {gene_name} with mutations={include_mutations}"
        )

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
                    gene_name,
                    feature_info["start"],
                    feature_info["end"],
                    impact_types,
                )

                if not mutations.empty:
                    self._debug_print(
                        f"Found {len(mutations)} mutations in alternative region"
                    )

                    # Generate mutation variants applied to canonical sequence
                    for _, mutation in mutations.iterrows():
                        try:
                            # Parse mutation
                            genomic_pos = int(mutation["position"])
                            hgvsc = str(mutation.get("hgvsc", ""))
                            hgvsp = str(mutation.get("hgvsp", ""))

                            # Parse HGVS to get alleles
                            ref_allele, alt_allele = self._parse_hgvs_to_alleles(hgvsc)
                            if not ref_allele or not alt_allele:
                                continue

                            # Apply mutation to canonical sequence
                            mutated_coding_sequence = (
                                self._apply_mutation_to_canonical_sequence(
                                    pair["transcript_id"],
                                    genomic_pos,
                                    ref_allele,
                                    alt_allele,
                                )
                            )

                            if mutated_coding_sequence:
                                # Translate mutated sequence
                                if len(mutated_coding_sequence) >= 3:
                                    mutated_protein = str(
                                        Seq(mutated_coding_sequence).translate()
                                    )
                                else:
                                    mutated_protein = ""

                                mutation_result = {
                                    "coding_sequence": mutated_coding_sequence,
                                    "protein": mutated_protein,
                                    "transcript_id": pair["transcript_id"],
                                    "mutation": {
                                        "position": genomic_pos,
                                        "reference": ref_allele,
                                        "alternate": alt_allele,
                                        "hgvsc": hgvsc,
                                        "hgvsp": hgvsp,
                                        "impact": mutation.get("impact", ""),
                                        "variant_id": mutation.get("variant_id", ""),
                                        "source": mutation.get("source", ""),
                                    },
                                }

                                enhanced_pair["alternative_mutations"].append(
                                    mutation_result
                                )
                                self._debug_print(
                                    f"Successfully created mutation variant"
                                )

                        except Exception as e:
                            self._debug_print(f"Error creating mutation variant: {e}")
                            continue

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
            gene_name: Name of the gene
            output_format: Format to save sequences ('fasta', 'csv', or None for no save)
            exclude_canonical: Whether to exclude canonical transcripts

        Returns:
            Dictionary with transcript IDs as keys, containing both canonical and
            alternative protein sequences for each transcript
        """
        result = {}

        try:
            # Get alternative features - transcript IDs are embedded
            alt_features = self.alt_isoforms.get_visualization_features(gene_name)
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
        impact_types: List[str] = None,
        output_format: str = "fasta,csv",
        min_length: int = 10,
        max_length: int = 100000,
    ) -> pd.DataFrame:
        """Create a dataset of protein sequences with optional mutation integration."""
        all_sequences = []
        successful_genes = 0
        skipped_genes = []

        if impact_types is None:
            impact_types = ["missense variant"]

        # Process each gene
        for gene_idx, gene_name in enumerate(gene_list, 1):
            try:
                self._debug_print(
                    f"Processing gene {gene_idx}/{len(gene_list)}: {gene_name}"
                )

                if include_mutations and self.mutation_handler:
                    # Use mutation-enhanced extraction
                    enhanced_pairs = await self.extract_gene_proteins_with_mutations(
                        gene_name,
                        include_mutations,
                        impact_types,
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
                                # Empty mutation fields for canonical
                                "mutation_position": None,
                                "mutation_change": None,
                                "aa_change": None,
                                "hgvsc": None,
                                "hgvsp": None,
                                "mutation_impact": None,
                                "mutation_source": None,
                                "clinvar_variant_id": None,
                            }
                        )

                        # Add base alternative sequence
                        all_sequences.append(
                            {
                                "gene": gene_name,
                                "transcript_id": transcript_id,
                                "variant_id": feature_id,
                                "sequence": alternative_protein,
                                "length": len(alternative_protein),
                                "variant_type": region_type,
                                "region_type": region_type,
                                # Empty mutation fields for base alternative
                                "mutation_position": None,
                                "mutation_change": None,
                                "aa_change": None,
                                "hgvsc": None,
                                "hgvsp": None,
                                "mutation_impact": None,
                                "mutation_source": None,
                                "clinvar_variant_id": None,
                            }
                        )

                        # Add mutation variants with full information
                        for mut_idx, mut_variant in enumerate(
                            pair["alternative_mutations"]
                        ):
                            mutated_protein = mut_variant["protein"]
                            if (
                                min_length <= len(mutated_protein) <= max_length
                                and mutated_protein != canonical_protein
                            ):
                                mut_info = mut_variant["mutation"]

                                # Calculate amino acid difference
                                aa_difference = self._calculate_aa_difference(
                                    canonical_protein, mutated_protein
                                )

                                # Create enhanced variant ID with AA change
                                base_variant_id = f"canonical_mut_{mut_info['position']}_{mut_info['reference']}>{mut_info['alternate']}"
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
                                        "variant_type": "canonical_mutated",
                                        "region_type": region_type,
                                        # Full mutation information
                                        "mutation_position": mut_info["position"],
                                        "mutation_change": f"{mut_info['reference']}>{mut_info['alternate']}",
                                        "aa_change": aa_difference,
                                        "hgvsc": mut_info.get("hgvsc", ""),
                                        "hgvsp": mut_info.get("hgvsp", ""),
                                        "mutation_impact": mut_info.get("impact", ""),
                                        "mutation_source": mut_info.get("source", ""),
                                        "clinvar_variant_id": mut_info.get(
                                            "variant_id", ""
                                        ),
                                    }
                                )
                else:
                    # Standard extraction without mutations
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
                            }
                        )

                successful_genes += 1
            except Exception as e:
                self._debug_print(f"Error processing gene {gene_name}: {str(e)}")
                skipped_genes.append(gene_name)

        # Create dataset
        dataset = pd.DataFrame(all_sequences)

        # Save dataset if requested
        if not dataset.empty and output_format:
            if "fasta" in output_format.lower():
                self._save_dataset_fasta(dataset, include_mutations)
            if "csv" in output_format.lower():
                self._save_dataset_csv(dataset, include_mutations)

        # Print summary
        self._print_dataset_summary(
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
            gene_list: List of gene names to process
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

    # ===== HELPER METHODS =====

    def _save_sequences(
        self,
        gene_name: str,
        sequences: Dict[str, Dict[str, str]],
        output_format: str = "fasta",
        exclude_canonical: bool = False,
    ) -> None:
        """Save generated sequences to files."""
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
        """Save dataset as FASTA file."""
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
        """Save dataset as CSV file."""
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
        """Print summary of dataset generation."""
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
