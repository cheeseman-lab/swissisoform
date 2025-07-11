"""Enhanced TruncatedProteinGenerator with validated mutation integration.

This extends the existing TruncatedProteinGenerator to incorporate mutations
found in truncation regions using a validated genomic position approach.
"""

import asyncio
import re
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd
from Bio.Seq import Seq
from swissisoform.translation import TruncatedProteinGenerator
from swissisoform.mutations import MutationHandler


class ValidatedMutationIntegratedProteinGenerator(TruncatedProteinGenerator):
    """Enhanced protein generator with corrected mutation integration.
    
    This class applies mutations to the canonical sequence first, then applies
    truncation logic to generate alternative isoforms from the mutated canonical.
    """
    
    def __init__(self, genome_handler, alt_isoform_handler, output_dir, 
                 mutation_handler=None, debug=False):
        super().__init__(genome_handler, alt_isoform_handler, output_dir)
        self.mutation_handler = mutation_handler or MutationHandler()
        self.debug = debug
        
    def _debug_print(self, message: str):
        """Print debug message if debug mode is enabled."""
        if self.debug:
            print(f"DEBUG: {message}")
    
    def _parse_hgvs_to_alleles(self, hgvsc: str) -> Tuple[Optional[str], Optional[str]]:
        """Parse HGVS coding notation to extract reference and alternate alleles."""
        if not hgvsc or pd.isna(hgvsc) or str(hgvsc) == 'nan':
            return None, None
            
        hgvsc = str(hgvsc).strip()
        
        # Handle simple substitutions: c.76C>T, c.52C>T, etc.
        if '>' in hgvsc:
            try:
                parts = hgvsc.split('>')
                if len(parts) == 2:
                    left_part = parts[0].strip()
                    alt_allele = parts[1].strip()
                    
                    # Extract reference allele: c.76C -> C
                    match = re.search(r'[ATCG]$', left_part, re.IGNORECASE)
                    if match:
                        ref_allele = match.group(0).upper()
                        alt_allele = alt_allele.upper()
                        return ref_allele, alt_allele
            except Exception as e:
                self._debug_print(f"Error parsing substitution HGVS '{hgvsc}': {e}")
                return None, None
                    
        # Handle deletions, insertions, etc. (for future expansion)
        elif 'del' in hgvsc.lower():
            self._debug_print(f"Deletion mutations not yet supported: {hgvsc}")
            return None, None
        elif 'ins' in hgvsc.lower():
            self._debug_print(f"Insertion mutations not yet supported: {hgvsc}")
            return None, None
            
        self._debug_print(f"Unsupported HGVS format: {hgvsc}")
        return None, None
    
    def _parse_protein_change(self, hgvsp: str) -> Optional[str]:
        """Parse protein change notation like 'V10L' to extract change."""
        if not hgvsp or pd.isna(hgvsp):
            return None
            
        # Handle formats like 'V10L', 'p.V10L', etc.
        match = re.search(r'([A-Z*])(\d+)([A-Z*])', str(hgvsp).upper())
        if match:
            from_aa, position, to_aa = match.groups()
            return f"{from_aa}{position}{to_aa}"
        return None
    
    def _compare_proteins(self, protein1: str, protein2: str) -> Optional[str]:
        """Compare two proteins and return the first change found."""
        if len(protein1) != len(protein2):
            return "length_change"
            
        changes = []
        for i, (aa1, aa2) in enumerate(zip(protein1, protein2)):
            if aa1 != aa2:
                changes.append(f"{aa1}{i+1}{aa2}")
        
        if len(changes) == 1:
            return changes[0]
        elif len(changes) == 0:
            return "no_change"
        else:
            return f"multiple_changes_{len(changes)}"
    
    async def _get_mutations_in_canonical_region(self, gene_name: str, 
                                               transcript_id: str,
                                               impact_types: List[str] = None) -> pd.DataFrame:
        """Get mutations within the canonical transcript region.
        
        Args:
            gene_name: Gene name
            transcript_id: Canonical transcript ID
            impact_types: List of impact types to filter by
            
        Returns:
            DataFrame of mutations in the canonical transcript
        """
        self._debug_print(f"Fetching mutations for canonical transcript {transcript_id}")
        
        # Get the full canonical transcript coordinates
        features = self.genome.get_transcript_features(transcript_id)
        if features.empty:
            self._debug_print(f"No features found for transcript {transcript_id}")
            return pd.DataFrame()
            
        # Get the transcript span
        transcript_start = features['start'].min()
        transcript_end = features['end'].max()
        chromosome = features.iloc[0]['chromosome']
        
        self._debug_print(f"Canonical transcript region: {chromosome}:{transcript_start}-{transcript_end}")
        
        # Create transcript features for mutation handler
        transcript_features = pd.DataFrame([{
            'start': transcript_start,
            'end': transcript_end,
            'chromosome': chromosome,
            'name': f'canonical_{transcript_id}'
        }])
        
        # Get mutations in the canonical transcript region
        mutations = await self.mutation_handler.get_visualization_ready_mutations(
            gene_name=gene_name,
            alt_features=transcript_features,
            sources=['clinvar']
        )
        
        if mutations is None or mutations.empty:
            self._debug_print(f"No mutations found in canonical transcript")
            return pd.DataFrame()
        
        self._debug_print(f"Found {len(mutations)} raw mutations in canonical transcript")
        
        # Filter by impact types if specified
        if impact_types:
            before_filter = len(mutations)
            mutations = mutations[mutations['impact'].isin(impact_types)]
            after_filter = len(mutations)
            self._debug_print(f"After impact filtering: {after_filter}/{before_filter} mutations remain")
            
        return mutations
    
    def _apply_mutation_to_canonical_sequence(self, transcript_id: str, 
                                            genomic_pos: int, 
                                            ref_allele: str, 
                                            alt_allele: str) -> Optional[str]:
        """Apply mutation to the canonical transcript sequence.
        
        Args:
            transcript_id: Canonical transcript ID
            genomic_pos: Genomic position of mutation
            ref_allele: Reference allele (HGVS coding strand)
            alt_allele: Alternate allele (HGVS coding strand)
            
        Returns:
            Mutated canonical coding sequence or None if failed
        """
        self._debug_print(f"Applying mutation to canonical transcript {transcript_id}")
        
        # Get transcript data
        transcript_data = self.genome.get_transcript_features_with_sequence(transcript_id)
        if not transcript_data or "sequence" not in transcript_data:
            self._debug_print(f"Could not get transcript sequence data")
            return None

        chromosome = transcript_data["sequence"]["chromosome"]
        strand = transcript_data["sequence"]["strand"]
        
        # Get all CDS regions for canonical transcript
        features = self.genome.get_transcript_features(transcript_id)
        cds_regions = features[features["feature_type"] == "CDS"].copy()
        
        if cds_regions.empty:
            self._debug_print(f"No CDS regions found for canonical transcript")
            return None
        
        # Sort CDS regions by genomic position
        if strand == "+":
            cds_regions = cds_regions.sort_values('start')
        else:
            cds_regions = cds_regions.sort_values('start', ascending=False)
        
        # Convert coding strand alleles to genomic coordinates for mutation application
        if strand == "-":
            complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            genomic_ref_allele = complement_map.get(ref_allele, ref_allele)
            genomic_alt_allele = complement_map.get(alt_allele, alt_allele)
            self._debug_print(f"Negative strand: {ref_allele}>{alt_allele} -> {genomic_ref_allele}>{genomic_alt_allele}")
        else:
            genomic_ref_allele = ref_allele
            genomic_alt_allele = alt_allele
        
        # Validate reference allele at genomic position
        actual_base = self.genome.get_sequence(chromosome, genomic_pos, genomic_pos, "+")
        actual_base = str(actual_base).upper()
        
        if actual_base != genomic_ref_allele.upper():
            self._debug_print(f"Reference mismatch at {genomic_pos}: expected {genomic_ref_allele}, found {actual_base}")
            return None
        
        # Extract canonical coding sequence with mutation applied
        canonical_coding_sequence = ""
        mutation_applied = False
        
        for _, cds in cds_regions.iterrows():
            # Get sequence for this CDS region
            cds_seq = self.genome.get_sequence(chromosome, cds["start"], cds["end"], strand)
            cds_seq = str(cds_seq)
            
            # Check if mutation falls within this CDS region
            if cds["start"] <= genomic_pos <= cds["end"]:
                # Apply mutation at exact genomic position
                if strand == "+":
                    seq_pos = genomic_pos - cds["start"]
                else:
                    seq_pos = cds["end"] - genomic_pos
                
                if 0 <= seq_pos < len(cds_seq):
                    # Apply mutation using genomic coordinates
                    original_base = self.genome.get_sequence(chromosome, genomic_pos, genomic_pos, "+")
                    if str(original_base).upper() == genomic_ref_allele.upper():
                        # For negative strand, we need to apply the complement
                        if strand == "-":
                            # Get the position in the CDS sequence and apply complement of alt allele
                            cds_seq = cds_seq[:seq_pos] + genomic_alt_allele.upper() + cds_seq[seq_pos + 1:]
                        else:
                            cds_seq = cds_seq[:seq_pos] + genomic_alt_allele.upper() + cds_seq[seq_pos + 1:]
                        mutation_applied = True
                        self._debug_print(f"Applied mutation at genomic pos {genomic_pos} (CDS pos {seq_pos})")
                    else:
                        self._debug_print(f"Reference mismatch in CDS sequence")
                        return None
            
            canonical_coding_sequence += cds_seq
        
        if not mutation_applied:
            self._debug_print(f"Warning: mutation not applied to canonical sequence")
            return None
        
        self._debug_print(f"Generated mutated canonical coding sequence length: {len(canonical_coding_sequence)}")
        return canonical_coding_sequence
    
    def _extract_truncated_from_mutated_canonical(self, transcript_id: str,
                                                truncation_feature: pd.Series,
                                                mutated_canonical_sequence: str) -> Optional[Dict]:
        """Extract truncated protein from mutated canonical sequence.
        
        This applies the same truncation logic as extract_truncated_protein,
        but starts from the mutated canonical sequence instead of the wild-type.
        
        Args:
            transcript_id: Transcript ID
            truncation_feature: Truncation feature data
            mutated_canonical_sequence: Mutated canonical coding sequence
            
        Returns:
            Dict with truncated protein from mutated canonical
        """
        self._debug_print(f"Extracting truncated protein from mutated canonical")
        
        # Get features and transcript data
        features = self.genome.get_transcript_features(transcript_id)
        transcript_data = self.genome.get_transcript_features_with_sequence(transcript_id)

        if not transcript_data or "sequence" not in transcript_data:
            return None

        strand = transcript_data["sequence"]["strand"]
        chromosome = transcript_data["sequence"]["chromosome"]

        # Get truncation positions
        truncation_start = truncation_feature["start"]
        truncation_end = truncation_feature["end"]

        # Find the original stop codon
        stop_codons = features[features["feature_type"] == "stop_codon"]
        stop_codon_start = None
        stop_codon_end = None
        if not stop_codons.empty:
            stop_codon_start = stop_codons.iloc[0]["start"]
            stop_codon_end = stop_codons.iloc[0]["end"]

        # Apply truncation logic to determine what part of canonical sequence to use
        if strand == "+":
            alt_start_pos = truncation_end + 1
            extract_end = stop_codon_end if stop_codon_end else None
        else:
            alt_start_pos = truncation_start
            extract_end = stop_codon_start if stop_codon_start else None

        # Get CDS regions that overlap with our truncation range
        cds_regions = features[features["feature_type"] == "CDS"].copy()
        
        # Calculate which part of the mutated canonical sequence corresponds to the truncated region
        if strand == "+":
            canonical_cds_regions = cds_regions.sort_values('start')
        else:
            canonical_cds_regions = cds_regions.sort_values('start', ascending=False)
        
        # Find the positions in the canonical sequence that correspond to our truncated region
        canonical_pos = 0
        truncated_sequence = ""
        
        for _, cds in canonical_cds_regions.iterrows():
            cds_start = cds["start"]
            cds_end = cds["end"]
            cds_length = cds_end - cds_start + 1
            
            # Determine if this CDS region contributes to the truncated region
            if strand == "+":
                # For positive strand, we want sequence from alt_start_pos onward
                if cds_end < alt_start_pos:
                    # This CDS is before our truncation start - skip it
                    canonical_pos += cds_length
                    continue
                if extract_end and cds_start > extract_end:
                    # This CDS is after our truncation end - stop
                    break
                
                # This CDS contributes to truncated region
                effective_start = max(cds_start, alt_start_pos)
                effective_end = min(cds_end, extract_end) if extract_end else cds_end
                
                # Calculate positions in canonical sequence
                start_offset = max(0, effective_start - cds_start)
                end_offset = min(cds_length, effective_end - cds_start + 1)
                
            else:
                # For negative strand, we want sequence from extract_end to alt_start_pos
                if cds_start > alt_start_pos:
                    # This CDS is before our truncation start - skip it
                    canonical_pos += cds_length
                    continue
                if extract_end and cds_end < extract_end:
                    # This CDS is after our truncation end - stop
                    break
                
                # This CDS contributes to truncated region
                effective_start = max(cds_start, extract_end) if extract_end else cds_start
                effective_end = min(cds_end, alt_start_pos)
                
                # Calculate positions in canonical sequence
                start_offset = max(0, effective_start - cds_start)
                end_offset = min(cds_length, effective_end - cds_start + 1)
            
            if end_offset > start_offset:
                # Extract the corresponding part from mutated canonical sequence
                seq_start = canonical_pos + start_offset
                seq_end = canonical_pos + end_offset
                
                if seq_end <= len(mutated_canonical_sequence):
                    truncated_sequence += mutated_canonical_sequence[seq_start:seq_end]
                    self._debug_print(f"Added {seq_end - seq_start} bases from canonical position {seq_start}-{seq_end}")
            
            canonical_pos += cds_length

        # Translate the truncated sequence
        if len(truncated_sequence) >= 3:
            remainder = len(truncated_sequence) % 3
            if remainder > 0:
                truncated_sequence = truncated_sequence[:-remainder]
            protein = str(Seq(truncated_sequence).translate())
        else:
            protein = ""

        self._debug_print(f"Generated truncated protein from mutated canonical, length: {len(protein)}")

        return {
            "coding_sequence": truncated_sequence,
            "protein": protein,
            "strand": strand,
            "transcript_id": transcript_id,
            "truncation_start": truncation_start,
            "truncation_end": truncation_end,
            "alternative_start_pos": alt_start_pos,
        }
    
    def extract_corrected_mutated_protein(self, transcript_id: str, 
                                        truncation_feature: pd.Series,
                                        mutation: pd.Series) -> Optional[Dict]:
        """Extract truncated protein with corrected mutation workflow.
        
        This method follows the corrected biological workflow:
        1. Apply mutation to canonical transcript
        2. Extract truncated protein from mutated canonical
        3. Validate against expected protein change
        
        Args:
            transcript_id: Canonical transcript ID
            truncation_feature: Truncation feature data
            mutation: Mutation data
            
        Returns:
            Dict with corrected mutation results or None if failed
        """
        self._debug_print(f"Starting corrected mutation workflow for {transcript_id}")
        
        # Extract mutation data
        genomic_pos = int(mutation['position'])
        hgvsc = str(mutation.get('hgvsc', ''))
        hgvsp = str(mutation.get('hgvsp', ''))
        
        # Parse HGVS to get alleles
        ref_allele, alt_allele = self._parse_hgvs_to_alleles(hgvsc)
        if not ref_allele or not alt_allele:
            self._debug_print(f"Could not parse HGVS notation: {hgvsc}")
            return None
        
        # Step 1: Apply mutation to canonical transcript
        mutated_canonical_sequence = self._apply_mutation_to_canonical_sequence(
            transcript_id, genomic_pos, ref_allele, alt_allele
        )
        if not mutated_canonical_sequence:
            self._debug_print(f"Failed to apply mutation to canonical sequence")
            return None
        
        # Step 2: Extract truncated protein from mutated canonical
        truncated_result = self._extract_truncated_from_mutated_canonical(
            transcript_id, truncation_feature, mutated_canonical_sequence
        )
        if not truncated_result:
            self._debug_print(f"Failed to extract truncated protein from mutated canonical")
            return None
        
        # Step 3: Generate comparison proteins for validation
        # Get original canonical protein (for comparison)
        original_canonical_result = self.extract_canonical_protein(transcript_id)
        original_truncated_result = self.extract_truncated_protein(transcript_id, truncation_feature)
        
        validation_result = {"reference_match": True}
        
        if original_canonical_result and original_truncated_result:
            # Compare canonical proteins
            canonical_change = self._compare_proteins(
                original_canonical_result['protein'],
                str(Seq(mutated_canonical_sequence).translate()) if len(mutated_canonical_sequence) >= 3 else ""
            )
            
            # Compare truncated proteins
            truncated_change = self._compare_proteins(
                original_truncated_result['protein'],
                truncated_result['protein']
            )
            
            validation_result.update({
                "canonical_protein_change": canonical_change,
                "truncated_protein_change": truncated_change,
            })
            
            # Validate against expected protein change if provided
            if hgvsp and str(hgvsp) != 'nan' and hgvsp.strip():
                expected_protein_change = self._parse_protein_change(hgvsp)
                validation_result.update({
                    "protein_change_expected": expected_protein_change,
                    "canonical_change_match": expected_protein_change == canonical_change,
                    "truncated_change_match": expected_protein_change == truncated_change,
                })
        
        return {
            'coding_sequence': truncated_result['coding_sequence'],
            'protein': truncated_result['protein'],
            'strand': truncated_result['strand'],
            'transcript_id': transcript_id,
            'mutated_canonical_sequence': mutated_canonical_sequence,
            'mutation': {
                'genomic_position': genomic_pos,
                'reference': ref_allele,
                'alternate': alt_allele,
                'hgvsc': hgvsc,
                'hgvsp': hgvsp,
                'impact': mutation.get('impact', ''),
                'variant_id': mutation.get('variant_id', ''),
                'source': mutation.get('source', '')
            },
            'validation': validation_result
        }
    
    async def extract_gene_proteins_with_corrected_mutations(self, gene_name: str, 
                                                           preferred_transcripts: Optional[Set[str]] = None,
                                                           include_mutations: bool = True,
                                                           impact_types: List[str] = None) -> Optional[List[Dict]]:
        """Extract proteins with corrected mutation workflow.
        
        This method applies mutations to canonical transcripts first, then
        generates truncated proteins from the mutated canonical sequences.
        
        Args:
            gene_name: Name of the gene
            preferred_transcripts: Set of preferred transcript IDs
            include_mutations: Whether to include mutation variants
            impact_types: Mutation impact types to include
            
        Returns:
            List of enhanced protein pair dictionaries
        """
        self._debug_print(f"Processing gene {gene_name} with corrected mutation workflow")
        
        # Get base pairs using the same logic as extract_gene_proteins
        base_pairs = self.extract_gene_proteins(gene_name, preferred_transcripts)
        if not base_pairs:
            self._debug_print(f"No base pairs found for {gene_name}")
            return None
            
        self._debug_print(f"Found {len(base_pairs)} base transcript-truncation pairs")
        
        enhanced_pairs = []
        
        for pair in base_pairs:
            enhanced_pair = {
                'gene_name': gene_name,
                'transcript_id': pair['transcript_id'],
                'truncation_id': pair['truncation_id'],
                'canonical': pair['canonical'],
                'truncated_base': pair['truncated'],
                'truncated_mutations': []
            }
            
            if include_mutations:
                self._debug_print(f"Processing mutations for canonical transcript {pair['transcript_id']}")
                
                # Get mutations in the CANONICAL transcript (not the truncation region)
                mutations = await self._get_mutations_in_canonical_region(
                    gene_name, 
                    pair['transcript_id'],
                    impact_types
                )
                
                if not mutations.empty:
                    self._debug_print(f"Found {len(mutations)} mutations in canonical transcript")
                    
                    # Generate corrected mutation variants
                    for _, mutation in mutations.iterrows():
                        try:
                            mutated_result = self.extract_corrected_mutated_protein(
                                pair['transcript_id'], 
                                pair['truncation'],
                                mutation
                            )
                            if mutated_result:
                                enhanced_pair['truncated_mutations'].append(mutated_result)
                                self._debug_print(f"Successfully created corrected mutation variant")
                        except Exception as e:
                            self._debug_print(f"Error creating corrected mutation variant: {e}")
                            continue
                else:
                    self._debug_print(f"No mutations found in canonical transcript")
                    
            enhanced_pairs.append(enhanced_pair)
            
        return enhanced_pairs