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
    """Enhanced protein generator with validated mutation integration.
    
    This class follows the same patterns as TruncatedProteinGenerator but adds
    validated mutation integration using genomic position validation.
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
        """Parse HGVS coding notation to extract reference and alternate alleles.
        
        Args:
            hgvsc: HGVS coding sequence notation (e.g., 'c.76C>T')
            
        Returns:
            Tuple of (reference_allele, alternate_allele) or (None, None) if parsing fails
        """
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
            
        self._debug_print(f"Could not parse HGVS notation: {hgvsc}")
        return None, None
    
    def _parse_protein_change(self, hgvsp: str) -> Optional[str]:
        """Parse protein change notation like 'V10L' to extract change.
        
        Args:
            hgvsp: Protein change notation
            
        Returns:
            Standardized protein change string or None
        """
        if not hgvsp or pd.isna(hgvsp):
            return None
            
        # Handle formats like 'V10L', 'p.V10L', etc.
        match = re.search(r'([A-Z*])(\d+)([A-Z*])', str(hgvsp).upper())
        if match:
            from_aa, position, to_aa = match.groups()
            return f"{from_aa}{position}{to_aa}"
        return None
    
    def _compare_proteins(self, protein1: str, protein2: str) -> Optional[str]:
        """Compare two proteins and return the first change found.
        
        Args:
            protein1: Original protein sequence
            protein2: Modified protein sequence
            
        Returns:
            Change notation (e.g., 'V10L') or special codes
        """
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
    
    async def _get_mutations_in_region(self, gene_name: str, start: int, end: int, 
                                      impact_types: List[str] = None) -> pd.DataFrame:
        """Get mutations within a specific genomic region.
        
        Args:
            gene_name: Gene name
            start: Start position of region
            end: End position of region
            impact_types: List of impact types to filter by
            
        Returns:
            DataFrame of mutations in the region
        """
        self._debug_print(f"Fetching mutations for {gene_name} in region {start}-{end}")
        self._debug_print(f"Looking for impact types: {impact_types}")
        
        # Create region features for mutation handler
        region_features = pd.DataFrame([{
            'start': start,
            'end': end,
            'chromosome': 'chr1',  # Will be corrected by mutation handler
            'name': f'region_{start}_{end}'
        }])
        
        # Get mutations in this region
        mutations = await self.mutation_handler.get_visualization_ready_mutations(
            gene_name=gene_name,
            alt_features=region_features,
            sources=['clinvar']
        )
        
        if mutations is None or mutations.empty:
            self._debug_print(f"No mutations found in region")
            return pd.DataFrame()
        
        self._debug_print(f"Found {len(mutations)} raw mutations")
        
        # Filter by impact types if specified
        if impact_types:
            before_filter = len(mutations)
            mutations = mutations[mutations['impact'].isin(impact_types)]
            after_filter = len(mutations)
            self._debug_print(f"After impact filtering: {after_filter}/{before_filter} mutations remain")
            
        return mutations
    
    def extract_validated_mutated_protein(self, transcript_id: str, 
                                        mutation: pd.Series,
                                        apply_to_canonical: bool = True) -> Optional[Dict]:
        """Extract canonical protein with validated mutation applied.
        
        This method follows the validated approach:
        1. Validate reference allele at genomic position
        2. Apply mutation at exact genomic position to canonical sequence
        3. Validate protein change against expected HGVSP
        
        Args:
            transcript_id: Transcript ID to process
            mutation: Mutation data with position, HGVS notation
            apply_to_canonical: If True, apply to canonical sequence; if False, apply to truncated
            
        Returns:
            Dict with validated mutation results or None if validation fails
        """
        self._debug_print(f"Starting validated mutation for transcript {transcript_id}")
        
        # Extract mutation data
        genomic_pos = int(mutation['position'])
        hgvsc = str(mutation.get('hgvsc', ''))
        hgvsp = str(mutation.get('hgvsp', ''))
        
        self._debug_print(f"Mutation: pos={genomic_pos}, HGVSC={hgvsc}, HGVSP={hgvsp}")
        
        # Step 1: Parse HGVS to get expected reference and alternate
        ref_allele, alt_allele = self._parse_hgvs_to_alleles(hgvsc)
        if not ref_allele or not alt_allele:
            self._debug_print(f"Could not parse HGVS notation: {hgvsc}")
            return None
            
        self._debug_print(f"Parsed alleles: {ref_allele}>{alt_allele}")
        
        # Get transcript features
        features = self.genome.get_transcript_features(transcript_id)
        if features.empty:
            self._debug_print(f"No features found for transcript")
            return None
            
        chromosome = features.iloc[0]['chromosome']
        strand = features.iloc[0]['strand']
        
        # Step 2: Validate reference allele at genomic position
        # CRITICAL: For negative strand genes, HGVS refers to the coding strand,
        # but genomic coordinates are on the reference (+) strand
        actual_base = self.genome.get_sequence(chromosome, genomic_pos, genomic_pos, "+")
        actual_base = str(actual_base).upper()
        
        # For negative strand genes, we need to compare against the complement
        if strand == "-":
            complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            expected_genomic_base = complement_map.get(ref_allele, ref_allele)
            self._debug_print(f"Negative strand: HGVS ref={ref_allele} -> genomic ref={expected_genomic_base}")
        else:
            expected_genomic_base = ref_allele
        
        self._debug_print(f"Reference validation: expected={expected_genomic_base}, found={actual_base}")
        
        if actual_base != expected_genomic_base:
            self._debug_print(f"Reference mismatch - skipping mutation")
            return None
        
        # Step 3: Generate original canonical protein (for comparison)
        original_result = self.extract_canonical_protein(transcript_id)
        if not original_result:
            self._debug_print(f"Could not extract original canonical protein")
            return None
            
        original_protein = original_result['protein']
        
        # Step 4: Generate mutated canonical protein using SIMPLE APPROACH
        mutated_coding_sequence = self._apply_simple_mutation_to_canonical_coding_sequence(
            transcript_id, genomic_pos, ref_allele, alt_allele
        )
        if not mutated_coding_sequence:
            self._debug_print(f"Could not apply mutation to canonical sequence")
            return None
        
        # Translate the mutated coding sequence
        if len(mutated_coding_sequence) >= 3:
            mutated_protein = str(Seq(mutated_coding_sequence).translate())
        else:
            mutated_protein = ""
        
        # Step 5: Validate protein change (if HGVSP provided)
        validation_result = {"reference_match": True}
        
        if hgvsp and str(hgvsp) != 'nan' and hgvsp.strip():
            expected_protein_change = self._parse_protein_change(hgvsp)
            actual_protein_change = self._compare_proteins(original_protein, mutated_protein)
            
            self._debug_print(f"Protein change: expected={expected_protein_change}, actual={actual_protein_change}")
            
            validation_result.update({
                "protein_change_expected": expected_protein_change,
                "protein_change_actual": actual_protein_change,
                "protein_change_match": expected_protein_change == actual_protein_change
            })
        
        return {
            'coding_sequence': mutated_coding_sequence,
            'protein': mutated_protein,
            'strand': strand,
            'transcript_id': transcript_id,
            'mutation': {
                'genomic_position': genomic_pos,
                'reference': ref_allele,  # Keep original HGVS reference
                'alternate': alt_allele,  # Keep original HGVS alternate
                'genomic_reference': expected_genomic_base,  # Add genomic reference
                'genomic_alternate': alt_allele,  # Keep simple for now
                'hgvsc': hgvsc,
                'hgvsp': hgvsp,
                'impact': mutation.get('impact', ''),
                'variant_id': mutation.get('variant_id', ''),
                'source': mutation.get('source', ''),
                'position': genomic_pos
            },
            'validation': validation_result
        }
    
    def _apply_simple_mutation_to_canonical_coding_sequence(self, transcript_id: str, 
                                                          genomic_pos: int, 
                                                          ref_allele: str, 
                                                          alt_allele: str) -> Optional[str]:
        """Apply mutation using a simple approach: modify the coding sequence directly.
        
        This avoids complex coordinate mapping by working with the extracted coding sequence.
        """
        self._debug_print(f"Applying simple mutation {ref_allele}>{alt_allele} at position {genomic_pos}")
        
        # Get the original canonical protein result
        canonical_result = self.extract_canonical_protein(transcript_id)
        if not canonical_result:
            self._debug_print(f"Could not extract canonical protein")
            return None
        
        original_coding_sequence = canonical_result['coding_sequence']
        strand = canonical_result['strand']
        
        # Get transcript features to map genomic position to coding sequence position
        features = self.genome.get_transcript_features(transcript_id)
        transcript_data = self.genome.get_transcript_features_with_sequence(transcript_id)
        
        if not transcript_data:
            return None
        
        chromosome = transcript_data["sequence"]["chromosome"]
        
        # Validate reference allele at genomic position
        actual_base = self.genome.get_sequence(chromosome, genomic_pos, genomic_pos, "+")
        actual_base = str(actual_base).upper()
        
        # Convert HGVS alleles to genomic coordinates for validation
        if strand == "-":
            complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            genomic_ref_allele = complement_map.get(ref_allele, ref_allele)
            genomic_alt_allele = complement_map.get(alt_allele, alt_allele)
        else:
            genomic_ref_allele = ref_allele
            genomic_alt_allele = alt_allele
        
        if actual_base != genomic_ref_allele.upper():
            self._debug_print(f"Reference mismatch at {genomic_pos}: expected {genomic_ref_allele}, found {actual_base}")
            return None
        
        # SIMPLE APPROACH: Find the position in the coding sequence by mapping through CDS regions
        start_codons = features[features["feature_type"] == "start_codon"]
        if start_codons.empty:
            return None
        
        start_codon_start = start_codons.iloc[0]["start"]
        cds_regions = features[features["feature_type"] == "CDS"].copy()
        
        if cds_regions.empty:
            return None
        
        # Sort CDS regions
        if strand == "+":
            cds_regions = cds_regions.sort_values('start')
        else:
            cds_regions = cds_regions.sort_values('start', ascending=False)
        
        # Find which CDS contains our mutation and the position within the coding sequence
        coding_pos = 0
        found_position = False
        
        for _, cds in cds_regions.iterrows():
            if cds["start"] <= genomic_pos <= cds["end"]:
                # Found the CDS containing our mutation
                if strand == "+":
                    # For positive strand
                    if cds["start"] >= start_codon_start:  # Only count CDS after start codon
                        offset_in_cds = genomic_pos - cds["start"]
                        final_coding_pos = coding_pos + offset_in_cds
                        found_position = True
                        break
                else:
                    # For negative strand - this is where it gets tricky
                    # The coding sequence is already in 5'->3' orientation
                    # We need to find where this genomic position maps to in the coding sequence
                    
                    # Calculate offset from the end of this CDS (since strand is -)
                    offset_from_end = cds["end"] - genomic_pos
                    final_coding_pos = coding_pos + offset_from_end
                    found_position = True
                    break
            
            # Add length of this CDS to coding position counter
            # Only count CDS regions that are part of the canonical protein
            if strand == "+":
                if cds["end"] >= start_codon_start:
                    if cds["start"] < start_codon_start:
                        # Partial CDS - only count part after start codon
                        coding_pos += cds["end"] - start_codon_start + 1
                    else:
                        # Full CDS
                        coding_pos += cds["end"] - cds["start"] + 1
            else:
                # For negative strand, we're iterating from highest to lowest genomic position
                # Add the length of each CDS
                coding_pos += cds["end"] - cds["start"] + 1
        
        if not found_position:
            self._debug_print(f"Could not map genomic position {genomic_pos} to coding sequence position")
            return None
        
        self._debug_print(f"Mapped genomic position {genomic_pos} to coding position {final_coding_pos}")
        
        # Apply the mutation to the coding sequence
        if final_coding_pos >= len(original_coding_sequence):
            self._debug_print(f"Coding position {final_coding_pos} is beyond sequence length {len(original_coding_sequence)}")
            return None
        
        # For the actual mutation in the coding sequence, use the HGVS alleles directly
        # (since the coding sequence is already in the correct orientation)
        current_base = original_coding_sequence[final_coding_pos]
        
        # Verify the reference matches what we expect in the coding sequence
        if current_base.upper() != ref_allele.upper():
            self._debug_print(f"Reference mismatch in coding sequence at pos {final_coding_pos}: expected {ref_allele}, found {current_base}")
            return None
        
        # Apply the mutation
        mutated_coding_sequence = (
            original_coding_sequence[:final_coding_pos] + 
            alt_allele.upper() + 
            original_coding_sequence[final_coding_pos + 1:]
        )
        
        self._debug_print(f"Applied mutation {ref_allele}>{alt_allele} at coding position {final_coding_pos}")
        self._debug_print(f"Generated mutated coding sequence length: {len(mutated_coding_sequence)}")
        
        return mutated_coding_sequence
    
    async def extract_gene_proteins_with_validated_mutations(self, gene_name: str, 
                                                           preferred_transcripts: Optional[Set[str]] = None,
                                                           include_mutations: bool = True,
                                                           impact_types: List[str] = None) -> Optional[List[Dict]]:
        """Extract proteins for all transcript-truncation pairs with validated mutations.
        
        This method follows the same pattern as extract_gene_proteins but adds
        validated mutation integration to canonical sequences.
        
        Args:
            gene_name: Name of the gene
            preferred_transcripts: Set of preferred transcript IDs
            include_mutations: Whether to include mutation variants
            impact_types: Mutation impact types to include
            
        Returns:
            List of enhanced protein pair dictionaries
        """
        self._debug_print(f"Processing gene {gene_name} with validated mutations")
        
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
                'truncated_mutations': []  # This will now contain canonical sequences with mutations
            }
            
            if include_mutations:
                self._debug_print(f"Processing mutations for {pair['transcript_id']} × {pair['truncation_id']}")
                
                # Get mutations in this truncation region
                truncation_feature = pair['truncation']
                mutations = await self._get_mutations_in_region(
                    gene_name, 
                    truncation_feature['start'], 
                    truncation_feature['end'],
                    impact_types
                )
                
                if not mutations.empty:
                    self._debug_print(f"Found {len(mutations)} mutations in truncation region")
                    
                    # Generate validated mutation variants applied to CANONICAL sequence
                    for _, mutation in mutations.iterrows():
                        try:
                            mutated_result = self.extract_validated_mutated_protein(
                                pair['transcript_id'], 
                                mutation,
                                apply_to_canonical=True  # Apply to canonical, not truncated
                            )
                            if mutated_result:
                                enhanced_pair['truncated_mutations'].append(mutated_result)
                                self._debug_print(f"Successfully created validated mutation variant")
                        except Exception as e:
                            self._debug_print(f"Error creating mutation variant: {e}")
                            continue
                else:
                    self._debug_print(f"No mutations found in truncation region")
                    
            enhanced_pairs.append(enhanced_pair)
            
        return enhanced_pairs