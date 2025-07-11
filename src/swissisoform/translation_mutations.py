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
            
        self._debug_print(f"Unsupported HGVS format: {hgvsc}")
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
                                        truncation_feature: pd.Series,
                                        mutation: pd.Series) -> Optional[Dict]:
        """Extract truncated protein with validated mutation applied.
        
        This method follows the validated approach:
        1. Validate reference allele at genomic position
        2. Apply mutation at exact genomic position
        3. Validate protein change against expected HGVSP
        
        Args:
            transcript_id: Transcript ID to process
            truncation_feature: Truncation feature data
            mutation: Mutation data with position, HGVS notation
            
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
        
        # For negative strand genes, convert alternate allele too
        if strand == "-":
            complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            genomic_alt_allele = complement_map.get(alt_allele, alt_allele)
        else:
            genomic_alt_allele = alt_allele
        
        # Step 3: Generate original truncated protein (for comparison)
        original_result = self.extract_truncated_protein(transcript_id, truncation_feature)
        if not original_result:
            self._debug_print(f"Could not extract original truncated protein")
            return None
            
        original_protein = original_result['protein']
        
        # Step 4: Generate mutated truncated protein
        mutated_result = self._apply_validated_genomic_mutation(
            transcript_id, truncation_feature, genomic_pos, 
            expected_genomic_base, alt_allele, strand
        )
        if not mutated_result:
            self._debug_print(f"Could not apply mutation")
            return None
            
        mutated_protein = mutated_result['protein']
        
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
            'coding_sequence': mutated_result['coding_sequence'],
            'protein': mutated_protein,
            'strand': strand,
            'transcript_id': transcript_id,
            'mutation': {
                'genomic_position': genomic_pos,
                'reference': ref_allele,  # Keep original HGVS reference
                'alternate': alt_allele,  # Keep original HGVS alternate
                'genomic_reference': expected_genomic_base,  # Add genomic reference
                'genomic_alternate': genomic_alt_allele if strand == "-" else alt_allele,  # Add genomic alternate
                'hgvsc': hgvsc,
                'hgvsp': hgvsp,
                'impact': mutation.get('impact', ''),
                'variant_id': mutation.get('variant_id', ''),
                'source': mutation.get('source', '')
            },
            'validation': validation_result
        }
    
    def _apply_validated_genomic_mutation(self, transcript_id: str, 
                                        truncation_feature: pd.Series,
                                        genomic_pos: int, genomic_ref_allele: str, 
                                        coding_alt_allele: str, strand: str) -> Optional[Dict]:
        """Apply mutation at specific genomic position following truncation logic.
        
        This method applies the same truncation logic as extract_truncated_protein
        but with the mutation applied at the exact genomic position.
        
        Args:
            transcript_id: Transcript ID
            truncation_feature: Truncation feature data
            genomic_pos: Genomic position of mutation
            genomic_ref_allele: Reference allele in genomic coordinates (+ strand)
            coding_alt_allele: Alternate allele from HGVS (coding strand)
            strand: Gene strand
            
        Returns:
            Dict with mutated protein data or None if failed
        """
        # Convert coding strand alternate allele to genomic coordinates
        if strand == "-":
            complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            genomic_alt_allele = complement_map.get(coding_alt_allele, coding_alt_allele)
            self._debug_print(f"Negative strand: coding alt={coding_alt_allele} -> genomic alt={genomic_alt_allele}")
        else:
            genomic_alt_allele = coding_alt_allele
        # Get features and transcript data (same as extract_truncated_protein)
        features = self.genome.get_transcript_features(transcript_id)
        transcript_data = self.genome.get_transcript_features_with_sequence(transcript_id)

        if not transcript_data or "sequence" not in transcript_data:
            return None

        transcript_start = transcript_data["sequence"]["start"]
        transcript_end = transcript_data["sequence"]["end"]
        strand = transcript_data["sequence"]["strand"]
        chromosome = transcript_data["sequence"]["chromosome"]

        # Get truncation positions (same logic as extract_truncated_protein)
        truncation_start = truncation_feature["start"]
        truncation_end = truncation_feature["end"]

        # Find the original stop codon
        stop_codons = features[features["feature_type"] == "stop_codon"]
        stop_codon_start = None
        stop_codon_end = None
        if not stop_codons.empty:
            stop_codon_start = stop_codons.iloc[0]["start"]
            stop_codon_end = stop_codons.iloc[0]["end"]

        # Apply truncation logic (same as extract_truncated_protein)
        if strand == "+":
            alt_start_pos = truncation_end + 1
            extract_end = stop_codon_end if stop_codon_end else None
        else:
            alt_start_pos = truncation_start
            extract_end = stop_codon_start if stop_codon_start else None

        # Get CDS regions that overlap with our new range (same logic)
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

            if effective_end >= effective_start:
                overlapping_cds.append(
                    {
                        "start": effective_start,
                        "end": effective_end,
                        "length": effective_end - effective_start + 1,
                    }
                )

        # Sort CDS regions properly for extraction (same logic)
        if strand == "+":
            overlapping_cds.sort(key=lambda x: x["start"])
        else:
            overlapping_cds.sort(key=lambda x: x["start"], reverse=True)

        # Extract sequence from each CDS region WITH MUTATION APPLIED
        coding_sequence = ""
        mutation_applied = False
        
        for cds in overlapping_cds:
            # Get sequence for this CDS region
            cds_seq = self.genome.get_sequence(
                chromosome, cds["start"], cds["end"], strand
            )
            cds_seq = str(cds_seq)
            
            # Check if mutation falls within this CDS region
            if cds["start"] <= genomic_pos <= cds["end"]:
                # Apply mutation at exact genomic position
                seq_pos = genomic_pos - cds["start"]
                if seq_pos < len(cds_seq) and cds_seq[seq_pos].upper() == genomic_ref_allele.upper():
                    cds_seq = cds_seq[:seq_pos] + genomic_alt_allele.upper() + cds_seq[seq_pos + 1:]
                    mutation_applied = True
                    self._debug_print(f"Applied {genomic_ref_allele}>{genomic_alt_allele} at genomic pos {genomic_pos} (CDS pos {seq_pos})")
                else:
                    current_base = cds_seq[seq_pos] if seq_pos < len(cds_seq) else 'OOB'
                    self._debug_print(f"Reference mismatch in CDS at pos {seq_pos}: expected {genomic_ref_allele}, found {current_base}")
                    # Continue anyway since we already validated at the genomic level
            
            coding_sequence += cds_seq

        if not mutation_applied:
            self._debug_print(f"Warning: mutation not applied to any CDS region")
            # Continue anyway to see if we can still generate a valid protein

        # Translate the sequence (same logic as extract_truncated_protein)
        if len(coding_sequence) >= 3:
            remainder = len(coding_sequence) % 3
            if remainder > 0:
                coding_sequence = coding_sequence[:-remainder]
            protein = str(Seq(coding_sequence).translate())
        else:
            protein = ""

        self._debug_print(f"Generated mutated protein length: {len(protein)}")

        return {
            "coding_sequence": coding_sequence,
            "protein": protein,
            "strand": strand,
            "transcript_id": transcript_id,
            "truncation_start": truncation_start,
            "truncation_end": truncation_end,
            "alternative_start_pos": alt_start_pos,
            "total_cds_regions": len(overlapping_cds),
            "mutation_applied": mutation_applied
        }
    
    async def extract_gene_proteins_with_validated_mutations(self, gene_name: str, 
                                                           preferred_transcripts: Optional[Set[str]] = None,
                                                           include_mutations: bool = True,
                                                           impact_types: List[str] = None) -> Optional[List[Dict]]:
        """Extract proteins for all transcript-truncation pairs with validated mutations.
        
        This method follows the same pattern as extract_gene_proteins but adds
        validated mutation integration.
        
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
                'truncated_mutations': []
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
                    
                    # Generate validated mutation variants
                    for _, mutation in mutations.iterrows():
                        try:
                            mutated_result = self.extract_validated_mutated_protein(
                                pair['transcript_id'], 
                                truncation_feature,
                                mutation
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