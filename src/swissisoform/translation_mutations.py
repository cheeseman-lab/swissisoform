"""Enhanced TruncatedProteinGenerator with mutation integration capability.

This extends the existing TruncatedProteinGenerator to incorporate mutations
found in truncation regions into the generated protein sequences.
"""

import asyncio
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd
import re
from Bio.Seq import Seq
from swissisoform.translation import TruncatedProteinGenerator
from swissisoform.mutations import MutationHandler


class MutationIntegratedProteinGenerator(TruncatedProteinGenerator):
    """Enhanced protein generator that can incorporate mutations into sequences."""
    
    def __init__(self, genome_handler, alt_isoform_handler, output_dir, mutation_handler=None):
        super().__init__(genome_handler, alt_isoform_handler, output_dir)
        self.mutation_handler = mutation_handler or MutationHandler()
        
    def _apply_mutation_to_sequence(self, sequence: str, position: int, 
                                  ref_allele: str, alt_allele: str, 
                                  sequence_start: int) -> str:
        """Apply a single mutation to a DNA sequence.
        
        Args:
            sequence: DNA sequence to modify
            position: Genomic position of mutation
            ref_allele: Reference allele
            alt_allele: Alternate allele  
            sequence_start: Genomic start position of the sequence
            
        Returns:
            Modified DNA sequence
        """
        # Convert genomic position to sequence position
        seq_position = position - sequence_start
        
        # Validate position is within sequence
        if seq_position < 0 or seq_position >= len(sequence):
            return sequence
            
        # Handle different mutation types
        if len(ref_allele) == 1 and len(alt_allele) == 1:
            # SNV
            if sequence[seq_position].upper() == ref_allele.upper():
                sequence = sequence[:seq_position] + alt_allele.upper() + sequence[seq_position + 1:]
        elif len(ref_allele) > len(alt_allele):
            # Deletion
            del_length = len(ref_allele) - len(alt_allele)
            if sequence[seq_position:seq_position + del_length + 1].upper() == ref_allele.upper():
                sequence = sequence[:seq_position] + alt_allele.upper() + sequence[seq_position + del_length + 1:]
        elif len(alt_allele) > len(ref_allele):
            # Insertion
            if sequence[seq_position].upper() == ref_allele.upper():
                sequence = sequence[:seq_position] + alt_allele.upper() + sequence[seq_position + 1:]
                
        return sequence
    
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
        # Create a fake alt_features DataFrame for the region
        region_features = pd.DataFrame([{
            'start': start,
            'end': end,
            'chromosome': 'chr1',  # Will be corrected by mutation handler
            'name': f'region_{start}_{end}'
        }])
        
        # Get mutations in this region (now properly awaited)
        mutations = await self.mutation_handler.get_visualization_ready_mutations(
            gene_name=gene_name,
            alt_features=region_features,
            sources=['clinvar']  # Can be parameterized
        )
        
        if mutations is None or mutations.empty:
            return pd.DataFrame()
            
        # Filter by impact types if specified
        if impact_types:
            mutations = mutations[mutations['impact'].isin(impact_types)]
            
        return mutations
    
    async def extract_truncated_protein_with_mutations(self, transcript_id: str, 
                                                      truncation_feature: pd.Series,
                                                      include_mutations: bool = True,
                                                      impact_types: List[str] = None) -> Dict:
        """Extract truncated protein with optional mutation integration.
        
        Args:
            transcript_id: Transcript ID
            truncation_feature: Truncation feature data
            include_mutations: Whether to include mutations
            impact_types: Types of mutations to include
            
        Returns:
            Dict with canonical and mutated protein sequences
        """
        # Get base truncated protein (without mutations)
        base_result = self.extract_truncated_protein(transcript_id, truncation_feature)
        if not base_result:
            return None
            
        result = {
            'transcript_id': transcript_id,
            'truncation_start': truncation_feature['start'],
            'truncation_end': truncation_feature['end'],
            'canonical_truncated': base_result,
            'mutated_variants': []
        }
        
        if not include_mutations:
            return result
            
        # Get gene name from transcript
        features = self.genome.get_transcript_features(transcript_id)
        if features.empty:
            return result
            
        gene_name = features.iloc[0].get('gene_name', '')
        if not gene_name:
            return result
            
        # Get mutations in the truncation region (now properly awaited)
        mutations = await self._get_mutations_in_region(
            gene_name, 
            truncation_feature['start'], 
            truncation_feature['end'],
            impact_types
        )
        
        if mutations.empty:
            return result
            
        # Generate mutated variants
        for _, mutation in mutations.iterrows():
            try:
                mutated_result = self._create_mutated_protein_variant(
                    base_result, mutation, transcript_id, truncation_feature
                )
                if mutated_result:
                    result['mutated_variants'].append(mutated_result)
            except Exception as e:
                print(f"Error creating mutated variant: {e}")
                continue
                
        return result
    
    def _create_mutated_protein_variant(self, base_result: Dict, mutation: pd.Series,
                                      transcript_id: str, truncation_feature: pd.Series) -> Dict:
        """Create a protein variant with a specific mutation applied.
        
        Args:
            base_result: Base truncated protein result
            mutation: Mutation data
            transcript_id: Transcript ID
            truncation_feature: Truncation feature
            
        Returns:
            Dict with mutated protein information
        """
        # Get the coding sequence regions for this transcript
        features = self.genome.get_transcript_features(transcript_id)
        cds_features = features[features['feature_type'] == 'CDS']
        
        if cds_features.empty:
            return None
            
        # Reconstruct the coding sequence with mutation
        strand = base_result['strand']
        chromosome = features.iloc[0]['chromosome']
        
        # Get mutation details
        mut_position = mutation['position']
        ref_allele = str(mutation.get('reference', ''))
        alt_allele = str(mutation.get('alternate', ''))
        
        if not ref_allele or not alt_allele or ref_allele == '.' or alt_allele == '.':
            return None
            
        # Rebuild coding sequence with mutation
        coding_sequence = ""
        mutation_applied = False
        
        # Sort CDS regions based on strand
        if strand == "+":
            sorted_cds = cds_features.sort_values('start')
        else:
            sorted_cds = cds_features.sort_values('start', ascending=False)
            
        for _, cds in sorted_cds.iterrows():
            # Get sequence for this CDS region
            cds_seq = self.genome.get_sequence(
                chromosome, cds['start'], cds['end'], strand
            )
            
            # Check if mutation falls within this CDS
            if cds['start'] <= mut_position <= cds['end'] and not mutation_applied:
                # Apply mutation to this CDS sequence
                cds_seq = self._apply_mutation_to_sequence(
                    str(cds_seq), mut_position, ref_allele, alt_allele, cds['start']
                )
                mutation_applied = True
                
            coding_sequence += str(cds_seq)
            
        # Apply truncation logic (start from after truncation)
        if strand == "+":
            alt_start_pos = truncation_feature['end'] + 1
        else:
            alt_start_pos = truncation_feature['start']
            
        # Find where to start in the coding sequence
        # This is simplified - you might need more sophisticated coordinate mapping
        truncated_coding_seq = coding_sequence  # Placeholder for truncation logic
        
        # Translate to protein
        if len(truncated_coding_seq) >= 3:
            remainder = len(truncated_coding_seq) % 3
            if remainder > 0:
                truncated_coding_seq = truncated_coding_seq[:-remainder]
            protein = str(Seq(truncated_coding_seq).translate())
        else:
            protein = ""
            
        return {
            'coding_sequence': truncated_coding_seq,
            'protein': protein,
            'mutation': {
                'position': mut_position,
                'reference': ref_allele,
                'alternate': alt_allele,
                'impact': mutation.get('impact', ''),
                'variant_id': mutation.get('variant_id', ''),
                'source': mutation.get('source', '')
            },
            'mutation_applied': mutation_applied
        }
    
    async def extract_gene_proteins_with_mutations(self, gene_name: str, 
                                                  preferred_transcripts: Optional[Set[str]] = None,
                                                  include_mutations: bool = True,
                                                  impact_types: List[str] = None) -> List[Dict]:
        """Extract proteins for all transcript-truncation pairs with mutation variants.
        
        Args:
            gene_name: Gene name
            preferred_transcripts: Preferred transcript IDs
            include_mutations: Whether to generate mutation variants
            impact_types: Mutation impact types to include
            
        Returns:
            List of enhanced protein pair dictionaries
        """
        # Get base pairs (canonical + truncated without mutations)
        base_pairs = self.extract_gene_proteins(gene_name, preferred_transcripts)
        if not base_pairs:
            return []
            
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
                # Get mutations for this truncation region (now properly awaited)
                truncation_feature = pair['truncation']
                mutated_result = await self.extract_truncated_protein_with_mutations(
                    pair['transcript_id'], 
                    truncation_feature,
                    include_mutations=True,
                    impact_types=impact_types
                )
                
                if mutated_result and mutated_result.get('mutated_variants'):
                    enhanced_pair['truncated_mutations'] = mutated_result['mutated_variants']
                    
            enhanced_pairs.append(enhanced_pair)
            
        return enhanced_pairs
        
    async def create_mutation_integrated_dataset(self, gene_list: List[str],
                                                preferred_transcripts: Optional[Set[str]] = None,
                                                include_mutations: bool = True,
                                                impact_types: List[str] = None,
                                                min_length: int = 10,
                                                max_length: int = 100000) -> pd.DataFrame:
        """Create dataset with canonical, truncated, and mutated protein sequences.
        
        Args:
            gene_list: List of gene names
            preferred_transcripts: Preferred transcript IDs
            include_mutations: Whether to include mutation variants
            impact_types: Mutation impact types to include
            min_length: Minimum protein length
            max_length: Maximum protein length
            
        Returns:
            Enhanced dataset DataFrame
        """
        all_sequences = []
        
        for gene_name in gene_list:
            try:
                enhanced_pairs = await self.extract_gene_proteins_with_mutations(
                    gene_name, preferred_transcripts, include_mutations, impact_types
                )
                
                for pair in enhanced_pairs:
                    transcript_id = pair['transcript_id']
                    truncation_id = pair['truncation_id']
                    
                    # Add canonical sequence
                    canonical_protein = pair['canonical']['protein']
                    if min_length <= len(canonical_protein) <= max_length:
                        all_sequences.append({
                            'gene': gene_name,
                            'transcript_id': transcript_id,
                            'truncation_id': truncation_id,
                            'variant_id': 'canonical',
                            'sequence': canonical_protein,
                            'length': len(canonical_protein),
                            'variant_type': 'canonical',
                            'mutation_info': None
                        })
                    
                    # Add base truncated sequence
                    truncated_protein = pair['truncated_base']['protein']
                    if min_length <= len(truncated_protein) <= max_length:
                        all_sequences.append({
                            'gene': gene_name,
                            'transcript_id': transcript_id,
                            'truncation_id': truncation_id,
                            'variant_id': f'{truncation_id}_base',
                            'sequence': truncated_protein,
                            'length': len(truncated_protein),
                            'variant_type': 'truncated',
                            'mutation_info': None
                        })
                    
                    # Add mutated truncated sequences
                    for i, mutated_variant in enumerate(pair['truncated_mutations']):
                        mutated_protein = mutated_variant['protein']
                        if min_length <= len(mutated_protein) <= max_length:
                            mutation_info = mutated_variant['mutation']
                            variant_id = f"{truncation_id}_mut_{mutation_info['position']}_{mutation_info['alternate']}"
                            
                            all_sequences.append({
                                'gene': gene_name,
                                'transcript_id': transcript_id,
                                'truncation_id': truncation_id,
                                'variant_id': variant_id,
                                'sequence': mutated_protein,
                                'length': len(mutated_protein),
                                'variant_type': 'truncated_mutated',
                                'mutation_info': mutation_info
                            })
                            
            except Exception as e:
                print(f"Error processing gene {gene_name}: {e}")
                continue
                
        return pd.DataFrame(all_sequences)