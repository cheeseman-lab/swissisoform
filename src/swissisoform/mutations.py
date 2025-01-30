from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport
import pandas as pd
from typing import Dict, Optional
import requests
import xml.etree.ElementTree as ET
from urllib.parse import urlencode
import json
import asyncio
from datetime import datetime

class MutationHandler:
    """
    A class to handle variant data from multiple mutation databases.
    Provides methods to fetch, parse, and analyze mutation data from:
    - gnomAD
    - ClinVar (to be implemented with correct API)
    - Other databases can be added...
    """
    
    def __init__(self):
        """Initialize the mutation handler with API endpoints"""
        # Setup gnomAD GraphQL client
        transport = AIOHTTPTransport(url="https://gnomad.broadinstitute.org/api")
        self.gnomad_client = Client(transport=transport, fetch_schema_from_transport=True)
        self.cached_data = {}
        
        # ClinVar E-utilities base URL and endpoints
        self.clinvar_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.clinvar_endpoints = {
            'esearch': f"{self.clinvar_base}/esearch.fcgi",
            'esummary': f"{self.clinvar_base}/esummary.fcgi",
            'efetch': f"{self.clinvar_base}/efetch.fcgi",
            'elink': f"{self.clinvar_base}/elink.fcgi"
        }
        
    async def get_gnomad_variants(self, gene_name: str) -> pd.DataFrame:
        """Get processed variant data from gnomAD"""
        gnomad_data = await self.fetch_gnomad_data(gene_name)
        return self.process_gnomad_variants(gnomad_data)
    
    async def get_gnomad_summary(self, gene_name: str) -> Dict:
        """Get summary statistics for gnomAD variants"""
        gnomad_df = await self.get_gnomad_variants(gene_name)
        
        if gnomad_df.empty:
            return {
                'total_variants': 0,
                'consequence_types': {},
                'mean_allele_frequency': 0.0,
                'variants_by_impact': {
                    'with_hgvsc': 0,
                    'with_hgvsp': 0
                }
            }
        
        summary = {
            'total_variants': len(gnomad_df),
            'consequence_types': gnomad_df['consequence'].value_counts().to_dict(),
            'mean_allele_frequency': gnomad_df['allele_frequency'].mean(),
            'variants_by_impact': {
                'with_hgvsc': gnomad_df['hgvsc'].notna().sum(),
                'with_hgvsp': gnomad_df['hgvsp'].notna().sum()
            }
        }
        return summary
        
    async def fetch_gnomad_data(self, gene_name: str, reference_genome: str = "GRCh38") -> Dict:
        """
        Fetch variant data from gnomAD API for a specific gene using GraphQL.
        
        Args:
            gene_name (str): Gene symbol (e.g., 'BRCA1')
            reference_genome (str): Reference genome version (default: 'GRCh38')
            
        Returns:
            Dict: Processed gnomAD data for the gene
        """
        cache_key = f"gnomad_{gene_name}_{reference_genome}"
        if cache_key in self.cached_data:
            return self.cached_data[cache_key]
            
        # Define the GraphQL query
        query = gql("""
            query VariantsInGene($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
              gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
                gene_id
                symbol
                variants(dataset: gnomad_r4) {
                  variant_id
                  pos
                  ref
                  alt
                  consequence
                  transcript_id
                  transcript_version
                  hgvs
                  hgvsc
                  hgvsp
                  flags
                  exome {
                    ac
                    ac_hemi
                    ac_hom
                    an
                    af
                    filters
                    populations {
                      id
                      ac
                      an
                      ac_hom
                      ac_hemi
                    }
                  }
                  genome {
                    ac
                    ac_hemi
                    ac_hom
                    an
                    af
                    filters
                    populations {
                      id
                      ac
                      an
                      ac_hom
                      ac_hemi
                    }
                  }
                }
              }
            }
        """)
        
        # Set up variables for the query
        variables = {
            "geneSymbol": gene_name,
            "referenceGenome": reference_genome
        }
        
        try:
            # Execute the query
            result = await self.gnomad_client.execute_async(query, variable_values=variables)
            
            # Cache the response
            if result.get('gene') is not None:
                self.cached_data[cache_key] = result
            return result
            
        except Exception as e:
            raise ConnectionError(f"Failed to fetch gnomAD data: {str(e)}")
    
    def process_gnomad_variants(self, data: Dict, filter_field: Optional[str] = None, filter_value: Optional[str] = None) -> pd.DataFrame:
        """
        Process gnomAD variant data into a pandas DataFrame.
        
        Args:
            data (Dict): Raw gnomAD API response
            filter_field (str, optional): Field to filter on (e.g., 'consequence')
            filter_value (str, optional): Value to filter for (e.g., '5_prime_UTR_variant')
            
        Returns:
            pd.DataFrame: Processed variant data
        """
        if not data or 'gene' not in data or not data['gene'] or 'variants' not in data['gene']:
            return pd.DataFrame()
            
        variants = data['gene']['variants']
        processed_variants = []
        
        for variant in variants:
            # Skip this variant if it doesn't match the filter condition
            if filter_field and filter_value:
                variant_field_value = variant.get(filter_field)
                if variant_field_value != filter_value:
                    continue
            
            # Get frequency data from exome or genome
            exome_data = variant.get('exome', {})
            genome_data = variant.get('genome', {})
            
            # Prefer exome data if available, otherwise use genome data
            freq_data = exome_data or genome_data
            
            if not freq_data:  # Skip if no frequency data available
                continue
            
            variant_info = {
                'position': variant['pos'],
                'variant_id': variant['variant_id'],
                'reference': variant['ref'],
                'alternate': variant['alt'],
                'consequence': variant['consequence'],
                'transcript_id': variant.get('transcript_id'),
                'transcript_version': variant.get('transcript_version'),
                'hgvs': variant.get('hgvs'),
                'hgvsc': variant.get('hgvsc'),
                'hgvsp': variant.get('hgvsp'),
                'flags': ','.join(variant.get('flags', [])),
                
                'allele_count': freq_data.get('ac'),
                'allele_count_hom': freq_data.get('ac_hom'),
                'allele_count_hemi': freq_data.get('ac_hemi'),
                'allele_number': freq_data.get('an'),
                'allele_frequency': freq_data.get('af'),
                'filters': ','.join(freq_data.get('filters', [])),
                
                'freq_source': 'exome' if exome_data else 'genome'
            }
            
            # Add population-specific data if available
            if 'populations' in freq_data:
                for pop in freq_data['populations']:
                    pop_id = pop['id']
                    variant_info[f'{pop_id}_ac'] = pop.get('ac')
                    variant_info[f'{pop_id}_an'] = pop.get('an')
                    variant_info[f'{pop_id}_ac_hom'] = pop.get('ac_hom')
                    variant_info[f'{pop_id}_ac_hemi'] = pop.get('ac_hemi')
            
            processed_variants.append(variant_info)
            
        return pd.DataFrame(processed_variants)
    
    async def get_clinvar_variants(self, gene_name: str) -> pd.DataFrame:
        """
        Get all ClinVar variants for a gene.
        
        Args:
            gene_name (str): Gene symbol (e.g., 'BRCA1')
            
        Returns:
            pd.DataFrame: Processed ClinVar variant data
        """
        cache_key = f"clinvar_{gene_name}"
        if cache_key in self.cached_data:
            return self.cached_data[cache_key]

        # Get variant IDs from esearch
        search_params = {
            'db': 'clinvar',
            'term': f"{gene_name}[gene] AND single_gene[prop]",
            'retmax': 500,  # Adjust as needed
            'retmode': 'json',
            'tool': 'swissisoform',
            'email': 'your.email@example.com'  # Replace with actual email
        }
        
        try:
            # First get the variant IDs
            response = await self._async_get_request(
                self.clinvar_endpoints['esearch'],
                params=search_params
            )
            await asyncio.sleep(0.35)  # Rate limit: max 3 requests per second
            
            data = json.loads(response)
            if 'esearchresult' not in data or 'idlist' not in data['esearchresult']:
                return pd.DataFrame()
                
            variant_ids = data['esearchresult']['idlist']
            if not variant_ids:
                return pd.DataFrame()
            
            # Process variants in batches
            batch_size = 50
            all_variants = []
            
            for i in range(0, len(variant_ids), batch_size):
                batch_ids = variant_ids[i:i + batch_size]
                
                # Get summaries for batch
                summary_params = {
                    'db': 'clinvar',
                    'id': ','.join(batch_ids),
                    'retmode': 'json',
                    'tool': 'swissisoform',
                    'email': 'your.email@example.com'  # Replace with actual email
                }
                
                summary_response = await self._async_get_request(
                    self.clinvar_endpoints['esummary'],
                    params=summary_params
                )
                await asyncio.sleep(0.35)  # Rate limit: max 3 requests per second
                
                summary_data = json.loads(summary_response)
                batch_variants = self._process_clinvar_variants(summary_data)
                
                if not batch_variants.empty:
                    all_variants.append(batch_variants)
            
            if not all_variants:
                return pd.DataFrame()
                
            final_df = pd.concat(all_variants, ignore_index=True)
            self.cached_data[cache_key] = final_df
            return final_df
            
        except Exception as e:
            print(f"Error fetching ClinVar data: {str(e)}")
            return pd.DataFrame()

    def _process_clinvar_variants(self, data: Dict) -> pd.DataFrame:
        """
        Process ClinVar variant data into a pandas DataFrame.
        
        Args:
            data (Dict): Raw ClinVar API response
            
        Returns:
            pd.DataFrame: Processed variant data with key fields
        """
        if 'result' not in data or 'uids' not in data['result']:
            return pd.DataFrame()
            
        variants = []
        for uid in data['result']['uids']:
            variant_data = data['result'][uid]
            
            # Extract basic variant information
            variant_info = {
                'variant_id': uid,
                'obj_type': variant_data.get('obj_type', ''),
                'accession': variant_data.get('accession', ''),
                'title': variant_data.get('title', ''),
                'gene_name': '',  # Will be populated from genes list
                'review_status': '',  # Will be populated from classifications
                'clinical_significance': '',  # Will be populated from classifications
                'last_evaluated': '',  # Will be populated from classifications
                'molecular_consequences': '',
                'protein_change': variant_data.get('protein_change', ''),
                'chromosome': '',  # Will be populated from variation_set
                'start': '',  # Will be populated from variation_set
                'stop': '',  # Will be populated from variation_set
                'ref_allele': '',  # Will be populated from variation_set
                'alt_allele': '',  # Will be populated from variation_set
                'submission_count': len(variant_data.get('supporting_submissions', {}).get('scv', []))
            }
            
            # Get gene information
            genes = variant_data.get('genes', [])
            if genes:
                variant_info['gene_name'] = genes[0].get('symbol', '')
                variant_info['gene_id'] = genes[0].get('geneid', '')
                variant_info['gene_strand'] = genes[0].get('strand', '')
                
            # Get clinical significance information
            # Check all possible classification types
            for class_type in ['germline_classification', 'clinical_impact_classification', 'oncogenicity_classification']:
                if class_type in variant_data:
                    classification = variant_data[class_type]
                    if classification.get('description'):  # Only use if there's a description
                        variant_info['clinical_significance'] = classification.get('description', '')
                        variant_info['review_status'] = classification.get('review_status', '')
                        variant_info['last_evaluated'] = classification.get('last_evaluated', '')
                        break  # Use the first classification that has a description
            
            # Get molecular consequences
            mol_cons = variant_data.get('molecular_consequence_list', [])
            variant_info['molecular_consequences'] = ';'.join(mol_cons) if mol_cons else ''
            
            # Get variant location information from variation_set
            variation_sets = variant_data.get('variation_set', [])
            if variation_sets:
                var_set = variation_sets[0]  # Use first variation set
                
                # Get variant details
                variant_info['variant_name'] = var_set.get('variation_name', '')
                variant_info['cdna_change'] = var_set.get('cdna_change', '')
                variant_info['canonical_spdi'] = var_set.get('canonical_spdi', '')
                
                # Get location information
                locations = var_set.get('variation_loc', [])
                for loc in locations:
                    if loc.get('assembly_name') == 'GRCh38' and loc.get('status') == 'current':
                        variant_info['chromosome'] = loc.get('chr', '')
                        variant_info['start'] = loc.get('start', '')
                        variant_info['stop'] = loc.get('stop', '')
                        variant_info['cytogenic_location'] = loc.get('band', '')
                        break
            
            variants.append(variant_info)
        
        df = pd.DataFrame(variants)
        
        # Convert relevant columns to appropriate types
        numeric_cols = ['start', 'stop', 'submission_count']
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors='ignore')
        
        return df

    async def get_clinvar_summary(self, gene_name: str) -> Dict:
        """
        Get comprehensive summary statistics for ClinVar variants.
        
        Args:
            gene_name (str): Gene symbol
            
        Returns:
            Dict: Summary of variant information including:
                - Overall counts
                - Clinical significance distribution
                - Variant type distribution
                - Molecular consequences
                - Location statistics
                - Review status metrics
        """
        variants_df = await self.get_clinvar_variants(gene_name)
        
        if variants_df.empty:
            return {
                'total_variants': 0,
                'clinical_significance': {},
                'variant_types': {},
                'molecular_consequences': {},
                'chromosome_distribution': {},
                'review_status': {},
                'metadata': {
                    'submission_count': {'total': 0, 'mean': 0, 'max': 0},
                    'unique_accessions': 0
                }
            }

        summary = {
            'total_variants': len(variants_df),
            
            # Clinical significance distribution
            'clinical_significance': variants_df['clinical_significance']
                .value_counts()
                .to_dict(),
                
            # Variant types
            'variant_types': variants_df['obj_type']
                .value_counts()
                .to_dict(),
                
            # Molecular consequences (split by semicolon and count unique)
            'molecular_consequences': variants_df['molecular_consequences']
                .str.split(';')
                .explode()
                .value_counts()
                .to_dict(),
                
            # Chromosome distribution
            'chromosome_distribution': variants_df['chromosome']
                .value_counts()
                .to_dict(),
                
            # Review status metrics
            'review_status': variants_df['review_status']
                .value_counts()
                .to_dict(),
                
            # Additional metadata
            'metadata': {
                'submission_count': {
                    'total': variants_df['submission_count'].sum(),
                    'mean': round(variants_df['submission_count'].mean(), 2),
                    'max': variants_df['submission_count'].max()
                },
                'unique_accessions': variants_df['accession'].nunique()
            }
        }
        
        # Add pathogenicity categorization
        pathogenicity_mapping = {
            'likely_pathogenic_or_pathogenic': variants_df['clinical_significance']
                .str.contains('pathogenic', case=False, na=False) & 
                ~variants_df['clinical_significance'].str.contains('conflicting', case=False, na=False),
            'benign_or_likely_benign': variants_df['clinical_significance']
                .str.contains('benign', case=False, na=False) & 
                ~variants_df['clinical_significance'].str.contains('conflicting', case=False, na=False),
            'uncertain_significance': variants_df['clinical_significance']
                .str.contains('uncertain|conflicting', case=False, na=False),
            'not_provided': variants_df['clinical_significance'].isna() | 
                (variants_df['clinical_significance'] == '')
        }
        
        summary['pathogenicity_categories'] = {
            category: int(variants_df[mask].shape[0])
            for category, mask in pathogenicity_mapping.items()
        }
        
        # Add variant location statistics if coordinates are present
        if variants_df['start'].notna().any():
            summary['location_stats'] = {
                'min_position': int(variants_df['start'].min()),
                'max_position': int(variants_df['stop'].max()),
                'span': int(variants_df['stop'].max() - variants_df['start'].min())
            }
        
        return summary
    
    def _process_hgvs_notation(self, notation: str) -> str:
        """
        Process HGVS notation to remove transcript/protein IDs.
        
        Args:
            notation (str): Raw HGVS notation (e.g., "ENST00000368235.8:c.8G>C")
            
        Returns:
            str: Cleaned HGVS notation (e.g., "c.8G>C")
        """
        if pd.isna(notation) or not isinstance(notation, str):
            return None
        
        # Split on colon and take the last part
        parts = notation.split(':')
        return parts[-1] if parts else None

    def _parse_variant_id(self, var_id: str) -> Dict[str, str]:
        """
        Parse variant ID in format "chr:position:ref:alt".
        
        Args:
            var_id (str): Variant ID (e.g., "1:156595364:G:A")
            
        Returns:
            Dict: Dictionary containing chromosome, position, reference, and alternate alleles
        """
        if pd.isna(var_id) or not isinstance(var_id, str):
            return {
                'chromosome': None,
                'position': None,
                'reference': None,
                'alternate': None
            }
        
        try:
            parts = var_id.split(':')
            if len(parts) == 4:
                return {
                    'chromosome': parts[0],
                    'position': int(parts[1]),
                    'reference': parts[2],
                    'alternate': parts[3]
                }
        except (ValueError, IndexError):
            pass
        
        return {
            'chromosome': None,
            'position': None,
            'reference': None,
            'alternate': None
        }

    def get_aggregator_variants(self, aggregator_csv_path: str) -> pd.DataFrame:
        """
        Load and parse aggregator variant data from a user-specified CSV.
        
        Args:
            aggregator_csv_path (str): Path to the aggregator CSV file
        Returns:
            pd.DataFrame: Parsed aggregator variant data with processed columns
        """
        cache_key = f"aggregator_{aggregator_csv_path}"
        if cache_key in self.cached_data:
            return self.cached_data[cache_key]
        
        try:
            # Read the CSV file
            df = pd.read_csv(aggregator_csv_path)
            
            # Process HGVS notations
            df['hgvsc'] = df['hgvsc'].apply(self._process_hgvs_notation)
            df['hgvsp'] = df['hgvsp'].apply(self._process_hgvs_notation)
            
            # Parse variant IDs
            var_info = df['varId'].apply(self._parse_variant_id)
            var_info_df = pd.DataFrame(var_info.tolist())
            
            # Add parsed information to the DataFrame
            df = pd.concat([df, var_info_df], axis=1)
            
            # Convert numeric columns
            numeric_cols = ['allelecount', 'allelnumber', 'allelefrequency', 
                        'homozygouscount', 'gnomADg_AC', 'gnomADg_AN', 'gnomADg_AF']
            for col in numeric_cols:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            
            self.cached_data[cache_key] = df
            return df
            
        except Exception as e:
            print(f"Error processing aggregator CSV: {str(e)}")
            return pd.DataFrame()

    def get_aggregator_summary(self, aggregator_csv_path: str) -> Dict:
        """
        Generate a summary of aggregator variant data from a user-specified CSV.
        
        Args:
            aggregator_csv_path (str): Path to the aggregator CSV file
        Returns:
            Dict: Summary of aggregator variant data including impact distribution,
                allele frequencies, and position information
        """
        aggregator_df = self.get_aggregator_variants(aggregator_csv_path)
        
        if aggregator_df.empty:
            return {
                "csv_path": aggregator_csv_path,
                "total_variants": 0,
                "impact_counts": {},
                "position_stats": {},
                "frequency_stats": {}
            }
        
        summary = {
            "csv_path": aggregator_csv_path,
            "total_variants": len(aggregator_df),
            
            # Impact and consequence distributions
            "impact_counts": aggregator_df["impact"].value_counts().to_dict(),
            "consequence_counts": aggregator_df["consequence"].value_counts().to_dict(),
            
            # Position statistics
            "position_stats": {
                "chromosomes": aggregator_df["chromosome"].unique().tolist(),
                "min_position": int(aggregator_df["position"].min()) if not pd.isna(aggregator_df["position"].min()) else None,
                "max_position": int(aggregator_df["position"].max()) if not pd.isna(aggregator_df["position"].max()) else None
            },
            
            # Frequency statistics
            "frequency_stats": {
                "mean_frequency": aggregator_df["allelefrequency"].mean(),
                "median_frequency": aggregator_df["allelefrequency"].median(),
                "max_frequency": aggregator_df["allelefrequency"].max()
            },
            
            # HGVS notation statistics
            "notation_stats": {
                "with_hgvsc": aggregator_df["hgvsc"].notna().sum(),
                "with_hgvsp": aggregator_df["hgvsp"].notna().sum()
            }
        }
        
        return summary

    async def _async_get_request(self, url: str, params: Dict = None) -> str:
        """
        Make an async GET request.
        
        Args:
            url (str): URL to request
            params (Dict, optional): Query parameters
            
        Returns:
            str: Response text
        """
        loop = asyncio.get_event_loop()
        def _make_request():
            response = requests.get(url, params=params)
            response.raise_for_status()
            return response.text
            
        return await loop.run_in_executor(None, _make_request)
        
    def clear_cache(self):
        """Clear the cached API responses"""
        self.cached_data = {}

    def standardize_mutation_data(self, variants_df: pd.DataFrame, source: str) -> pd.DataFrame:
        """
        Standardize mutation data from different sources into a consistent format.
        
        Args:
            variants_df (pd.DataFrame): Raw variants DataFrame from any source
            source (str): Source of the data ('gnomad', 'clinvar', or 'aggregator')
            
        Returns:
            pd.DataFrame: Standardized mutation data with consistent columns
        """
        # Define standard columns for visualization
        standard_columns = [
            'position',          # Genomic position
            'variant_id',        # Unique identifier
            'reference',         # Reference allele
            'alternate',         # Alternate allele
            'source',           # Data source
            'impact',           # Functional impact/consequence
            'hgvsc',            # Coding sequence change
            'hgvsp',            # Protein change
            'allele_frequency', # Frequency in population
            'clinical_significance'  # Clinical interpretation
        ]
        
        standardized_df = pd.DataFrame(columns=standard_columns)
        
        if variants_df.empty:
            return standardized_df
            
        if source.lower() == 'gnomad':
            standardized_df = pd.DataFrame({
                'position': variants_df['position'],
                'variant_id': variants_df['variant_id'],
                'reference': variants_df['reference'],
                'alternate': variants_df['alternate'],
                'source': 'gnomAD',
                'impact': variants_df['consequence'],
                'hgvsc': variants_df['hgvsc'],
                'hgvsp': variants_df['hgvsp'],
                'allele_frequency': variants_df['allele_frequency'],
                'clinical_significance': None
            })
            
        elif source.lower() == 'clinvar':
            standardized_df = pd.DataFrame({
                'position': variants_df['start'],
                'variant_id': variants_df['accession'],
                'reference': variants_df['ref_allele'],
                'alternate': variants_df['alt_allele'],
                'source': 'ClinVar',
                'impact': variants_df['molecular_consequences'],
                'hgvsc': variants_df['cdna_change'],
                'hgvsp': variants_df['protein_change'],
                'allele_frequency': None,
                'clinical_significance': variants_df['clinical_significance']
            })
            
        elif source.lower() == 'aggregator':
            standardized_df = pd.DataFrame({
                'position': variants_df['position'],  # Now comes from variant ID parsing
                'variant_id': variants_df['varId'],
                'reference': variants_df['reference'],  # Now comes from variant ID parsing
                'alternate': variants_df['alternate'],  # Now comes from variant ID parsing
                'source': 'Aggregator',
                'impact': variants_df['consequence'],  # Changed from impact to consequence to match gnomAD
                'hgvsc': variants_df['hgvsc'],
                'hgvsp': variants_df['hgvsp'],
                'allele_frequency': variants_df['allelefrequency'],
                'clinical_significance': None
            })
        
        return self._clean_mutation_data(standardized_df)

    def _parse_hgvsc(self, hgvsc: str) -> dict:
        """
        Parse position and alleles from HGVS coding sequence notation.
        
        Args:
            hgvsc (str): HGVS coding sequence notation
            
        Returns:
            dict: Dictionary containing position, reference, and alternate alleles
        """
        if not isinstance(hgvsc, str):
            return {'position': None, 'reference': None, 'alternate': None}
        
        try:
            # Extract position and change from strings like "c.123A>G" or "c.123del"
            import re
            match = re.search(r'c\.(\d+)([ACGT]?)>?([ACGT]|del|ins|dup)?', hgvsc)
            if match:
                position = int(match.group(1))
                reference = match.group(2) if match.group(2) else None
                alternate = match.group(3) if match.group(3) else None
                return {
                    'position': position,
                    'reference': reference,
                    'alternate': alternate
                }
        except:
            pass
        
        return {'position': None, 'reference': None, 'alternate': None}

    def _clean_mutation_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Clean and validate the standardized mutation data.
        
        Args:
            df (pd.DataFrame): Standardized mutation DataFrame
            
        Returns:
            pd.DataFrame: Cleaned mutation data
        """
        # Convert position to integer if possible
        df['position'] = pd.to_numeric(df['position'], errors='coerce')
        
        # Clean frequency data
        df['allele_frequency'] = pd.to_numeric(df['allele_frequency'], errors='coerce')
        
        # Fill missing values
        df['clinical_significance'] = df['clinical_significance'].fillna('Unknown')
        df['impact'] = df['impact'].fillna('Unknown')
        
        # Clean string columns
        string_columns = ['hgvsc', 'hgvsp', 'reference', 'alternate']
        for col in string_columns:
            df[col] = df[col].fillna('.')
            df[col] = df[col].astype(str)
        
        # Sort by position if available
        if df['position'].notna().any():
            df = df.sort_values('position')
        
        return df.reset_index(drop=True)        
    
    async def get_visualization_ready_mutations(self, gene_name: str, 
                                            aggregator_csv_path: Optional[str] = None) -> pd.DataFrame:
        """
        Get mutation data from all sources in a format ready for visualization.
        
        Args:
            gene_name (str): Gene symbol
            aggregator_csv_path (str, optional): Path to aggregator CSV
            
        Returns:
            pd.DataFrame: Combined and standardized mutation data
        """
        # Get data from each source
        try:
            gnomad_df = await self.get_gnomad_variants(gene_name)
            gnomad_standardized = self.standardize_mutation_data(gnomad_df, 'gnomad')
        except Exception as e:
            print(f"Error fetching gnomAD data: {str(e)}")
            gnomad_standardized = pd.DataFrame()
        
        try:
            clinvar_df = await self.get_clinvar_variants(gene_name)
            clinvar_standardized = self.standardize_mutation_data(clinvar_df, 'clinvar')
        except Exception as e:
            print(f"Error fetching ClinVar data: {str(e)}")
            clinvar_standardized = pd.DataFrame()
        
        # Initialize list to hold all DataFrames
        dfs_to_combine = [df for df in [gnomad_standardized, clinvar_standardized] if not df.empty]
        
        # Add aggregator data if provided
        if aggregator_csv_path:
            try:
                agg_df = self.get_aggregator_variants(aggregator_csv_path)
                agg_standardized = self.standardize_mutation_data(agg_df, 'aggregator')
                if not agg_standardized.empty:
                    dfs_to_combine.append(agg_standardized)
            except Exception as e:
                print(f"Error processing aggregator data: {str(e)}")
        
        if not dfs_to_combine:
            return pd.DataFrame()
        
        # Combine all DataFrames
        combined_df = pd.concat(dfs_to_combine, ignore_index=True)
        
        # Remove duplicates (same position, ref, and alt)
        combined_df = combined_df.drop_duplicates(
            subset=['position', 'reference', 'alternate'],
            keep='first'
        )
        
        # Calculate relative position based on transcript start
        # This would be helpful for visualization
        if 'position' in combined_df.columns and combined_df['position'].notna().any():
            min_pos = combined_df['position'].min()
            combined_df['relative_position'] = combined_df['position'] - min_pos
        
        return combined_df