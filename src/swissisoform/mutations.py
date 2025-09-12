"""Mutation data handling and analysis module.

This module provides the MutationHandler class for fetching, processing, and
analyzing mutation data from various sources including gnomAD and ClinVar.
"""

from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport
import pandas as pd
from typing import Dict, Optional, List, Any, Set
import requests
import json
import asyncio
import logging
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)


class MutationHandler:
    """Handles variant data from multiple mutation databases.
    
    Provides methods to fetch, parse, and analyze mutation data from:
    - gnomAD
    - ClinVar  
    - COSMIC
    """
    
    def __init__(self, api_key: Optional[str] = "26214372a9ab53b26a43ba8546f4a5195308"):
        """Initialize the mutation handler with API endpoints.
        
        Args:
            api_key: Optional API key for NCBI E-utilities
        """
        self._setup_gnomad_client()
        self._setup_clinvar_endpoints(api_key)
        self._setup_mutation_categories()
        self.cached_data = {}
    
    def _setup_gnomad_client(self):
        """Initialize gnomAD GraphQL client."""
        transport = AIOHTTPTransport(url="https://gnomad.broadinstitute.org/api")
        self.gnomad_client = Client(
            transport=transport, 
            fetch_schema_from_transport=True
        )
    
    def _setup_clinvar_endpoints(self, api_key: str):
        """Initialize ClinVar endpoints and API key."""
        self.clinvar_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.clinvar_endpoints = {
            "esearch": f"{self.clinvar_base}/esearch.fcgi",
            "esummary": f"{self.clinvar_base}/esummary.fcgi", 
            "efetch": f"{self.clinvar_base}/efetch.fcgi",
            "elink": f"{self.clinvar_base}/elink.fcgi",
        }
        self.api_key = api_key
    
    def _setup_mutation_categories(self):
        """Initialize mutation impact categories mapping."""
        self.mutation_categories = {
            "missense": ["missense", "missense_variant", "missense variant"],
            "synonymous": ["synonymous", "synonymous_variant", "synonymous variant"], 
            "nonsense": ["nonsense", "stop_gained", "stop gained", "nonsense_variant"],
            "inframe_del": ["inframe_deletion", "inframe deletion", "inframe_del", "inframe_indel"],
            "inframe_ins": ["inframe_insertion", "inframe insertion", "inframe_ins"],
            "frameshift": ["frameshift", "frameshift_variant", "frameshift variant"],
            "splice": ["splice", "splice_variant", "splice_region_variant", "splice_donor_variant", "splice_acceptor_variant", "splice donor variant", "splice acceptor variant"],
            "start_lost": ["start_lost", "start lost", "start_loss"],
            "utr_5prime": ["5_prime_utr", "5 prime utr", "5_prime_utr_variant", "5 prime UTR variant"],
            "utr_3prime": ["3_prime_utr", "3 prime utr", "3_prime_utr_variant", "3 prime UTR variant"],
            "intronic": ["intron_variant", "intron variant"],
            "stop_gained": ["stop_gained", "stop gained"],
            "stop_lost": ["stop_lost", "stop lost"],
        }
    
    async def get_gnomad_variants(self, gene_name: str) -> pd.DataFrame:
        """Get processed variant data from gnomAD.
        
        Args:
            gene_name: Gene symbol to query
            
        Returns:
            DataFrame of gnomAD variants
        """
        gnomad_data = await self.fetch_gnomad_data(gene_name)
        return self.process_gnomad_variants(gnomad_data)
    
    async def get_gnomad_summary(self, gene_name: str) -> Dict:
        """Get summary statistics for gnomAD variants.
        
        Args:
            gene_name: Gene symbol to query
            
        Returns:
            Dictionary of summary statistics
        """
        gnomad_df = await self.get_gnomad_variants(gene_name)

        # Log detailed raw data structure
        if not gnomad_df.empty:
            self._log_raw_dataframe_structure(gnomad_df, "gnomAD", gene_name)
        else:
            logger.info(f"No gnomAD variants found for {gene_name}")
        
        # Log raw data before processing
        if not gnomad_df.empty:
            logger.info(f"Raw gnomAD data for {gene_name}:")
            logger.info(f"  Total variants: {len(gnomad_df)}")
            logger.info(f"  Columns: {list(gnomad_df.columns)}")
            logger.info(f"  Consequence types: {gnomad_df['consequence'].value_counts().to_dict()}")
            logger.info(f"  Sample rows:\n{gnomad_df.head()}")
        else:
            logger.info(f"No gnomAD variants found for {gene_name}")
        
        if gnomad_df.empty:
            return self._get_empty_gnomad_summary()
        
        return {
            "total_variants": len(gnomad_df),
            "consequence_types": gnomad_df["consequence"].value_counts().to_dict(),
            "mean_allele_frequency": gnomad_df["allele_frequency"].mean(),
            "variants_by_impact": {
                "with_hgvsc": gnomad_df["hgvsc"].notna().sum(),
                "with_hgvsp": gnomad_df["hgvsp"].notna().sum(),
            },
        }
    
    def _get_empty_gnomad_summary(self) -> Dict:
        """Return empty gnomAD summary structure."""
        return {
            "total_variants": 0,
            "consequence_types": {},
            "mean_allele_frequency": 0.0,
            "variants_by_impact": {"with_hgvsc": 0, "with_hgvsp": 0},
        }
    
    async def fetch_gnomad_data(
        self, gene_name: str, reference_genome: str = "GRCh38"
    ) -> Dict:
        """Fetch variant data from gnomAD API for a specific gene using GraphQL.
        
        Args:
            gene_name: Gene symbol (e.g., 'BRCA1')
            reference_genome: Reference genome version (default: 'GRCh38')
            
        Returns:
            Processed gnomAD data for the gene
        """
        cache_key = f"gnomad_{gene_name}_{reference_genome}"
        if cache_key in self.cached_data:
            return self.cached_data[cache_key]
        
        query = self._get_gnomad_query()
        variables = {"geneSymbol": gene_name, "referenceGenome": reference_genome}
        
        try:
            result = await self.gnomad_client.execute_async(
                query, variable_values=variables
            )
            
            if result.get("gene") is not None:
                self.cached_data[cache_key] = result
            return result
            
        except Exception as e:
            logger.error(f"Failed to fetch gnomAD data: {str(e)}")
            raise ConnectionError(f"Failed to fetch gnomAD data: {str(e)}")
    
    def _get_gnomad_query(self):
        """Get the GraphQL query for gnomAD data."""
        return gql("""
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

    def process_gnomad_variants(
        self,
        data: Dict,
        filter_field: Optional[str] = None,
        filter_value: Optional[str] = None,
    ) -> pd.DataFrame:
        """Process gnomAD variant data into a pandas DataFrame.
        
        Args:
            data: Raw gnomAD API response
            filter_field: Field to filter on (e.g., 'consequence')
            filter_value: Value to filter for (e.g., '5_prime_UTR_variant')
            
        Returns:
            DataFrame of processed variant data
        """
        if not self._validate_gnomad_data(data):
            return pd.DataFrame()
        
        variants = data["gene"]["variants"]
        processed_variants = []
        
        for variant in variants:
            if self._should_skip_variant(variant, filter_field, filter_value):
                continue
                
            variant_info = self._process_gnomad_variant(variant)
            if variant_info:
                processed_variants.append(variant_info)
        
        return pd.DataFrame(processed_variants)
    
    def _validate_gnomad_data(self, data: Dict) -> bool:
        """Validate gnomAD data structure."""
        return (
            data and 
            "gene" in data and 
            data["gene"] and 
            "variants" in data["gene"]
        )
    
    def _should_skip_variant(self, variant: Dict, filter_field: str, filter_value: str) -> bool:
        """Check if variant should be skipped based on filters."""
        if filter_field and filter_value:
            variant_field_value = variant.get(filter_field)
            if variant_field_value != filter_value:
                return True
        return False
    
    def _process_gnomad_variant(self, variant: Dict) -> Optional[Dict]:
        """Process a single gnomAD variant."""
        exome_data = variant.get("exome", {})
        genome_data = variant.get("genome", {})
        freq_data = exome_data or genome_data
        
        if not freq_data:
            return None
        
        variant_info = {
            "position": variant["pos"],
            "variant_id": variant["variant_id"],
            "reference": variant["ref"],
            "alternate": variant["alt"],
            "consequence": variant["consequence"],
            "transcript_id": variant.get("transcript_id"),
            "transcript_version": variant.get("transcript_version"),
            "hgvs": variant.get("hgvs"),
            "hgvsc": variant.get("hgvsc"),
            "hgvsp": variant.get("hgvsp"),
            "flags": ",".join(variant.get("flags", [])),
            "allele_count": freq_data.get("ac"),
            "allele_count_hom": freq_data.get("ac_hom"),
            "allele_count_hemi": freq_data.get("ac_hemi"),
            "allele_number": freq_data.get("an"),
            "allele_frequency": freq_data.get("af"),
            "filters": ",".join(freq_data.get("filters", [])),
            "freq_source": "exome" if exome_data else "genome",
        }
        
        # Add population-specific data
        self._add_population_data(variant_info, freq_data)
        return variant_info
    
    def _add_population_data(self, variant_info: Dict, freq_data: Dict):
        """Add population-specific frequency data to variant info."""
        if "populations" in freq_data:
            for pop in freq_data["populations"]:
                pop_id = pop["id"]
                variant_info[f"{pop_id}_ac"] = pop.get("ac")
                variant_info[f"{pop_id}_an"] = pop.get("an")
                variant_info[f"{pop_id}_ac_hom"] = pop.get("ac_hom")
                variant_info[f"{pop_id}_ac_hemi"] = pop.get("ac_hemi")

    async def get_clinvar_variants(self, gene_name: str) -> pd.DataFrame:
        """Get all ClinVar variants for a gene.
        
        Args:
            gene_name: Gene symbol (e.g., 'BRCA1')
            
        Returns:
            DataFrame of processed ClinVar variant data
        """
        cache_key = f"clinvar_{gene_name}"
        if cache_key in self.cached_data:
            return self.cached_data[cache_key]
        
        try:
            variant_ids = await self._fetch_clinvar_variant_ids(gene_name)
            if not variant_ids:
                return pd.DataFrame()
            
            all_variants = await self._process_clinvar_batches(variant_ids)
            
            if all_variants:
                final_df = pd.concat(all_variants, ignore_index=True)
                self.cached_data[cache_key] = final_df
                return final_df
            
            return pd.DataFrame()
            
        except Exception as e:
            logger.error(f"Error fetching ClinVar data: {str(e)}")
            return pd.DataFrame()
    
    async def _fetch_clinvar_variant_ids(self, gene_name: str) -> List[str]:
        """Fetch ClinVar variant IDs for a gene."""
        search_params = {
            "db": "clinvar",
            "term": f"{gene_name}[gene]",
            "retmax": 5000,
            "retmode": "json",
            "tool": "swissisoform",
            "email": "mdiberna@wi.mit.edu",
        }
        
        if self.api_key:
            search_params["api_key"] = self.api_key
        
        response = await self._async_get_request(
            self.clinvar_endpoints["esearch"], params=search_params
        )
        await asyncio.sleep(0.35)
        
        data = json.loads(response)
        if "esearchresult" not in data or "idlist" not in data["esearchresult"]:
            return []
        
        return data["esearchresult"]["idlist"]
    
    async def _process_clinvar_batches(self, variant_ids: List[str]) -> List[pd.DataFrame]:
        """Process ClinVar variants in batches."""
        batch_size = 50
        all_variants = []
        
        for i in range(0, len(variant_ids), batch_size):
            batch_ids = variant_ids[i : i + batch_size]
            batch_variants = await self._process_clinvar_batch(batch_ids)
            
            if not batch_variants.empty:
                all_variants.append(batch_variants)
        
        return all_variants
    
    async def _process_clinvar_batch(self, batch_ids: List[str]) -> pd.DataFrame:
        """Process a single batch of ClinVar variants."""
        summary_params = {
            "db": "clinvar",
            "id": ",".join(batch_ids),
            "retmode": "json",
            "tool": "swissisoform",
            "email": "mdiberna@wi.mit.edu",
        }
        
        if self.api_key:
            summary_params["api_key"] = self.api_key
        
        summary_response = await self._async_get_request(
            self.clinvar_endpoints["esummary"], params=summary_params
        )
        await asyncio.sleep(0.35)
        
        summary_data = json.loads(summary_response)
        return self._process_clinvar_variants(summary_data)

    def _process_clinvar_variants(self, data: Dict) -> pd.DataFrame:
        """Process ClinVar variant data into a pandas DataFrame.
        
        Args:
            data: Raw ClinVar API response
            
        Returns:
            DataFrame of processed variant data
        """
        if "result" not in data or "uids" not in data["result"]:
            return pd.DataFrame()
        
        variants = []
        for uid in data["result"]["uids"]:
            variant_data = data["result"][uid]
            variant_info = self._extract_clinvar_variant_info(uid, variant_data)
            variants.append(variant_info)
        
        df = pd.DataFrame(variants)
        return self._convert_clinvar_columns(df)
    
    def _extract_clinvar_variant_info(self, uid: str, variant_data: Dict) -> Dict:
        """Extract information from a single ClinVar variant."""
        variant_info = {
            "variant_id": uid,
            "obj_type": variant_data.get("obj_type", ""),
            "accession": variant_data.get("accession", ""),
            "title": variant_data.get("title", ""),
            "gene_name": "",
            "review_status": "",
            "clinical_significance": "",
            "last_evaluated": "",
            "molecular_consequences": "",
            "protein_change": variant_data.get("protein_change", ""),
            "chromosome": "",
            "start": "",
            "stop": "",
            "ref_allele": "",
            "alt_allele": "",
            "submission_count": len(
                variant_data.get("supporting_submissions", {}).get("scv", [])
            ),
        }
        
        # Extract gene information
        self._extract_clinvar_gene_info(variant_info, variant_data)
        
        # Extract clinical significance
        self._extract_clinvar_clinical_significance(variant_info, variant_data)
        
        # Extract molecular consequences
        self._extract_clinvar_molecular_consequences(variant_info, variant_data)
        
        # Extract location information
        self._extract_clinvar_location_info(variant_info, variant_data)
        
        return variant_info
    
    def _extract_clinvar_gene_info(self, variant_info: Dict, variant_data: Dict):
        """Extract gene information from ClinVar variant."""
        genes = variant_data.get("genes", [])
        if genes:
            variant_info["gene_name"] = genes[0].get("symbol", "")
            variant_info["gene_id"] = genes[0].get("geneid", "")
            variant_info["gene_strand"] = genes[0].get("strand", "")
    
    def _extract_clinvar_clinical_significance(self, variant_info: Dict, variant_data: Dict):
        """Extract clinical significance from ClinVar variant."""
        classification_types = [
            "germline_classification",
            "clinical_impact_classification", 
            "oncogenicity_classification",
        ]
        
        for class_type in classification_types:
            if class_type in variant_data:
                classification = variant_data[class_type]
                if classification.get("description"):
                    variant_info["clinical_significance"] = classification.get("description", "")
                    variant_info["review_status"] = classification.get("review_status", "")
                    variant_info["last_evaluated"] = classification.get("last_evaluated", "")
                    break
    
    def _extract_clinvar_molecular_consequences(self, variant_info: Dict, variant_data: Dict):
        """Extract molecular consequences from ClinVar variant."""
        mol_cons = variant_data.get("molecular_consequence_list", [])
        variant_info["molecular_consequences"] = ";".join(mol_cons) if mol_cons else ""
    
    def _extract_clinvar_location_info(self, variant_info: Dict, variant_data: Dict):
        """Extract location information from ClinVar variant."""
        variation_sets = variant_data.get("variation_set", [])
        if not variation_sets:
            return
            
        var_set = variation_sets[0]
        variant_info["variant_name"] = var_set.get("variation_name", "")
        variant_info["cdna_change"] = var_set.get("cdna_change", "")
        variant_info["canonical_spdi"] = var_set.get("canonical_spdi", "")
        
        # Extract location from GRCh38
        locations = var_set.get("variation_loc", [])
        for loc in locations:
            if (loc.get("assembly_name") == "GRCh38" and 
                loc.get("status") == "current"):
                variant_info["chromosome"] = loc.get("chr", "")
                variant_info["start"] = loc.get("start", "")
                variant_info["stop"] = loc.get("stop", "")
                variant_info["cytogenic_location"] = loc.get("band", "")
                break
    
    def _convert_clinvar_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert ClinVar columns to appropriate types."""
        numeric_cols = ["start", "stop", "submission_count"]
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors="ignore")
        return df


    async def get_clinvar_summary(self, gene_name: str) -> Dict:
        """Get comprehensive summary statistics for ClinVar variants.
        
        Args:
            gene_name: Gene symbol
            
        Returns:
            Dictionary of summary statistics
        """
        variants_df = await self.get_clinvar_variants(gene_name)

        # Log detailed raw data structure
        if not variants_df.empty:
            self._log_raw_dataframe_structure(variants_df, "ClinVar", gene_name)
        else:
            logger.info(f"No ClinVar variants found for {gene_name}")
        
        # Log raw data before processing
        if not variants_df.empty:
            logger.info(f"Raw ClinVar data for {gene_name}:")
            logger.info(f"  Total variants: {len(variants_df)}")
            logger.info(f"  Columns: {list(variants_df.columns)}")
            logger.info(f"  Clinical significance types: {variants_df['clinical_significance'].value_counts().to_dict()}")
            logger.info(f"  Molecular consequences: {variants_df['molecular_consequences'].value_counts().to_dict()}")
            logger.info(f"  Sample rows:\n{variants_df.head()}")
        else:
            logger.info(f"No ClinVar variants found for {gene_name}")
        
        if variants_df.empty:
            return self._get_empty_clinvar_summary()
        
        summary = {
            "total_variants": len(variants_df),
            "clinical_significance": variants_df["clinical_significance"].value_counts().to_dict(),
            "variant_types": variants_df["obj_type"].value_counts().to_dict(),
            "molecular_consequences": self._get_molecular_consequences_counts(variants_df),
            "chromosome_distribution": variants_df["chromosome"].value_counts().to_dict(),
            "review_status": variants_df["review_status"].value_counts().to_dict(),
            "metadata": self._get_clinvar_metadata(variants_df),
        }
    
        # Add pathogenicity categorization
        summary["pathogenicity_categories"] = self._get_pathogenicity_categories(variants_df)
        return summary
    
    def _get_empty_clinvar_summary(self) -> Dict:
        """Return empty ClinVar summary structure."""
        return {
            "total_variants": 0,
            "clinical_significance": {},
            "variant_types": {},
            "molecular_consequences": {},
            "chromosome_distribution": {},
            "review_status": {},
            "metadata": {
                "submission_count": {"total": 0, "mean": 0, "max": 0},
                "unique_accessions": 0,
            },
        }
    
    def _get_molecular_consequences_counts(self, variants_df: pd.DataFrame) -> Dict:
        """Get molecular consequences counts."""
        return (variants_df["molecular_consequences"]
                .str.split(";")
                .explode()
                .value_counts()
                .to_dict())
    
    def _get_clinvar_metadata(self, variants_df: pd.DataFrame) -> Dict:
        """Get ClinVar metadata statistics."""
        return {
            "submission_count": {
                "total": variants_df["submission_count"].sum(),
                "mean": round(variants_df["submission_count"].mean(), 2),
                "max": variants_df["submission_count"].max(),
            },
            "unique_accessions": variants_df["accession"].nunique(),
        }
    
    def _get_pathogenicity_categories(self, variants_df: pd.DataFrame) -> Dict:
        """Get pathogenicity category counts."""
        pathogenicity_mapping = {
            "likely_pathogenic_or_pathogenic": (
                variants_df["clinical_significance"].str.contains("pathogenic", case=False, na=False) &
                ~variants_df["clinical_significance"].str.contains("conflicting", case=False, na=False)
            ),
            "benign_or_likely_benign": (
                variants_df["clinical_significance"].str.contains("benign", case=False, na=False) &
                ~variants_df["clinical_significance"].str.contains("conflicting", case=False, na=False)
            ),
            "uncertain_significance": variants_df["clinical_significance"].str.contains(
                "uncertain|conflicting", case=False, na=False
            ),
            "not_provided": (
                variants_df["clinical_significance"].isna() |
                (variants_df["clinical_significance"] == "")
            ),
        }
        
        return {
            category: int(variants_df[mask].shape[0])
            for category, mask in pathogenicity_mapping.items()
        }

    async def get_cosmic_variants(self, gene_name: str) -> pd.DataFrame:
        """Get COSMIC variants for a gene using the downloaded TSV/Parquet database.
        
        Args:
            gene_name: Gene symbol (e.g., 'BRCA1')
            
        Returns:
            DataFrame of processed COSMIC variants
        """
        cache_key = f"cosmic_{gene_name}"
        if cache_key in self.cached_data:
            return self.cached_data[cache_key]
        
        try:
            loop = asyncio.get_event_loop()
            cosmic_data = await loop.run_in_executor(
                None, 
                self._fetch_cosmic_sync, 
                gene_name
            )
            
            self.cached_data[cache_key] = cosmic_data
            return cosmic_data
            
        except Exception as e:
            logger.error(f"Failed to fetch COSMIC data for {gene_name}: {str(e)}")
            return pd.DataFrame()
    
    def _find_cosmic_database(self) -> Dict[str, Optional[str]]:
        """Find COSMIC Parquet files and return paths for different datasets."""
        search_paths = [
            "../data/mutation_data",
            "data/mutation_data", 
            ".",
            "../data/cosmic",
            "data/cosmic",
        ]
        
        cosmic_files = {'cancer_census': None, 'noncoding': None}
        
        for search_path in search_paths:
            path = Path(search_path)
            if not path.exists():
                continue
                
            self._find_cosmic_files_in_path(path, cosmic_files)
        
        return cosmic_files
    
    def _find_cosmic_files_in_path(self, path: Path, cosmic_files: Dict[str, Optional[str]]):
        """Find COSMIC files in a specific path."""
        # Look for Cancer Mutation Census files
        if not cosmic_files['cancer_census']:
            cosmic_files['cancer_census'] = self._find_file_with_fallback(
                path, 
                "*CancerMutationCensus_AllData*.parquet",
                "*CancerMutationCensus_AllData*.tsv",
                "COSMIC Cancer Census"
            )
        
        # Look for NonCoding variants files
        if not cosmic_files['noncoding']:
            cosmic_files['noncoding'] = self._find_file_with_fallback(
                path,
                "*Cosmic_NonCodingVariants*.parquet", 
                "*Cosmic_NonCodingVariants*.tsv",
                "COSMIC NonCoding variants"
            )
    
    def _find_file_with_fallback(self, path: Path, parquet_pattern: str, 
                                tsv_pattern: str, file_type: str) -> Optional[str]:
        """Find file with Parquet preference and TSV fallback."""
        # Try Parquet first
        parquet_files = list(path.glob(parquet_pattern))
        if parquet_files:
            logger.info(f"Found {file_type} (Parquet): {parquet_files[0]}")
            return str(parquet_files[0])
        
        # Fallback to TSV
        tsv_files = list(path.glob(tsv_pattern))
        if tsv_files:
            logger.info(f"Found {file_type} (TSV): {tsv_files[0]}")
            return str(tsv_files[0])
        
        return None
    
    def _fetch_cosmic_sync(self, gene_name: str) -> pd.DataFrame:
        """Synchronous wrapper that combines both COSMIC datasets."""
        try:
            cosmic_files = self._find_cosmic_database()
            combined_results = []
            
            # Process both datasets
            for dataset_name, file_path in cosmic_files.items():
                if file_path:
                    result = self._process_cosmic_dataset(file_path, gene_name, dataset_name)
                    if result is not None and not result.empty:
                        combined_results.append(result)
            
            # Combine results
            if combined_results:
                combined_df = pd.concat(combined_results, ignore_index=True, sort=False)
                logger.info(f"Combined {len(combined_df)} total COSMIC variants for {gene_name}")
                # Print column names
                return self._process_cosmic_variants(combined_df)
            else:
                logger.warning(f"No COSMIC data found for {gene_name}")
                return pd.DataFrame()
                
        except Exception as e:
            logger.error(f"COSMIC data fetch failed: {str(e)}")
            return pd.DataFrame()
    
    def _process_cosmic_dataset(self, file_path: str, gene_name: str, dataset_name: str) -> pd.DataFrame:
        """Process a single COSMIC dataset."""
        try:
            logger.info(f"Querying {dataset_name} for {gene_name}")
            
            if file_path.endswith('.parquet'):
                result = self._query_parquet_by_gene(file_path, gene_name, dataset_name)
            else:
                gene_column = 'GENE_NAME' if dataset_name == 'cancer_census' else 'GENE_SYMBOL'
                result = self._query_tsv_by_gene(file_path, gene_name, gene_column)
            
            if result is not None and not result.empty:
                result['cosmic_dataset'] = dataset_name
                logger.info(f"Found {len(result)} variants in {dataset_name}")
                return result
                
        except Exception as e:
            logger.warning(f"{dataset_name} query failed: {str(e)}")
            
        return pd.DataFrame()

    def _process_cosmic_variants(self, cosmic_df: pd.DataFrame) -> pd.DataFrame:
        """Process raw COSMIC data into standardized format, handling both datasets"""
        if cosmic_df.empty:
            return pd.DataFrame()
        
        processed_variants = []
        
        for _, variant in cosmic_df.iterrows():
            # Extract position from genomic coordinate (GRCh37 or GRCh38)
            position = self._extract_cosmic_position(variant)
            if position is None:
                continue  # Skip variants without clear genomic position
                
            # Base variant info that works for both datasets
            variant_info = {
                'position': position,
                'cosmic_id': variant.get('GENOMIC_MUTATION_ID', variant.get('LEGACY_MUTATION_ID', '')),
                'gene_name': variant.get('GENE_NAME', variant.get('GENE_SYMBOL', '')),  # Handle both datasets
                'accession_number': variant.get('ACCESSION_NUMBER', variant.get('TRANSCRIPT_ACCESSION', '')),
                'chromosome': self._extract_cosmic_chromosome(variant),
                'cosmic_dataset': variant.get('cosmic_dataset', 'unknown'),
                
                # Genomic coordinates - handle both formats
                'hgvs_genomic_grch37': variant.get('Mutation genome position GRCh37', ''),
                'hgvs_genomic_grch38': variant.get('Mutation genome position GRCh38', variant.get('HGVSG', '')),
                
                # Allele information - handle both column name formats
                'genomic_wt_allele': variant.get('GENOMIC_WT_ALLELE_SEQ', variant.get('GENOMIC_WT_ALLELE', '')),
                'genomic_mut_allele': variant.get('GENOMIC_MUT_ALLELE_SEQ', variant.get('GENOMIC_MUT_ALLELE', '')),
                'reference': variant.get('GENOMIC_WT_ALLELE_SEQ', variant.get('GENOMIC_WT_ALLELE', '')),
                'alternate': variant.get('GENOMIC_MUT_ALLELE_SEQ', variant.get('GENOMIC_MUT_ALLELE', '')),
                
                # Common annotations (mostly from Cancer Census)
                'oncogene_tsg': variant.get('ONC_TSG', ''),
                'cgc_tier': variant.get('CGC_TIER', ''),
                'disease': variant.get('DISEASE', ''),
                'clinvar_significance': variant.get('CLINVAR_CLNSIG', ''),
                'clinvar_trait': variant.get('CLINVAR_TRAIT', ''),
                
                # Sample information - different column names between datasets
                'sample_name': variant.get('SAMPLE_NAME', ''),
                'cosmic_sample_id': variant.get('COSMIC_SAMPLE_ID', ''),
                'zygosity': variant.get('ZYGOSITY', ''),
                'somatic_status': variant.get('MUTATION_SOMATIC_STATUS', ''),
                
                # Study information
                'cosmic_study_id': variant.get('COSMIC_STUDY_ID', ''),
                'pubmed_pmid': variant.get('PUBMED_PMID', ''),
                
                # Frequency data (mainly in Cancer Census)
                'gnomad_exomes_af': variant.get('GNOMAD_EXOMES_AF', ''),
                'gnomad_genomes_af': variant.get('GNOMAD_GENOMES_AF', ''),
                
                # Functional scores (mainly in Cancer Census)
                'gerp_score': variant.get('GERP++_RS', ''),
                'sift_score': variant.get('MIN_SIFT_SCORE', ''),
                'sift_prediction': variant.get('MIN_SIFT_PRED', ''),
                
                # URLs
                'mutation_url': variant.get('MUTATION_URL', ''),
            }
            
            # Dataset-specific fields
            if variant.get('cosmic_dataset') == 'cancer_census':
                # Cancer Census specific fields (using actual column names with underscores)
                variant_info.update({
                    'mutation_description': variant.get('Mutation_CDS', ''),
                    'mutation_aa': variant.get('Mutation_AA', ''),
                    'mutation_type_cds': variant.get('Mutation_Description_CDS', ''),
                    'mutation_type_aa': variant.get('Mutation_Description_AA', ''),
                    'consequence': variant.get('Mutation_Description_AA', ''),
                    'hgvs_coding': variant.get('Mutation_CDS', ''),
                    'hgvs_protein': variant.get('Mutation_AA', ''),
                    'aa_mut_start': variant.get('AA_MUT_START', ''),
                    'aa_mut_stop': variant.get('AA_MUT_STOP', ''),
                    'shared_aa': variant.get('SHARED_AA', ''),
                    'cosmic_samples_tested': variant.get('COSMIC_SAMPLE_TESTED', ''),
                    'cosmic_samples_mutated': variant.get('COSMIC_SAMPLE_MUTATED', ''),
                    'mutation_significance_tier': variant.get('MUTATION_SIGNIFICANCE_TIER', ''),
                    'ontology_code': variant.get('ONTOLOGY_MUTATION_CODE', ''),
                    'aa_wt_allele': variant.get('AA_WT_ALLELE_SEQ', ''),
                    'aa_mut_allele': variant.get('AA_MUT_ALLELE_SEQ', ''),
                    'dnds_disease_qval': variant.get('DNDS_DISEASE_QVAL', ''),
                    'wgs_disease': variant.get('WGS_DISEASE', ''),
                })
            else:
                # NonCoding variants - fill with empty values for consistency
                variant_info.update({
                    'mutation_description': '',
                    'mutation_aa': '',
                    'mutation_type_cds': '',
                    'mutation_type_aa': 'noncoding',
                    'consequence': 'noncoding_variant',
                    'hgvs_coding': '',
                    'hgvs_protein': '',
                    'aa_mut_start': '',
                    'aa_mut_stop': '',
                    'shared_aa': '',
                    'cosmic_samples_tested': '',
                    'cosmic_samples_mutated': '',
                    'mutation_significance_tier': '',
                    'ontology_code': '',
                    'aa_wt_allele': '',
                    'aa_mut_allele': '',
                    'dnds_disease_qval': '',
                    'wgs_disease': '',
                    # NonCoding specific fields
                    'mutation_nc_id': variant.get('MUTATION_NC_ID', ''),
                    'cosmic_phenotype_id': variant.get('COSMIC_PHENOTYPE_ID', ''),
                    'genome_stop': variant.get('GENOME_STOP', ''),
                })
            
            processed_variants.append(variant_info)
        
        result_df = pd.DataFrame(processed_variants)
        logger.info(f"Processed {len(result_df)} combined COSMIC variants")
        
        # Add summary of datasets
        if not result_df.empty:
            dataset_counts = result_df['cosmic_dataset'].value_counts()
            logger.info(f"Dataset breakdown: {dataset_counts.to_dict()}")
        
        return result_df
    
    def _extract_cosmic_position(self, variant: pd.Series) -> Optional[int]:
        """Extract genomic position from COSMIC variant data (both datasets)"""
        try:
            # For NonCoding variants: use GENOME_START
            if 'GENOME_START' in variant and pd.notna(variant.get('GENOME_START')):
                return int(variant.get('GENOME_START'))
            
            # For Cancer Census: try the mutation position fields
            pos_fields = [
                'Mutation genome position GRCh38',
                'Mutation genome position GRCh37',
                'Mutation_genome_position_GRCh38',
                'Mutation_genome_position_GRCh37',
            ]
            
            for field in pos_fields:
                if field in variant and pd.notna(variant.get(field)):
                    pos_str = str(variant.get(field)).strip()
                    
                    if not pos_str or pos_str.lower() in ['nan', 'none', '']:
                        continue
                        
                    if ':' in pos_str:
                        try:
                            pos_part = pos_str.split(':')[1]
                            if '-' in pos_part:
                                return int(pos_part.split('-')[0])
                            else:
                                return int(pos_part)
                        except (ValueError, IndexError):
                            continue
                    elif pos_str.isdigit():
                        return int(pos_str)
            
            return None
            
        except Exception:
            return None

    def _extract_cosmic_chromosome(self, variant: pd.Series) -> str:
        """Extract chromosome from COSMIC variant data (both datasets)"""
        
        # For NonCoding variants: direct chromosome column
        if 'CHROMOSOME' in variant and pd.notna(variant.get('CHROMOSOME')):
            return str(variant.get('CHROMOSOME'))
        
        # For Cancer Census: extract from position strings
        pos_fields = [
            'Mutation genome position GRCh38',
            'Mutation genome position GRCh37',
            'Mutation_genome_position_GRCh38',
            'Mutation_genome_position_GRCh37',
        ]
        
        for field in pos_fields:
            if field in variant and pd.notna(variant.get(field)):
                pos_str = str(variant.get(field))
                if ':' in pos_str:
                    return pos_str.split(':')[0]
        
        return ''
    
    def _query_parquet_by_gene(self, file_path: str, gene_name: str, dataset_name: str) -> pd.DataFrame:
        """Query Parquet file by gene name."""
        import pyarrow.parquet as pq
        import pyarrow.compute as pc
        
        table = pq.read_table(file_path)
        gene_column = 'GENE_NAME' if dataset_name == 'cancer_census' else 'GENE_SYMBOL'
        
        mask = pc.match_substring(
            pc.utf8_upper(table[gene_column]), 
            gene_name.upper()
        )
        filtered_table = table.filter(mask)
        return filtered_table.to_pandas()
    
    def _query_tsv_by_gene(self, file_path: str, gene_name: str, gene_column: str) -> pd.DataFrame:
        """Query TSV file by gene name using chunked reading for memory efficiency"""
        import pandas as pd
        
        results = []
        chunk_size = 10000
        
        try:
            for chunk in pd.read_csv(file_path, sep='\t', chunksize=chunk_size, low_memory=False):
                if gene_column in chunk.columns:
                    gene_filter = chunk[gene_column].str.contains(gene_name, case=False, na=False)
                    matches = chunk[gene_filter]
                    if not matches.empty:
                        results.append(matches)
            
            if results:
                return pd.concat(results, ignore_index=True)
            else:
                return pd.DataFrame()
                
        except Exception as e:
            logger.error(f"Error reading TSV file {file_path}: {str(e)}")
            return pd.DataFrame()

    async def get_cosmic_summary(self, gene_name: str) -> Dict:
        """Get summary statistics for COSMIC variants."""
        cosmic_df = await self.get_cosmic_variants(gene_name)
        
        # Log detailed raw data structure
        if not cosmic_df.empty:
            self._log_raw_dataframe_structure(cosmic_df, "COSMIC", gene_name)
        else:
            logger.info(f"No COSMIC variants found for {gene_name}")

        # Log raw data before processing
        if not cosmic_df.empty:
            logger.info(f"Raw COSMIC data for {gene_name}:")
            logger.info(f"  Total variants: {len(cosmic_df)}")
            logger.info(f"  Columns: {list(cosmic_df.columns)}")
            if "mutation_type_cds" in cosmic_df.columns:
                logger.info(f"  Mutation types (CDS): {cosmic_df['mutation_type_cds'].value_counts().to_dict()}")
            if "mutation_type_aa" in cosmic_df.columns:
                logger.info(f"  Mutation types (AA): {cosmic_df['mutation_type_aa'].value_counts().to_dict()}")
            if "disease" in cosmic_df.columns:
                logger.info(f"  Disease distribution: {cosmic_df['disease'].value_counts().head().to_dict()}")
            logger.info(f"  Sample rows:\n{cosmic_df.head()}")
        else:
            logger.info(f"No COSMIC variants found for {gene_name}")
        
        if cosmic_df is None or cosmic_df.empty:
            return self._get_empty_cosmic_summary()
            
        summary = {
            "total_variants": len(cosmic_df),
            "mutation_types_cds": cosmic_df["mutation_type_cds"].value_counts(dropna=False).to_dict() if "mutation_type_cds" in cosmic_df else {},
            "mutation_types_aa": cosmic_df["mutation_type_aa"].value_counts(dropna=False).to_dict() if "mutation_type_aa" in cosmic_df else {},
            "diseases": cosmic_df["disease"].value_counts(dropna=False).to_dict() if "disease" in cosmic_df else {},
            "oncogene_tsg_roles": cosmic_df["oncogene_tsg"].value_counts(dropna=False).to_dict() if "oncogene_tsg" in cosmic_df else {},
            "cgc_tiers": cosmic_df["cgc_tier"].value_counts(dropna=False).to_dict() if "cgc_tier" in cosmic_df else {},
            "mutation_significance_tiers": cosmic_df["mutation_significance_tier"].value_counts(dropna=False).to_dict() if "mutation_significance_tier" in cosmic_df else {},
            "clinvar_significance": cosmic_df["clinvar_significance"].value_counts(dropna=False).to_dict() if "clinvar_significance" in cosmic_df else {},
            "sift_predictions": cosmic_df["sift_prediction"].value_counts(dropna=False).to_dict() if "sift_prediction" in cosmic_df else {},
            "cosmic_samples_data": {
                "total_tested": pd.to_numeric(cosmic_df["cosmic_samples_tested"], errors="coerce").sum() if "cosmic_samples_tested" in cosmic_df else 0,
                "total_mutated": pd.to_numeric(cosmic_df["cosmic_samples_mutated"], errors="coerce").sum() if "cosmic_samples_mutated" in cosmic_df else 0,
            }
        }
        return summary
        
    def _get_empty_cosmic_summary(self) -> Dict:
        """Return empty COSMIC summary structure."""
        return {
            'total_variants': 0,
            'mutation_types_cds': {},
            'mutation_types_aa': {},
            'diseases': {},
            'oncogene_tsg_roles': {},
            'cgc_tiers': {},
            'mutation_significance_tiers': {},
            'clinvar_significance': {},
            'sift_predictions': {},
            'cosmic_samples_data': {'total_tested': 0, 'total_mutated': 0}
        }
    
    def _standardize_hgvs(self, hgvs_str: str) -> str:
        """Standardize HGVS notation by extracting the core variant description.
        
        Args:
            hgvs_str: Raw HGVS string from any source
            
        Returns:
            Standardized HGVS string with just the variant description
        """
        if pd.isna(hgvs_str) or not isinstance(hgvs_str, str):
            return ""
        
        hgvs_str = hgvs_str.strip()
        if not hgvs_str:
            return ""
        
        # Handle ClinVar format: "NM_005228.5(EGFR):c.-216G>T" -> "c.-216G>T"
        if ":" in hgvs_str:
            parts = hgvs_str.split(":")
            if len(parts) >= 2:
                # Take the part after the colon (the actual HGVS description)
                hgvs_str = parts[-1].strip()
        
        # Remove any leading/trailing whitespace and common prefixes
        hgvs_str = hgvs_str.strip()
        
        # Handle protein descriptions that might have extra info after space
        # "p.Gly5Val some_extra_info" -> "p.Gly5Val"
        if hgvs_str.startswith("p.") and " " in hgvs_str:
            hgvs_str = hgvs_str.split(" ")[0]
        
        return hgvs_str

    def standardize_mutation_data(self, variants_df: pd.DataFrame, source: str) -> pd.DataFrame:
        """Standardize mutation data from different sources into a consistent format.
        
        Args:
            variants_df: Raw variants DataFrame from any source
            source: Source of the data ('gnomad', 'clinvar', or 'cosmic')
            
        Returns:
            Standardized mutation data
        """
        if variants_df.empty:
            return self._get_empty_standardized_df()
        
        source_lower = source.lower()
        
        if source_lower == "gnomad":
            standardized_df = self._standardize_gnomad_data(variants_df)
        elif source_lower == "clinvar":
            standardized_df = self._standardize_clinvar_data(variants_df)
        elif source_lower == "cosmic":
            standardized_df = self._standardize_cosmic_data(variants_df)
        else:
            return self._get_empty_standardized_df()
        
        return self._clean_mutation_data(standardized_df)
    

    def _get_empty_standardized_df(self) -> pd.DataFrame:
        """Return empty standardized DataFrame."""
        standard_columns = [
            "position", "variant_id", "reference", "alternate", "source",
            "impact", "hgvsc", "hgvsp", "allele_frequency", "allele_count_hom", "clinical_significance"
        ]
        return pd.DataFrame(columns=standard_columns)
    
    def _standardize_gnomad_data(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Standardize gnomAD data."""
        return pd.DataFrame({
            "position": variants_df["position"],
            "variant_id": variants_df["variant_id"],
            "reference": variants_df["reference"],
            "alternate": variants_df["alternate"],
            "source": "gnomAD",
            "impact": variants_df["consequence"].apply(self._standardize_impact_category),
            "hgvsc": variants_df["hgvsc"].apply(self._standardize_hgvs),
            "hgvsp": variants_df["hgvsp"].apply(self._standardize_hgvs),
            "allele_frequency": variants_df["allele_frequency"],
            "allele_count_hom": variants_df["allele_count_hom"],  
            "clinical_significance": None,
        })
        
    def _standardize_clinvar_data(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Standardize ClinVar data."""
        return pd.DataFrame({
            "position": variants_df["start"],
            "variant_id": variants_df["accession"],
            "reference": variants_df["ref_allele"],
            "alternate": variants_df["alt_allele"],
            "source": "ClinVar",
            "impact": variants_df["molecular_consequences"].apply(self._standardize_impact_category),
            "hgvsc": variants_df["cdna_change"].apply(self._standardize_hgvs),
            "hgvsp": variants_df["protein_change"].apply(self._standardize_hgvs),
            "allele_frequency": None,
            "allele_count_hom": None,  
            "clinical_significance": variants_df["clinical_significance"],
        })
        
    def _standardize_cosmic_data(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Standardize COSMIC data."""
        return pd.DataFrame({
            "position": variants_df["position"],
            "variant_id": variants_df["cosmic_id"],
            "reference": variants_df["reference"],
            "alternate": variants_df["alternate"],
            "source": "COSMIC",
            "impact": variants_df["consequence"].apply(self._standardize_cosmic_impact),
            "hgvsc": variants_df["hgvs_coding"].apply(self._standardize_hgvs),
            "hgvsp": variants_df["hgvs_protein"].apply(self._standardize_hgvs),
            "allele_frequency": pd.to_numeric(variants_df["gnomad_exomes_af"], errors='coerce'),
            "allele_count_hom": None,  # COSMIC doesn't have this field
            "clinical_significance": variants_df["clinvar_significance"],
        })
    
    def _standardize_impact_category(self, impact: str) -> str:
        """Standardize impact categories to a consistent set of terms."""
        if pd.isna(impact) or not isinstance(impact, str):
            return "unknown"
        
        impact = impact.lower().strip()
        
        # Check each category
        for category, terms in self.mutation_categories.items():
            if any(term in impact for term in terms):
                return self._get_standardized_impact_name(category)
        
        return "other variant"
    
    def _get_standardized_impact_name(self, category: str) -> str:
        """Get standardized impact name for category."""
        mapping = {
            "missense": "missense variant",
            "synonymous": "synonymous variant", 
            "nonsense": "nonsense variant",
            "inframe_del": "inframe deletion",
            "inframe_ins": "inframe insertion",
            "frameshift": "frameshift variant",
            "splice": "splice variant",
            "start_lost": "start lost variant",
            "utr_5prime": "5 prime UTR variant", 
            "utr_3prime": "3 prime UTR variant",
            "intronic": "intron variant",
            "stop_gained": "stop gained variant",
            "stop_lost": "stop lost variant",
        }
        return mapping.get(category, "other variant")
    
    def _standardize_cosmic_impact(self, cosmic_impact: str) -> str:
        """Convert COSMIC mutation classifications to standardized impact terms."""
        if pd.isna(cosmic_impact) or not isinstance(cosmic_impact, str):
            return "unknown"
        
        cosmic_impact = cosmic_impact.lower().strip()
        
        # Map COSMIC descriptions to standard terms
        cosmic_mapping = {
            'missense': 'missense variant',
            'nonsense': 'nonsense variant',
            'stop': 'nonsense variant',
            'frameshift': 'frameshift variant',
            'synonymous': 'synonymous variant',
            'silent': 'synonymous variant',
            'splice': 'splice variant',
            'noncoding': 'noncoding variant',
            'deletion': 'deletion',
            'insertion': 'insertion',
            'substitution': 'substitution',
        }
        
        for key, value in cosmic_mapping.items():
            if key in cosmic_impact:
                return value
        
        # Handle inframe variants
        if 'inframe' in cosmic_impact:
            if 'deletion' in cosmic_impact:
                return 'inframe deletion'
            elif 'insertion' in cosmic_impact:
                return 'inframe insertion'
            else:
                return 'inframe variant'
        
        return 'other variant'
    
    def _clean_mutation_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Clean and validate the standardized mutation data."""
        # Convert position to integer if possible
        df["position"] = pd.to_numeric(df["position"], errors="coerce")
        
        # Clean frequency data
        df["allele_frequency"] = pd.to_numeric(df["allele_frequency"], errors="coerce")
        
        # Fill missing values
        df["clinical_significance"] = df["clinical_significance"].fillna("Unknown")
        df["impact"] = df["impact"].fillna("Unknown")
        
        # Clean string columns
        string_columns = ["hgvsc", "hgvsp", "reference", "alternate"]
        for col in string_columns:
            df[col] = df[col].fillna(".")
            df[col] = df[col].astype(str)
        
        # Sort by position if available
        if df["position"].notna().any():
            df = df.sort_values("position")
        
        return df.reset_index(drop=True)
    
    def _log_raw_dataframe_structure(self, df: pd.DataFrame, source: str, gene_name: str):
        """Log detailed information about raw DataFrame structure."""
        logger.info(f"=== Raw {source} DataFrame structure for {gene_name} ===")
        logger.info(f"Shape: {df.shape}")
        logger.info(f"Columns ({len(df.columns)}): {list(df.columns)}")
        
        # Log data types
        logger.info("Data types:")
        for col in df.columns:
            dtype = df[col].dtype
            non_null_count = df[col].notna().sum()
            logger.info(f"  {col}: {dtype} ({non_null_count}/{len(df)} non-null)")
        
        # Log sample of actual data for key columns
        if not df.empty:
            sample_size = min(3, len(df))
            logger.info(f"Sample data (first {sample_size} rows):")
            for idx in range(sample_size):
                logger.info(f"  Row {idx}:")
                for col in df.columns:
                    value = df.iloc[idx][col]
                    if pd.isna(value):
                        logger.info(f"    {col}: <NA>")
                    else:
                        # Truncate long values
                        str_val = str(value)
                        if len(str_val) > 100:
                            str_val = str_val[:97] + "..."
                        logger.info(f"    {col}: {str_val}")
        logger.info("=" * 50)
    
    async def _async_get_request(self, url: str, params: Dict = None) -> str:
        """Make an async GET request.
        
        Args:
            url: URL to request
            params: Query parameters
            
        Returns:
            Response text
        """
        if params is None:
            params = {}
        
        if self.api_key:
            params["api_key"] = self.api_key
        
        loop = asyncio.get_event_loop()
        
        def _make_request():
            response = requests.get(url, params=params)
            response.raise_for_status()
            return response.text
        
        return await loop.run_in_executor(None, _make_request)
    
    def clear_cache(self):
        """Clear the cached API responses."""
        self.cached_data = {}
    
    async def get_visualization_ready_mutations(
        self,
        gene_name: str,
        alt_features: Optional[pd.DataFrame] = None,
        sources: Optional[List[str]] = None,
        aggregator_csv_path: Optional[str] = None,
        verbose: bool = False,
    ) -> pd.DataFrame:
        """Get mutation data from specified sources in a format ready for visualization.
        
        Args:
            gene_name: Name of the gene
            alt_features: DataFrame containing alternative isoform features
            sources: List of mutation sources to query
            aggregator_csv_path: Path to aggregator CSV file (unused in current implementation)
            verbose: Whether to print detailed information
            
        Returns:
            DataFrame containing mutations ready for visualization
        """
        if sources is None:
            sources = ["gnomad", "clinvar", "cosmic"]
        sources = [s.lower() for s in sources]
        
        # Fetch and combine data from sources
        dfs_to_combine = []
        if verbose:
            print(f"Fetching mutations from sources: {', '.join(sources)}...")
        
        for source in sources:
            try:
                df = await self._fetch_source_data(source, gene_name)
                standardized_df = self.standardize_mutation_data(df, source)
                
                if not standardized_df.empty:
                    dfs_to_combine.append(standardized_df)
                    if verbose:
                        print(f"  Found {len(standardized_df)} variants from {source}")
                        
            except Exception as e:
                logger.error(f"Error fetching {source} data: {str(e)}")
        
        if not dfs_to_combine:
            if verbose:
                print("No mutations found")
            return pd.DataFrame()
        
        # Process combined data
        combined_df = pd.concat(dfs_to_combine, ignore_index=True)
        combined_df["position"] = pd.to_numeric(combined_df["position"], errors="coerce")
        
        # Filter by alt_features if provided
        if alt_features is not None and not alt_features.empty:
            return self._filter_by_alt_features(combined_df, alt_features, verbose)
        else:
            # No filtering, just add relative positions
            if len(combined_df) > 0:
                min_pos = combined_df["position"].min()
                combined_df["relative_position"] = combined_df["position"] - min_pos
            return combined_df
    
    async def _fetch_source_data(self, source: str, gene_name: str) -> pd.DataFrame:
        """Fetch data from a specific source."""
        if source == "gnomad":
            return await self.get_gnomad_variants(gene_name)
        elif source == "clinvar":
            return await self.get_clinvar_variants(gene_name)
        elif source == "cosmic":
            return await self.get_cosmic_variants(gene_name)
        else:
            return pd.DataFrame()
    
    def _filter_by_alt_features(self, combined_df: pd.DataFrame, alt_features: pd.DataFrame, verbose: bool) -> pd.DataFrame:
        """Filter mutations by alternative features."""
        alt_features["start"] = pd.to_numeric(alt_features["start"])
        alt_features["end"] = pd.to_numeric(alt_features["end"])
        
        if verbose:
            print("\nAnalyzing mutations in alternative features:")
        
        filtered_dfs = []
        
        for idx, feature in alt_features.iterrows():
            start_pos = feature["start"]
            end_pos = feature["end"]
            name = feature.get("name", f"Feature {idx + 1}")
            
            # Filter mutations for this feature
            feature_mask = (
                (combined_df["position"] >= start_pos) & 
                (combined_df["position"] <= end_pos)
            )
            feature_mutations = combined_df[feature_mask].copy()
            
            # Add feature information
            if not feature_mutations.empty:
                feature_mutations["alt_feature_region"] = f"{int(start_pos)}-{int(end_pos)}"
                feature_mutations["alt_feature_name"] = name
                filtered_dfs.append(feature_mutations)
            
            if verbose:
                print(f"  Feature {name} ({start_pos}-{end_pos}): {len(feature_mutations)} mutations")
        
        if filtered_dfs:
            filtered_df = pd.concat(filtered_dfs, ignore_index=True)
            if verbose:
                print(f"\nTotal mutations in all features: {len(filtered_df)}")
            
            # Calculate relative positions
            if len(filtered_df) > 0:
                min_pos = filtered_df["position"].min()
                filtered_df["relative_position"] = filtered_df["position"] - min_pos
            return filtered_df
        else:
            if verbose:
                print("No mutations found in any feature")
            return pd.DataFrame()

    async def analyze_gene_mutations_comprehensive(
        self,
        gene_name: str,
        genome_handler,
        alt_isoform_handler,
        output_dir: str,
        visualize: bool = False,
        impact_types: Optional[Dict[str, List[str]]] = None,
        preferred_transcripts: Optional[Set[str]] = None,
    ) -> Dict:
        """Comprehensive mutation analysis for a gene with transcript-truncation pairs.

        This method performs the same analysis as the process_gene function in
        analyze_truncations.py but within the MutationHandler class.

        Args:
            gene_name: Name of the gene to analyze
            genome_handler: Initialized GenomeHandler instance
            alt_isoform_handler: Initialized AlternativeIsoform instance
            output_dir: Directory to save output files
            visualize: Whether to generate visualizations
            impact_types: Dict mapping sources to impact types to filter by
            preferred_transcripts: Set of transcript IDs to prioritize

        Returns:
            Dictionary containing comprehensive analysis results
        """
        try:
            # Get alternative isoform features
            print(f"  ├─ Getting alternative features...", end="", flush=True)
            alt_features = alt_isoform_handler.get_visualization_features(gene_name)
            if alt_features.empty:
                print(f"\r  ├─ No alternative features found")
                return {"gene_name": gene_name, "status": "no_features", "error": None}
            print(f"\r  ├─ Found {len(alt_features)} alternative features")

            # Get transcript information
            print(f"  ├─ Getting transcript information...", end="", flush=True)
            transcript_info = genome_handler.get_transcript_ids(gene_name)
            if transcript_info.empty:
                print(f"\r  ├─ No transcript info found")
                return {
                    "gene_name": gene_name,
                    "status": "no_transcripts",
                    "error": None,
                }

            # Filter by preferred transcripts if provided
            original_transcript_count = len(transcript_info)
            if preferred_transcripts and not transcript_info.empty:
                # First try exact matches
                filtered = transcript_info[
                    transcript_info["transcript_id"].isin(preferred_transcripts)
                ]

                # If nothing matched and versions might be present, try matching base IDs
                if filtered.empty and any("." in t for t in preferred_transcripts):
                    base_preferred = {t.split(".")[0] for t in preferred_transcripts}
                    filtered = transcript_info[
                        transcript_info["transcript_id"]
                        .str.split(".", expand=True)[0]
                        .isin(base_preferred)
                    ]

                # If we found matches, use the filtered set
                if not filtered.empty:
                    transcript_info = filtered
                    print(
                        f"\r  ├─ Filtered to {len(transcript_info)} preferred transcripts (out of {original_transcript_count})"
                    )
                else:
                    print(
                        f"\r  ├─ No preferred transcripts found for gene {gene_name}, using all {len(transcript_info)} transcripts"
                    )

            # Create transcript-truncation pairs based on overlap
            print(f"\r  ├─ Creating transcript-truncation pairs...", end="", flush=True)
            transcript_truncation_pairs = []

            for _, transcript in transcript_info.iterrows():
                transcript_id = transcript["transcript_id"]
                transcript_start = transcript["start"]
                transcript_end = transcript["end"]
                transcript_chromosome = transcript["chromosome"]
                transcript_strand = transcript["strand"]

                # Check each truncation feature for overlap with this transcript
                for idx, truncation in alt_features.iterrows():
                    trunc_start = truncation["start"]
                    trunc_end = truncation["end"]
                    trunc_chrom = truncation["chromosome"]

                    # Skip if chromosomes don't match
                    if transcript_chromosome != trunc_chrom:
                        continue

                    # Check for overlap
                    if not (
                        transcript_end < trunc_start or transcript_start > trunc_end
                    ):
                        # Create an entry for this transcript-truncation pair
                        truncation_id = f"trunc_{idx}"

                        # If we have more identifiable information about the truncation, use it
                        if "start_codon" in truncation and not pd.isna(
                            truncation["start_codon"]
                        ):
                            truncation_id = f"trunc_{truncation['start_codon']}_{trunc_start}_{trunc_end}"

                        transcript_truncation_pairs.append(
                            {
                                "transcript_id": transcript_id,
                                "truncation_id": truncation_id,
                                "truncation_idx": idx,
                                "transcript_start": transcript_start,
                                "transcript_end": transcript_end,
                                "transcript_strand": transcript_strand,
                                "truncation_start": trunc_start,
                                "truncation_end": trunc_end,
                            }
                        )

            if not transcript_truncation_pairs:
                print(f"\r  ├─ No transcripts overlap with truncation regions")
                return {
                    "gene_name": gene_name,
                    "status": "no_overlapping_transcripts",
                    "error": None,
                }

            print(
                f"\r  ├─ Found {len(transcript_truncation_pairs)} transcript-truncation pairs across {len(transcript_info)} transcripts"
            )

            # Get the list of desired impact types for each source
            desired_impact_types = []
            if impact_types:
                for source, impacts in impact_types.items():
                    desired_impact_types.extend(impacts)

            # Fetch raw mutations and filter them for each transcript-truncation pair
            print(f"  ├─ Fetching mutation data...", end="", flush=True)

            # Get raw mutation data from ClinVar (we're not filtering yet)
            all_mutations = await self.get_visualization_ready_mutations(
                gene_name=gene_name,
                alt_features=None,  # Don't filter by alt_features yet
                sources=["clinvar"],
                aggregator_csv_path=None,
            )

            if all_mutations is None or all_mutations.empty:
                print(f"\r  ├─ No mutations found for this gene")
                all_mutations = pd.DataFrame()
            else:
                print(
                    f"\r  ├─ Found {len(all_mutations)} total mutations for this gene"
                )

                # Apply impact type filtering if specified
                if impact_types:
                    print(f"  ├─ Filtering for impact types: {impact_types}")
                    for source, impacts in impact_types.items():
                        source_mutations = all_mutations[
                            all_mutations["source"].str.lower() == source.lower()
                        ]

                        if not source_mutations.empty:
                            filtered_mutations = source_mutations[
                                source_mutations["impact"].isin(impacts)
                            ]
                            # Replace the original mutations with filtered ones
                            all_mutations = all_mutations[
                                all_mutations["source"].str.lower() != source.lower()
                            ]
                            all_mutations = pd.concat(
                                [all_mutations, filtered_mutations]
                            )

            # Container for all transcript-truncation analysis results
            pair_results = []

            # Process each transcript-truncation pair
            print(f"  ├─ Analyzing mutations for each transcript-truncation pair...")

            for pair_idx, pair in enumerate(transcript_truncation_pairs, 1):
                transcript_id = pair["transcript_id"]
                truncation_id = pair["truncation_id"]
                truncation_idx = pair["truncation_idx"]

                # Get truncation start and end positions
                trunc_start = pair["truncation_start"]
                trunc_end = pair["truncation_end"]

                # Initialize default counts for all desired impact types
                mutation_categories = {}

                # Ensure all our desired impact types have columns, even if zero
                for impact_type in desired_impact_types:
                    category_key = f"mutations_{impact_type.replace(' ', '_').lower()}"
                    mutation_categories[category_key] = 0

                    # Also create empty columns for ClinVar IDs for each impact type
                    impact_key = f"clinvar_ids_{impact_type.replace(' ', '_').lower()}"
                    mutation_categories[impact_key] = ""

                # Create a truncation-specific filter for mutations
                if not all_mutations.empty:
                    # Filter mutations to only those in this truncation region
                    pair_mutations = all_mutations[
                        (all_mutations["position"] >= trunc_start)
                        & (all_mutations["position"] <= trunc_end)
                    ].copy()

                    pair_mutation_count = len(pair_mutations)

                    if pair_mutation_count > 0:
                        print(
                            f"  │  ├─ {transcript_id} × {truncation_id}: {pair_mutation_count} mutations"
                        )

                        clinvar_ids = []

                        # Count by impact category
                        for impact in pair_mutations["impact"].unique():
                            if pd.isna(impact):
                                continue

                            # Get mutations for this impact
                            impact_mutations = pair_mutations[
                                pair_mutations["impact"] == impact
                            ]
                            category_count = len(impact_mutations)

                            # Store count for this impact type
                            category_key = (
                                f"mutations_{impact.replace(' ', '_').lower()}"
                            )
                            mutation_categories[category_key] = category_count

                            # Store ClinVar IDs for this impact type
                            if "variant_id" in impact_mutations.columns:
                                impact_ids = (
                                    impact_mutations["variant_id"]
                                    .dropna()
                                    .unique()
                                    .tolist()
                                )
                                # Convert to strings (to handle numeric IDs), filter empty strings
                                impact_ids = [
                                    str(id).strip()
                                    for id in impact_ids
                                    if str(id).strip()
                                ]

                                # Add to the category-specific IDs
                                if impact_ids:
                                    impact_key = f"clinvar_ids_{impact.replace(' ', '_').lower()}"
                                    mutation_categories[impact_key] = ",".join(
                                        impact_ids
                                    )

                        # Collect all ClinVar IDs for this truncation region
                        if "variant_id" in pair_mutations.columns:
                            clinvar_ids = (
                                pair_mutations["variant_id"].dropna().unique().tolist()
                            )
                            clinvar_ids = [
                                str(id).strip() for id in clinvar_ids if str(id).strip()
                            ]

                        # Add results for this pair with detailed mutation categories
                        pair_results.append(
                            {
                                "transcript_id": transcript_id,
                                "truncation_id": truncation_id,
                                "truncation_start": trunc_start,
                                "truncation_end": trunc_end,
                                "mutation_count_total": pair_mutation_count,
                                "clinvar_variant_ids": ",".join(clinvar_ids)
                                if clinvar_ids
                                else "",
                                **mutation_categories,
                            }
                        )
                    else:
                        # No mutations for this pair, still record it with zeros
                        pair_results.append(
                            {
                                "transcript_id": transcript_id,
                                "truncation_id": truncation_id,
                                "truncation_start": trunc_start,
                                "truncation_end": trunc_end,
                                "mutation_count_total": 0,
                                "clinvar_variant_ids": "",
                                **mutation_categories,  # Include zero counts for all categories
                            }
                        )
                else:
                    # No mutations at all, record with zeros
                    pair_results.append(
                        {
                            "transcript_id": transcript_id,
                            "truncation_id": truncation_id,
                            "truncation_start": trunc_start,
                            "truncation_end": trunc_end,
                            "mutation_count_total": 0,
                            "clinvar_variant_ids": "",
                            **mutation_categories,  # Include zero counts for all categories
                        }
                    )

            # Generate visualizations if requested
            if visualize:
                from swissisoform.visualize import GenomeVisualizer

                visualizer = GenomeVisualizer(genome_handler)
                gene_dir = Path(output_dir) / gene_name
                gene_dir.mkdir(parents=True, exist_ok=True)

                print(
                    f"  ├─ Generating visualizations for each transcript-truncation pair:"
                )

                for pair_idx, pair in enumerate(transcript_truncation_pairs, 1):
                    transcript_id = pair["transcript_id"]
                    truncation_id = pair["truncation_id"]
                    truncation_idx = pair["truncation_idx"]

                    # Get the specific truncation for this pair
                    if truncation_idx in alt_features.index:
                        truncation_feature = alt_features.loc[[truncation_idx]].copy()
                    else:
                        print(
                            f"  │  ├─ Warning: Invalid truncation index {truncation_idx}, skipping visualization"
                        )
                        continue

                    print(
                        f"  │  ├─ Visualizing pair {pair_idx}/{len(transcript_truncation_pairs)}: {transcript_id} × {truncation_id}"
                    )

                    # Create directories organized by transcript and truncation
                    transcript_dir = gene_dir / transcript_id
                    transcript_dir.mkdir(exist_ok=True)

                    # Prepare output paths
                    pair_base_filename = f"{transcript_id}_{truncation_id}"

                    # Filter mutations for this specific truncation
                    if not all_mutations.empty:
                        trunc_start = pair["truncation_start"]
                        trunc_end = pair["truncation_end"]

                        pair_mutations = all_mutations[
                            (all_mutations["position"] >= trunc_start)
                            & (all_mutations["position"] <= trunc_end)
                        ].copy()

                        # Create the visualization using only this truncation feature
                        print(
                            f"  │  │  ├─ Creating view with {len(pair_mutations)} mutations"
                        )
                        visualizer.visualize_transcript(
                            gene_name=gene_name,
                            transcript_id=transcript_id,
                            alt_features=truncation_feature,
                            mutations_df=pair_mutations,
                            output_file=str(
                                transcript_dir / f"{pair_base_filename}_filtered.pdf"
                            ),
                        )

                        print(f"  │  │  └─ Creating zoomed view")
                        visualizer.visualize_transcript_zoomed(
                            gene_name=gene_name,
                            transcript_id=transcript_id,
                            alt_features=truncation_feature,
                            mutations_df=pair_mutations,
                            output_file=str(
                                transcript_dir
                                / f"{pair_base_filename}_filtered_zoom.pdf"
                            ),
                            padding=100,
                        )

            # Calculate total mutations across all pairs
            total_mutations = (
                sum(pair["mutation_count_total"] for pair in pair_results)
                if pair_results
                else 0
            )

            print("  └─ Processing complete")
            return {
                "gene_name": gene_name,
                "status": "success",
                "total_transcripts": len(transcript_info),
                "truncation_features": len(alt_features),
                "transcript_truncation_pairs": len(transcript_truncation_pairs),
                "mutations_filtered": total_mutations,
                "pair_results": pair_results,
                "error": None,
            }

        except Exception as e:
            logger.error(f"Error processing gene {gene_name}: {str(e)}")
            print(f"  └─ Error: {str(e)}")
            return {"gene_name": gene_name, "status": "error", "error": str(e)}
