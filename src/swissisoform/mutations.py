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
import gget
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)


class MutationHandler:
    """Handles variant data from multiple mutation databases.

    Provides methods to fetch, parse, and analyze mutation data from:
    - gnomAD
    - ClinVar
    """

    def __init__(self, api_key: Optional[str] = "26214372a9ab53b26a43ba8546f4a5195308"):
        """Initialize the mutation handler with API endpoints.

        Args:
            api_key: Optional API key for NCBI E-utilities
        """
        # Setup gnomAD GraphQL client
        transport = AIOHTTPTransport(url="https://gnomad.broadinstitute.org/api")
        self.gnomad_client = Client(
            transport=transport, fetch_schema_from_transport=True
        )
        self.cached_data = {}

        # ClinVar E-utilities base URL and endpoints
        self.clinvar_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.clinvar_endpoints = {
            "esearch": f"{self.clinvar_base}/esearch.fcgi",
            "esummary": f"{self.clinvar_base}/esummary.fcgi",
            "efetch": f"{self.clinvar_base}/efetch.fcgi",
            "elink": f"{self.clinvar_base}/elink.fcgi",
        }
        self.api_key = api_key

        # Basic mutation impact colors matching standardized categories
        self.mutation_categories = {
            "missense": ["missense", "missense_variant"],
            "synonymous": ["synonymous", "synonymous_variant"],
            "nonsense": ["nonsense", "stop_gained", "stop gained", "nonsense_variant"],
            "inframe_del": ["inframe_deletion", "inframe deletion", "inframe_del"],
            "inframe_ins": ["inframe_insertion", "inframe insertion", "inframe_ins"],
            "frameshift": ["frameshift", "frameshift_variant"],
            "splice": ["splice", "splice_variant", "splice_region_variant"],
            "start_lost": ["start_lost", "start lost", "start_loss"],
            "utr_5prime": ["5_prime_utr", "5 prime utr", "5_prime_utr_variant"],
            "utr_3prime": ["3_prime_utr", "3 prime utr", "3_prime_utr_variant"],
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

        if gnomad_df.empty:
            return {
                "total_variants": 0,
                "consequence_types": {},
                "mean_allele_frequency": 0.0,
                "variants_by_impact": {"with_hgvsc": 0, "with_hgvsp": 0},
            }

        summary = {
            "total_variants": len(gnomad_df),
            "consequence_types": gnomad_df["consequence"].value_counts().to_dict(),
            "mean_allele_frequency": gnomad_df["allele_frequency"].mean(),
            "variants_by_impact": {
                "with_hgvsc": gnomad_df["hgvsc"].notna().sum(),
                "with_hgvsp": gnomad_df["hgvsp"].notna().sum(),
            },
        }
        return summary

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
        variables = {"geneSymbol": gene_name, "referenceGenome": reference_genome}

        try:
            # Execute the query
            result = await self.gnomad_client.execute_async(
                query, variable_values=variables
            )

            # Cache the response
            if result.get("gene") is not None:
                self.cached_data[cache_key] = result
            return result

        except Exception as e:
            logger.error(f"Failed to fetch gnomAD data: {str(e)}")
            raise ConnectionError(f"Failed to fetch gnomAD data: {str(e)}")

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
        if (
            not data
            or "gene" not in data
            or not data["gene"]
            or "variants" not in data["gene"]
        ):
            return pd.DataFrame()

        variants = data["gene"]["variants"]
        processed_variants = []

        for variant in variants:
            # Skip this variant if it doesn't match the filter condition
            if filter_field and filter_value:
                variant_field_value = variant.get(filter_field)
                if variant_field_value != filter_value:
                    continue

            # Get frequency data from exome or genome
            exome_data = variant.get("exome", {})
            genome_data = variant.get("genome", {})

            # Prefer exome data if available, otherwise use genome data
            freq_data = exome_data or genome_data

            if not freq_data:  # Skip if no frequency data available
                continue

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

            # Add population-specific data if available
            if "populations" in freq_data:
                for pop in freq_data["populations"]:
                    pop_id = pop["id"]
                    variant_info[f"{pop_id}_ac"] = pop.get("ac")
                    variant_info[f"{pop_id}_an"] = pop.get("an")
                    variant_info[f"{pop_id}_ac_hom"] = pop.get("ac_hom")
                    variant_info[f"{pop_id}_ac_hemi"] = pop.get("ac_hemi")

            processed_variants.append(variant_info)

        return pd.DataFrame(processed_variants)

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

        # Get variant IDs from esearch
        search_params = {
            "db": "clinvar",
            "term": f"{gene_name}[gene]",
            "retmax": 5000,  # Adjust as needed
            "retmode": "json",
            "tool": "swissisoform",
            "email": "mdiberna@wi.mit.edu",  # Replace with actual email
        }

        if self.api_key:
            search_params["api_key"] = self.api_key

        try:
            # First get the variant IDs
            response = await self._async_get_request(
                self.clinvar_endpoints["esearch"], params=search_params
            )
            await asyncio.sleep(0.35)  # Rate limit: max 3 requests per second

            data = json.loads(response)
            if "esearchresult" not in data or "idlist" not in data["esearchresult"]:
                return pd.DataFrame()

            variant_ids = data["esearchresult"]["idlist"]
            if not variant_ids:
                return pd.DataFrame()

            # Process variants in batches
            batch_size = 50
            all_variants = []

            for i in range(0, len(variant_ids), batch_size):
                batch_ids = variant_ids[i : i + batch_size]

                # Get summaries for batch
                summary_params = {
                    "db": "clinvar",
                    "id": ",".join(batch_ids),
                    "retmode": "json",
                    "tool": "swissisoform",
                    "email": "mdiberna@wi.mit.edu",  # Replace with actual email
                }

                if self.api_key:
                    summary_params["api_key"] = self.api_key

                summary_response = await self._async_get_request(
                    self.clinvar_endpoints["esummary"], params=summary_params
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
            logger.error(f"Error fetching ClinVar data: {str(e)}")
            return pd.DataFrame()

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

            # Extract basic variant information
            variant_info = {
                "variant_id": uid,
                "obj_type": variant_data.get("obj_type", ""),
                "accession": variant_data.get("accession", ""),
                "title": variant_data.get("title", ""),
                "gene_name": "",  # Will be populated from genes list
                "review_status": "",  # Will be populated from classifications
                "clinical_significance": "",  # Will be populated from classifications
                "last_evaluated": "",  # Will be populated from classifications
                "molecular_consequences": "",
                "protein_change": variant_data.get("protein_change", ""),
                "chromosome": "",  # Will be populated from variation_set
                "start": "",  # Will be populated from variation_set
                "stop": "",  # Will be populated from variation_set
                "ref_allele": "",  # Will be populated from variation_set
                "alt_allele": "",  # Will be populated from variation_set
                "submission_count": len(
                    variant_data.get("supporting_submissions", {}).get("scv", [])
                ),
            }

            # Get gene information
            genes = variant_data.get("genes", [])
            if genes:
                variant_info["gene_name"] = genes[0].get("symbol", "")
                variant_info["gene_id"] = genes[0].get("geneid", "")
                variant_info["gene_strand"] = genes[0].get("strand", "")

            # Get clinical significance information
            # Check all possible classification types
            for class_type in [
                "germline_classification",
                "clinical_impact_classification",
                "oncogenicity_classification",
            ]:
                if class_type in variant_data:
                    classification = variant_data[class_type]
                    if classification.get(
                        "description"
                    ):  # Only use if there's a description
                        variant_info["clinical_significance"] = classification.get(
                            "description", ""
                        )
                        variant_info["review_status"] = classification.get(
                            "review_status", ""
                        )
                        variant_info["last_evaluated"] = classification.get(
                            "last_evaluated", ""
                        )
                        break  # Use the first classification that has a description

            # Get molecular consequences
            mol_cons = variant_data.get("molecular_consequence_list", [])
            variant_info["molecular_consequences"] = (
                ";".join(mol_cons) if mol_cons else ""
            )

            # Get variant location information from variation_set
            variation_sets = variant_data.get("variation_set", [])
            if variation_sets:
                var_set = variation_sets[0]  # Use first variation set

                # Get variant details
                variant_info["variant_name"] = var_set.get("variation_name", "")
                variant_info["cdna_change"] = var_set.get("cdna_change", "")
                variant_info["canonical_spdi"] = var_set.get("canonical_spdi", "")

                # Get location information
                locations = var_set.get("variation_loc", [])
                for loc in locations:
                    if (
                        loc.get("assembly_name") == "GRCh38"
                        and loc.get("status") == "current"
                    ):
                        variant_info["chromosome"] = loc.get("chr", "")
                        variant_info["start"] = loc.get("start", "")
                        variant_info["stop"] = loc.get("stop", "")
                        variant_info["cytogenic_location"] = loc.get("band", "")
                        break

            variants.append(variant_info)

        df = pd.DataFrame(variants)

        # Convert relevant columns to appropriate types
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

        if variants_df.empty:
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

        summary = {
            "total_variants": len(variants_df),
            # Clinical significance distribution
            "clinical_significance": variants_df["clinical_significance"]
            .value_counts()
            .to_dict(),
            # Variant types
            "variant_types": variants_df["obj_type"].value_counts().to_dict(),
            # Molecular consequences
            "molecular_consequences": variants_df["molecular_consequences"]
            .str.split(";")
            .explode()
            .value_counts()
            .to_dict(),
            # Chromosome distribution
            "chromosome_distribution": variants_df["chromosome"]
            .value_counts()
            .to_dict(),
            # Review status metrics
            "review_status": variants_df["review_status"].value_counts().to_dict(),
            # Additional metadata
            "metadata": {
                "submission_count": {
                    "total": variants_df["submission_count"].sum(),
                    "mean": round(variants_df["submission_count"].mean(), 2),
                    "max": variants_df["submission_count"].max(),
                },
                "unique_accessions": variants_df["accession"].nunique(),
            },
        }

        # Add pathogenicity categorization
        pathogenicity_mapping = {
            "likely_pathogenic_or_pathogenic": variants_df[
                "clinical_significance"
            ].str.contains("pathogenic", case=False, na=False)
            & ~variants_df["clinical_significance"].str.contains(
                "conflicting", case=False, na=False
            ),
            "benign_or_likely_benign": variants_df[
                "clinical_significance"
            ].str.contains("benign", case=False, na=False)
            & ~variants_df["clinical_significance"].str.contains(
                "conflicting", case=False, na=False
            ),
            "uncertain_significance": variants_df["clinical_significance"].str.contains(
                "uncertain|conflicting", case=False, na=False
            ),
            "not_provided": variants_df["clinical_significance"].isna()
            | (variants_df["clinical_significance"] == ""),
        }

        summary["pathogenicity_categories"] = {
            category: int(variants_df[mask].shape[0])
            for category, mask in pathogenicity_mapping.items()
        }

        return summary

    def _find_cosmic_database(self) -> Dict[str, Optional[str]]:
        """Find COSMIC Parquet files and return paths for different datasets"""
        import os
        from pathlib import Path
        
        search_paths = [
            "../data/mutation_data",
            "data/mutation_data", 
            ".",
            "../data/cosmic",
            "data/cosmic",
        ]
        
        cosmic_files = {
            'cancer_census': None,
            'noncoding': None
        }
        
        for search_path in search_paths:
            path = Path(search_path)
            if path.exists():
                # Look for Cancer Mutation Census Parquet files (prefer Parquet over TSV)
                if not cosmic_files['cancer_census']:
                    parquet_files = list(path.glob("*CancerMutationCensus_AllData*.parquet"))
                    if parquet_files:
                        cosmic_files['cancer_census'] = str(parquet_files[0])
                        logger.info(f"Found COSMIC Cancer Census (Parquet): {parquet_files[0]}")
                    else:
                        # Fallback to TSV if no Parquet
                        tsv_files = list(path.glob("*CancerMutationCensus_AllData*.tsv"))
                        if tsv_files:
                            cosmic_files['cancer_census'] = str(tsv_files[0])
                            logger.info(f"Found COSMIC Cancer Census (TSV): {tsv_files[0]}")
                
                # Look for NonCoding variants Parquet files (prefer Parquet over TSV)
                if not cosmic_files['noncoding']:
                    parquet_files = list(path.glob("*Cosmic_NonCodingVariants*.parquet"))
                    if parquet_files:
                        cosmic_files['noncoding'] = str(parquet_files[0])
                        logger.info(f"Found COSMIC NonCoding variants (Parquet): {parquet_files[0]}")
                    else:
                        # Fallback to TSV if no Parquet
                        tsv_files = list(path.glob("*Cosmic_NonCodingVariants*.tsv"))
                        if tsv_files:
                            cosmic_files['noncoding'] = str(tsv_files[0])
                            logger.info(f"Found COSMIC NonCoding variants (TSV): {tsv_files[0]}")
        
        return cosmic_files

    def _fetch_cosmic_sync(self, gene_name: str) -> pd.DataFrame:
        """Synchronous wrapper that combines both COSMIC datasets using fast Parquet/TSV reading"""
        try:
            cosmic_files = self._find_cosmic_database()
            combined_results = []
            
            # Process Cancer Mutation Census
            if cosmic_files['cancer_census']:
                try:
                    logger.info(f"Querying Cancer Mutation Census for {gene_name}")
                    file_path = cosmic_files['cancer_census']
                    
                    if file_path.endswith('.parquet'):
                        # Fast Parquet query
                        import pyarrow.parquet as pq
                        import pyarrow.compute as pc
                        
                        table = pq.read_table(file_path)
                        # Filter by gene name
                        mask = pc.match_substring(
                            pc.utf8_upper(table['GENE_NAME']), 
                            gene_name.upper()
                        )
                        filtered_table = table.filter(mask)
                        cancer_result = filtered_table.to_pandas()
                    else:
                        # TSV fallback with chunked reading for large files
                        cancer_result = self._query_tsv_by_gene(file_path, gene_name, 'GENE_NAME')
                    
                    if cancer_result is not None and not cancer_result.empty:
                        # Add dataset identifier
                        cancer_result['cosmic_dataset'] = 'cancer_census'
                        combined_results.append(cancer_result)
                        logger.info(f"Found {len(cancer_result)} variants in Cancer Census")
                        
                except Exception as e:
                    logger.warning(f"Cancer Mutation Census query failed: {str(e)}")
            
            # Process NonCoding variants
            if cosmic_files['noncoding']:
                try:
                    logger.info(f"Querying NonCoding variants for {gene_name}")
                    file_path = cosmic_files['noncoding']
                    
                    if file_path.endswith('.parquet'):
                        # Fast Parquet query
                        import pyarrow.parquet as pq
                        import pyarrow.compute as pc
                        
                        table = pq.read_table(file_path)
                        # Filter by gene symbol
                        mask = pc.match_substring(
                            pc.utf8_upper(table['GENE_SYMBOL']), 
                            gene_name.upper()
                        )
                        filtered_table = table.filter(mask)
                        noncoding_result = filtered_table.to_pandas()
                    else:
                        # TSV fallback with chunked reading
                        noncoding_result = self._query_tsv_by_gene(file_path, gene_name, 'GENE_SYMBOL')
                    
                    if noncoding_result is not None and not noncoding_result.empty:
                        # Add dataset identifier
                        noncoding_result['cosmic_dataset'] = 'noncoding'
                        combined_results.append(noncoding_result)
                        logger.info(f"Found {len(noncoding_result)} variants in NonCoding dataset")
                        
                except Exception as e:
                    logger.warning(f"NonCoding variants query failed: {str(e)}")
            
            # Combine results
            if combined_results:
                # Concatenate all datasets
                combined_df = pd.concat(combined_results, ignore_index=True, sort=False)
                logger.info(f"Combined {len(combined_df)} total COSMIC variants for {gene_name}")
                
                # Process the combined dataset
                return self._process_cosmic_variants(combined_df)
            else:
                logger.warning(f"No COSMIC data found for {gene_name}")
                return pd.DataFrame()
                
        except Exception as e:
            logger.error(f"COSMIC data fetch failed: {str(e)}")
            return pd.DataFrame()

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

    def _extract_cosmic_position(self, variant: pd.Series) -> Optional[int]:
        """Extract genomic position from COSMIC variant data (both datasets)"""
        try:
            # For NonCoding variants: use GENOME_START
            if 'GENOME_START' in variant and pd.notna(variant.get('GENOME_START')):
                return int(variant.get('GENOME_START'))
            
            # For Cancer Census: try the mutation position fields
            pos_fields = [
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
            'Mutation_genome_position_GRCh38', 
            'Mutation_genome_position_GRCh37', 
        ]
        
        for field in pos_fields:
            if field in variant and pd.notna(variant.get(field)):
                pos_str = str(variant.get(field))
                if ':' in pos_str:
                    return pos_str.split(':')[0]
        
        return ''

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
            # Run _fetch_cosmic_sync in executor since it's not async
            loop = asyncio.get_event_loop()
            cosmic_data = await loop.run_in_executor(
                None, 
                self._fetch_cosmic_sync, 
                gene_name
            )
            
            # Cache the result
            self.cached_data[cache_key] = cosmic_data
            return cosmic_data
            
        except Exception as e:
            logger.error(f"Failed to fetch COSMIC data for {gene_name}: {str(e)}")
            return pd.DataFrame()

    def _standardize_cosmic_impact(self, cosmic_impact: str) -> str:
        """Convert COSMIC mutation classifications to standardized impact terms"""
        if pd.isna(cosmic_impact) or not isinstance(cosmic_impact, str):
            return "unknown"
        
        cosmic_impact = cosmic_impact.lower().strip()
        
        # Map COSMIC descriptions to standard terms
        if 'missense' in cosmic_impact:
            return 'missense variant'
        elif 'nonsense' in cosmic_impact or 'stop' in cosmic_impact:
            return 'nonsense variant'
        elif 'frameshift' in cosmic_impact:
            return 'frameshift variant'
        elif 'synonymous' in cosmic_impact or 'silent' in cosmic_impact:
            return 'synonymous variant'
        elif 'splice' in cosmic_impact:
            return 'splice variant'
        elif 'inframe' in cosmic_impact:
            if 'deletion' in cosmic_impact:
                return 'inframe deletion'
            elif 'insertion' in cosmic_impact:
                return 'inframe insertion'
            else:
                return 'inframe variant'
        elif 'noncoding' in cosmic_impact:
            return 'noncoding variant'
        elif 'deletion' in cosmic_impact:
            return 'deletion'
        elif 'insertion' in cosmic_impact:
            return 'insertion'
        elif 'substitution' in cosmic_impact:
            return 'substitution'
        else:
            return 'other variant'

    async def get_cosmic_summary(self, gene_name: str) -> Dict:
        """Get summary statistics for COSMIC variants using actual data structure"""
        cosmic_df = await self.get_cosmic_variants(gene_name)
        
        if cosmic_df.empty:
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
        
        # Calculate COSMIC-specific frequencies
        total_tested = cosmic_df['cosmic_samples_tested'].apply(
            lambda x: int(x) if pd.notna(x) and str(x).isdigit() else 0
        ).sum()
        
        total_mutated = cosmic_df['cosmic_samples_mutated'].apply(
            lambda x: int(x) if pd.notna(x) and str(x).isdigit() else 0
        ).sum()
        
        summary = {
            'total_variants': len(cosmic_df),
            'mutation_types_cds': cosmic_df['mutation_type_cds'].value_counts().to_dict(),
            'mutation_types_aa': cosmic_df['mutation_type_aa'].value_counts().to_dict(),
            'diseases': cosmic_df['disease'].value_counts().head(10).to_dict(),  # Top 10 diseases
            'oncogene_tsg_roles': cosmic_df['oncogene_tsg'].value_counts().to_dict(),
            'cgc_tiers': cosmic_df['cgc_tier'].value_counts().to_dict(),
            'mutation_significance_tiers': cosmic_df['mutation_significance_tier'].value_counts().to_dict(),
            'clinvar_significance': cosmic_df['clinvar_significance'].value_counts().to_dict(),
            'sift_predictions': cosmic_df['sift_prediction'].value_counts().to_dict(),
            'cosmic_samples_data': {
                'total_tested': total_tested,
                'total_mutated': total_mutated,
                'frequency': total_mutated / total_tested if total_tested > 0 else 0
            },
            'gnomad_frequencies': {
                'with_exome_freq': cosmic_df['gnomad_exomes_af'].notna().sum(),
                'with_genome_freq': cosmic_df['gnomad_genomes_af'].notna().sum()
            }
        }
        
        return summary

    def _standardize_impact_category(self, impact: str) -> str:
        """Standardize impact categories to a consistent set of terms.

        Args:
            impact: Raw impact/consequence string

        Returns:
            Standardized impact category
        """
        if pd.isna(impact) or not isinstance(impact, str):
            return "unknown"

        impact = impact.lower().strip()

        # Check each category
        for category, terms in self.mutation_categories.items():
            if any(term in impact for term in terms):
                if category == "missense":
                    return "missense variant"
                elif category == "synonymous":
                    return "synonymous variant"
                elif category == "nonsense":
                    return "nonsense variant"
                elif category == "inframe_del":
                    return "inframe deletion"
                elif category == "inframe_ins":
                    return "inframe insertion"
                elif category == "frameshift":
                    return "frameshift variant"
                elif category == "splice":
                    return "splice variant"
                elif category == "start_lost":
                    return "start lost variant"
                elif category == "utr_5prime":
                    return "5 prime UTR variant"
                elif category == "utr_3prime":
                    return "3 prime UTR variant"

        return "other variant"

    def standardize_mutation_data(
        self, variants_df: pd.DataFrame, source: str
    ) -> pd.DataFrame:
        """Standardize mutation data from different sources into a consistent format.

        Args:
            variants_df: Raw variants DataFrame from any source
            source: Source of the data ('gnomad', 'clinvar', or 'aggregator')

        Returns:
            Standardized mutation data
        """
        # Define standard columns for visualization
        standard_columns = [
            "position",  # Genomic position
            "variant_id",  # Unique identifier
            "reference",  # Reference allele
            "alternate",  # Alternate allele
            "source",  # Data source
            "impact",  # Functional impact/consequence
            "hgvsc",  # Coding sequence change
            "hgvsp",  # Protein change
            "allele_frequency",  # Frequency in population
            "clinical_significance",  # Clinical interpretation
        ]

        standardized_df = pd.DataFrame(columns=standard_columns)

        if variants_df.empty:
            return standardized_df

        if source.lower() == "gnomad":
            standardized_df = pd.DataFrame(
                {
                    "position": variants_df["position"],
                    "variant_id": variants_df["variant_id"],
                    "reference": variants_df["reference"],
                    "alternate": variants_df["alternate"],
                    "source": "gnomAD",
                    "impact": variants_df["consequence"].apply(
                        self._standardize_impact_category
                    ),
                    "hgvsc": variants_df["hgvsc"],
                    "hgvsp": variants_df["hgvsp"],
                    "allele_frequency": variants_df["allele_frequency"],
                    "clinical_significance": None,
                }
            )

        elif source.lower() == "clinvar":
            standardized_df = pd.DataFrame(
                {
                    "position": variants_df["start"],
                    "variant_id": variants_df["accession"],
                    "reference": variants_df["ref_allele"],
                    "alternate": variants_df["alt_allele"],
                    "source": "ClinVar",
                    "impact": variants_df["molecular_consequences"].apply(
                        self._standardize_impact_category
                    ),
                    "hgvsc": variants_df["cdna_change"],
                    "hgvsp": variants_df["protein_change"],
                    "allele_frequency": None,
                    "clinical_significance": variants_df["clinical_significance"],
                }
            )

        elif source.lower() == "cosmic":
            standardized_df = pd.DataFrame({
                "position": variants_df["position"],
                "variant_id": variants_df["cosmic_id"],
                "reference": variants_df["reference"],  # GENOMIC_WT_ALLELE_SEQ
                "alternate": variants_df["alternate"],  # GENOMIC_MUT_ALLELE_SEQ
                "source": "COSMIC",
                "impact": variants_df["consequence"].apply(self._standardize_cosmic_impact),
                "hgvsc": variants_df["hgvs_coding"],  # Mutation CDS
                "hgvsp": variants_df["hgvs_protein"],  # Mutation AA
                "allele_frequency": pd.to_numeric(variants_df["gnomad_exomes_af"], errors='coerce'),  # Use gnomAD freq from COSMIC
                "clinical_significance": variants_df["clinvar_significance"],  # CLINVAR_CLNSIG from COSMIC
            })

        return self._clean_mutation_data(standardized_df)

    def _clean_mutation_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Clean and validate the standardized mutation data.

        Args:
            df: Standardized mutation DataFrame

        Returns:
            Cleaned mutation data
        """
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

    async def _async_get_request(self, url: str, params: Dict = None) -> str:
        """Make an async GET request.

        Args:
            url: URL to request
            params: Query parameters

        Returns:
            Response text
        """
        # Initialize params if None
        if params is None:
            params = {}

        # Add API key to params if available
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
            aggregator_csv_path: Path to aggregator CSV file
            verbose: Whether to print detailed information

        Returns:
            DataFrame containing mutations ready for visualization
        """
        if sources is None:
            sources = ["gnomad", "clinvar", "cosmic"]  # ADD "cosmic" here
        sources = [s.lower() for s in sources]

        # Fetch and combine data from sources
        dfs_to_combine = []

        if verbose:
            print(f"Fetching mutations from sources: {', '.join(sources)}...")

        for source in sources:
            try:
                if source == "gnomad":
                    df = await self.get_gnomad_variants(gene_name)
                    standardized_df = self.standardize_mutation_data(df, source)
                elif source == "clinvar":
                    df = await self.get_clinvar_variants(gene_name)
                    standardized_df = self.standardize_mutation_data(df, source)
                elif source == "cosmic":  
                    df = await self.get_cosmic_variants(gene_name)
                    standardized_df = self.standardize_mutation_data(df, source)
                else:
                    continue

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

        # Convert position to numeric
        combined_df["position"] = pd.to_numeric(
            combined_df["position"], errors="coerce"
        )

        # Filter by alt_features if provided
        if alt_features is not None and not alt_features.empty:
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
                feature_mask = (combined_df["position"] >= start_pos) & (
                    combined_df["position"] <= end_pos
                )
                feature_mutations = combined_df[feature_mask].copy()

                # Add feature information
                if not feature_mutations.empty:
                    feature_mutations["alt_feature_region"] = (
                        f"{int(start_pos)}-{int(end_pos)}"
                    )
                    feature_mutations["alt_feature_name"] = name
                    filtered_dfs.append(feature_mutations)

                if verbose:
                    print(
                        f"  Feature {name} ({start_pos}-{end_pos}): {len(feature_mutations)} mutations"
                    )

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

        else:
            # No filtering, just add relative positions
            if len(combined_df) > 0:
                min_pos = combined_df["position"].min()
                combined_df["relative_position"] = combined_df["position"] - min_pos

            return combined_df

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
