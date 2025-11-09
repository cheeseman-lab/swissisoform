"""Mutation data handling and analysis module.

This module provides the MutationHandler class for fetching, processing, and
analyzing mutation data from various sources including gnomAD, ClinVar, and COSMIC.
"""

from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport
import pandas as pd
from typing import Dict, Optional, List, Tuple
import requests
import json
import asyncio
import logging
from pathlib import Path
import aiohttp

from .config import Config

logger = logging.getLogger(__name__)


def is_start_codon_loss(ref_codon: str, alt_codon: str) -> bool:
    """Check if a codon change represents a start codon loss.

    Start loss rules:
    1. AUG → anything else = START LOSS
    2. CUG/UUG/GUG/ACG → anything except (CUG/UUG/GUG/ACG or AUG) = START LOSS

    Args:
        ref_codon: Reference codon (RNA notation: AUG, CUG, etc.)
        alt_codon: Alternate codon after mutation (RNA notation)

    Returns:
        True if this is a start codon loss, False otherwise
    """
    if not ref_codon or not alt_codon:
        return False

    ref_codon = ref_codon.upper()
    alt_codon = alt_codon.upper()

    # Define start codons
    canonical_start = "AUG"
    alternative_starts = {"CUG", "UUG", "GUG", "ACG"}
    all_starts = alternative_starts | {canonical_start}

    # Rule 1: AUG → anything else is always a loss
    if ref_codon == canonical_start and alt_codon != canonical_start:
        return True

    # Rule 2: Alternative start → non-start (excluding AUG) is a loss
    if ref_codon in alternative_starts and alt_codon not in all_starts:
        return True

    return False


class MutationHandler:
    """Handles variant data from multiple mutation databases.

    Provides methods to fetch, parse, and analyze mutation data from:
    - gnomAD
    - ClinVar
    - COSMIC
    """

    def __init__(self, api_key: Optional[str] = None, genome_handler=None):
        """Initialize the mutation handler with API endpoints.

        Args:
            api_key (Optional[str]): Optional API key for NCBI E-utilities.
                                    If not provided, loads from Config (environment variables or .env file).
            genome_handler (Optional[GenomeHandler]): Genome handler for sequence extraction (needed for ClinVar deletion inference).
        """
        self._setup_gnomad_client()
        self._setup_clinvar_endpoints(api_key or Config.NCBI_API_KEY)
        self._setup_mutation_categories()
        self.cached_data = {}
        self.genome_handler = genome_handler

    def _setup_gnomad_client(self):
        """Initialize the gnomAD GraphQL client with appropriate logging configuration."""
        # Suppress GQL library logging
        logging.getLogger("gql").setLevel(logging.WARNING)
        logging.getLogger("gql.transport").setLevel(logging.WARNING)
        logging.getLogger("gql.transport.aiohttp").setLevel(logging.WARNING)

        # Configure timeout for reliability (in seconds)
        transport = AIOHTTPTransport(
            url=Config.GNOMAD_API_URL, timeout=Config.GNOMAD_TIMEOUT
        )
        self.gnomad_client = Client(
            transport=transport, fetch_schema_from_transport=True
        )

    def _setup_clinvar_endpoints(self, api_key: Optional[str]):
        """Set up ClinVar API endpoint URLs and store the API key.

        Args:
            api_key (Optional[str]): API key for NCBI E-utilities.
        """
        self.clinvar_base = Config.CLINVAR_BASE_URL
        self.clinvar_endpoints = {
            "esearch": f"{self.clinvar_base}/esearch.fcgi",
            "esummary": f"{self.clinvar_base}/esummary.fcgi",
            "efetch": f"{self.clinvar_base}/efetch.fcgi",
            "elink": f"{self.clinvar_base}/elink.fcgi",
        }
        self.api_key = api_key

    def _setup_mutation_categories(self):
        """Initialize mutation category mappings for variant classification."""
        self.mutation_categories = {
            # UTR impact first
            "utr_5prime": [
                "5_prime_utr",
                "5 prime utr",
                "5_prime_utr_variant",
                "5 prime utr variant",
            ],
            "utr_3prime": [
                "3_prime_utr",
                "3 prime utr",
                "3_prime_utr_variant",
                "3 prime utr variant",
            ],
            # High impact second
            "nonsense": ["nonsense", "stop_gained", "stop gained", "nonsense_variant"],
            "frameshift": ["frameshift", "frameshift_variant", "frameshift variant"],
            "start_lost": ["start_lost", "start lost", "start_loss"],
            "stop_lost": ["stop_lost", "stop lost"],
            # Moderate impact
            "missense": ["missense", "missense_variant", "missense variant"],
            "inframe_del": [
                "inframe_deletion",
                "inframe deletion",
                "inframe_del",
                "inframe_indel",
            ],
            "inframe_ins": ["inframe_insertion", "inframe insertion", "inframe_ins"],
            "synonymous": ["synonymous", "synonymous_variant", "synonymous variant"],
            # Unclassified - requires validation to determine true impact
            "unclassified": [
                "unclassified",
                "unclassified variant",
                "unclassified_variant",
            ],
            # Low impact last
            "splice": [
                "splice",
                "splice_variant",
                "splice_region_variant",
                "splice_donor_variant",
                "splice_acceptor_variant",
            ],
            "intronic": ["intron_variant", "intron variant"],
        }

        # Define which categories to exclude by default
        # Note: intronic variants are now handled via bulk filtering before validation
        # Splice variants are still filtered here as they have VEP annotations
        self.low_impact_categories = ["splice"]

    async def get_gnomad_variants(self, gene_name: str) -> pd.DataFrame:
        """Get processed variant data from gnomAD.

        Args:
            gene_name (str): Gene symbol to query.

        Returns:
            pd.DataFrame: DataFrame of gnomAD variants with standardized columns.
        """
        gnomad_data = await self.fetch_gnomad_data(gene_name)
        return self.process_gnomad_variants(gnomad_data)

    async def get_gnomad_summary(self, gene_name: str) -> Dict:
        """Get summary statistics for gnomAD variants.

        Args:
            gene_name (str): Gene symbol to query.

        Returns:
            Dict: Dictionary of summary statistics including variant counts and frequencies.
        """
        gnomad_df = await self.get_gnomad_variants(gene_name)

        # Log summary
        if not gnomad_df.empty:
            logger.info(f"Found {len(gnomad_df)} gnomAD variants for {gene_name}")
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
        """Return empty gnomAD summary structure.

        Returns:
            Dict: Empty summary dictionary with default values.
        """
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

        Implements retry logic with exponential backoff for improved reliability.

        Args:
            gene_name (str): Gene symbol (e.g., 'BRCA1').
            reference_genome (str): Reference genome version. Defaults to 'GRCh38'.

        Returns:
            Dict: Processed gnomAD data for the gene. Returns empty dict on failure.
        """
        cache_key = f"gnomad_{gene_name}_{reference_genome}"
        if cache_key in self.cached_data:
            return self.cached_data[cache_key]

        query = self._get_gnomad_query()
        variables = {"geneSymbol": gene_name, "referenceGenome": reference_genome}

        # Retry logic with exponential backoff
        for attempt in range(Config.GNOMAD_MAX_RETRIES):
            try:
                result = await self.gnomad_client.execute_async(
                    query, variable_values=variables
                )

                if result.get("gene") is not None:
                    self.cached_data[cache_key] = result
                    if attempt > 0:
                        logger.info(
                            f"gnomAD request succeeded for {gene_name} on attempt {attempt + 1}"
                        )
                    return result
                else:
                    # Gene not found in gnomAD (not an error, just no data)
                    logger.debug(f"Gene {gene_name} not found in gnomAD")
                    return {}

            except asyncio.TimeoutError:
                wait_time = Config.GNOMAD_RETRY_DELAY * (2**attempt)
                if attempt < Config.GNOMAD_MAX_RETRIES - 1:
                    logger.warning(
                        f"gnomAD request timeout for {gene_name}, retrying in {wait_time}s (attempt {attempt + 1}/{Config.GNOMAD_MAX_RETRIES})"
                    )
                    await asyncio.sleep(wait_time)
                else:
                    logger.error(
                        f"gnomAD request timeout for {gene_name} after {Config.GNOMAD_MAX_RETRIES} attempts"
                    )
                    return {}

            except (aiohttp.ClientError, ConnectionError) as e:
                wait_time = Config.GNOMAD_RETRY_DELAY * (2**attempt)
                if attempt < Config.GNOMAD_MAX_RETRIES - 1:
                    logger.warning(
                        f"gnomAD connection error for {gene_name}: {str(e)}, retrying in {wait_time}s (attempt {attempt + 1}/{Config.GNOMAD_MAX_RETRIES})"
                    )
                    await asyncio.sleep(wait_time)
                else:
                    logger.error(
                        f"gnomAD connection failed for {gene_name} after {Config.GNOMAD_MAX_RETRIES} attempts: {str(e)}"
                    )
                    return {}

            except Exception as e:
                # Unexpected errors - log and return empty to avoid crashing pipeline
                logger.error(
                    f"Unexpected error fetching gnomAD data for {gene_name}: {str(e)}"
                )
                return {}

        return {}

    def _get_gnomad_query(self):
        """Get the GraphQL query for gnomAD data.

        Returns:
            gql.Document: GraphQL query document for fetching gnomAD variant data.
        """
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
            data (Dict): Raw gnomAD API response.
            filter_field (Optional[str]): Field to filter on.
            filter_value (Optional[str]): Value to filter for.

        Returns:
            pd.DataFrame: DataFrame of processed variant data.
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
        """Validate gnomAD data structure.

        Args:
            data (Dict): Raw gnomAD API response data.

        Returns:
            bool: True if data structure is valid, False otherwise.
        """
        return data and "gene" in data and data["gene"] and "variants" in data["gene"]

    def _should_skip_variant(
        self, variant: Dict, filter_field: str, filter_value: str
    ) -> bool:
        """Check if variant should be skipped based on filters.

        Args:
            variant (Dict): Variant data dictionary.
            filter_field (str): Field name to filter on.
            filter_value (str): Value to filter for.

        Returns:
            bool: True if variant should be skipped, False otherwise.
        """
        if filter_field and filter_value:
            variant_field_value = variant.get(filter_field)
            if variant_field_value != filter_value:
                return True
        return False

    def _process_gnomad_variant(self, variant: Dict) -> Optional[Dict]:
        """Process a single gnomAD variant.

        Args:
            variant (Dict): Raw variant data from gnomAD API.

        Returns:
            Optional[Dict]: Processed variant information, or None if no frequency data available.
        """
        exome_data = variant.get("exome", {})
        genome_data = variant.get("genome", {})
        freq_data = exome_data or genome_data

        if not freq_data:
            return None

        # Extract chromosome from variant_id (format: "chr-pos-ref-alt")
        variant_id = variant["variant_id"]
        chromosome = variant_id.split("-")[0] if "-" in variant_id else ""

        variant_info = {
            "position": variant["pos"],
            "variant_id": variant_id,
            "chromosome": chromosome,
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
        """Add population-specific frequency data to variant info.

        Args:
            variant_info (Dict): Variant information dictionary to update.
            freq_data (Dict): Frequency data containing population-specific information.
        """
        if "populations" in freq_data:
            for pop in freq_data["populations"]:
                pop_id = pop["id"]
                variant_info[f"{pop_id}_ac"] = pop.get("ac")
                variant_info[f"{pop_id}_an"] = pop.get("an")
                variant_info[f"{pop_id}_ac_hom"] = pop.get("ac_hom")
                variant_info[f"{pop_id}_ac_hemi"] = pop.get("ac_hemi")

    async def get_clinvar_variants(
        self, gene_name: str, refseq_id: str = None
    ) -> pd.DataFrame:
        """Get all ClinVar variants for a gene.

        Args:
            gene_name (str): Gene symbol (e.g., 'BRCA1').
            refseq_id (Optional[str]): RefSeq transcript ID to filter ClinVar variants (transcript-level query).

        Returns:
            pd.DataFrame: DataFrame of processed ClinVar variant data.
        """
        cache_key = f"clinvar_{gene_name}_{refseq_id or 'all'}"
        if cache_key in self.cached_data:
            return self.cached_data[cache_key]

        try:
            variant_ids = await self._fetch_clinvar_variant_ids(gene_name, refseq_id)
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

    async def _fetch_clinvar_variant_ids(
        self, gene_name: str, refseq_id: str = None
    ) -> List[str]:
        """Fetch ClinVar variant IDs for a gene, optionally filtered by RefSeq transcript."""
        if refseq_id and refseq_id != "NA":
            search_term = f"{gene_name}[gene] AND {refseq_id}[transcript]"
        else:
            search_term = f"{gene_name}[gene]"

        search_params = {
            "db": "clinvar",
            "term": search_term,
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

    async def _process_clinvar_batches(
        self, variant_ids: List[str]
    ) -> List[pd.DataFrame]:
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

    def _extract_clinvar_clinical_significance(
        self, variant_info: Dict, variant_data: Dict
    ):
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
                    variant_info["clinical_significance"] = classification.get(
                        "description", ""
                    )
                    variant_info["review_status"] = classification.get(
                        "review_status", ""
                    )
                    variant_info["last_evaluated"] = classification.get(
                        "last_evaluated", ""
                    )
                    break

    def _extract_clinvar_molecular_consequences(
        self, variant_info: Dict, variant_data: Dict
    ):
        """Extract molecular consequences from ClinVar variant.

        Takes only the first (primary) consequence when multiple are present.
        """
        mol_cons = variant_data.get("molecular_consequence_list", [])

        if mol_cons:
            # Store ONLY the primary (first) consequence
            variant_info["molecular_consequences"] = mol_cons[0].strip()
        else:
            variant_info["molecular_consequences"] = ""

    def _extract_clinvar_location_info(self, variant_info: Dict, variant_data: Dict):
        """Extract location information from ClinVar variant."""
        variation_sets = variant_data.get("variation_set", [])
        if not variation_sets:
            return

        var_set = variation_sets[0]
        variant_info["variant_name"] = var_set.get("variation_name", "")
        variant_info["cdna_change"] = var_set.get("cdna_change", "")
        variant_info["canonical_spdi"] = var_set.get("canonical_spdi", "")

        # Extract genomic alleles from canonical_spdi
        canonical_spdi = var_set.get("canonical_spdi", "")
        if canonical_spdi:
            spdi_parts = canonical_spdi.split(":")
            if len(spdi_parts) == 4:
                # SPDI format: sequence:position:deletion:insertion
                variant_info["ref_allele"] = spdi_parts[
                    2
                ]  # deletion (reference allele)
                variant_info["alt_allele"] = spdi_parts[
                    3
                ]  # insertion (alternate allele)

        # Extract location from GRCh38 (prioritize current, fallback to any GRCh38, then any location)
        locations = var_set.get("variation_loc", [])
        location_found = False

        # Try 1: GRCh38 with current status
        for loc in locations:
            if loc.get("assembly_name") == "GRCh38" and loc.get("status") == "current":
                variant_info["chromosome"] = loc.get("chr", "")
                variant_info["start"] = loc.get("start", "")
                variant_info["stop"] = loc.get("stop", "")
                variant_info["cytogenic_location"] = loc.get("band", "")
                location_found = True
                break

        # Try 2: Any GRCh38 location
        if not location_found:
            for loc in locations:
                if loc.get("assembly_name") == "GRCh38":
                    variant_info["chromosome"] = loc.get("chr", "")
                    variant_info["start"] = loc.get("start", "")
                    variant_info["stop"] = loc.get("stop", "")
                    variant_info["cytogenic_location"] = loc.get("band", "")
                    location_found = True
                    break

        # Try 3: Any location with chromosome info
        if not location_found:
            for loc in locations:
                if loc.get("chr"):
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
            gene_name (str): Gene symbol.

        Returns:
            Dict: Dictionary of summary statistics.
        """
        variants_df = await self.get_clinvar_variants(gene_name)

        # Log summary
        if not variants_df.empty:
            logger.info(f"Found {len(variants_df)} ClinVar variants for {gene_name}")
        else:
            logger.info(f"No ClinVar variants found for {gene_name}")

        if variants_df.empty:
            return self._get_empty_clinvar_summary()

        summary = {
            "total_variants": len(variants_df),
            "clinical_significance": variants_df["clinical_significance"]
            .value_counts()
            .to_dict(),
            "variant_types": variants_df["obj_type"].value_counts().to_dict(),
            "molecular_consequences": self._get_molecular_consequences_counts(
                variants_df
            ),
            "chromosome_distribution": variants_df["chromosome"]
            .value_counts()
            .to_dict(),
            "review_status": variants_df["review_status"].value_counts().to_dict(),
            "metadata": self._get_clinvar_metadata(variants_df),
        }

        # Add pathogenicity categorization
        summary["pathogenicity_categories"] = self._get_pathogenicity_categories(
            variants_df
        )
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
        return (
            variants_df["molecular_consequences"]
            .str.split(";")
            .explode()
            .value_counts()
            .to_dict()
        )

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
                variants_df["clinical_significance"].str.contains(
                    "pathogenic", case=False, na=False
                )
                & ~variants_df["clinical_significance"].str.contains(
                    "conflicting", case=False, na=False
                )
            ),
            "benign_or_likely_benign": (
                variants_df["clinical_significance"].str.contains(
                    "benign", case=False, na=False
                )
                & ~variants_df["clinical_significance"].str.contains(
                    "conflicting", case=False, na=False
                )
            ),
            "uncertain_significance": variants_df["clinical_significance"].str.contains(
                "uncertain|conflicting", case=False, na=False
            ),
            "not_provided": (
                variants_df["clinical_significance"].isna()
                | (variants_df["clinical_significance"] == "")
            ),
        }

        return {
            category: int(variants_df[mask].shape[0])
            for category, mask in pathogenicity_mapping.items()
        }

    async def get_cosmic_variants(self, gene_name: str) -> pd.DataFrame:
        """Get COSMIC variants for a gene using the downloaded Parquet database.

        Args:
            gene_name (str): Gene symbol (e.g., 'BRCA1').

        Returns:
            pd.DataFrame: DataFrame of processed COSMIC variants.
        """
        cache_key = f"cosmic_{gene_name}"
        if cache_key in self.cached_data:
            return self.cached_data[cache_key]

        try:
            loop = asyncio.get_event_loop()
            cosmic_data = await loop.run_in_executor(
                None, self._fetch_cosmic_sync, gene_name
            )

            self.cached_data[cache_key] = cosmic_data
            return cosmic_data

        except Exception as e:
            logger.error(f"Failed to fetch COSMIC data for {gene_name}: {str(e)}")
            return pd.DataFrame()

    def _find_cosmic_database(self) -> Optional[str]:
        """Find COSMIC variants database file and return path."""
        search_paths = [
            "../data/mutation_data",
            "data/mutation_data",
            ".",
        ]

        for search_path in search_paths:
            path = Path(search_path)
            if not path.exists():
                continue

            # Look for cosmic_variants_combined.parquet
            combined_file = path / "cosmic_variants_combined.parquet"
            if combined_file.exists():
                logger.debug(f"Found COSMIC variants database: {combined_file}")
                return str(combined_file)

        return None

    def _fetch_cosmic_sync(self, gene_name: str) -> pd.DataFrame:
        """Fetch COSMIC data for a gene."""
        try:
            cosmic_file = self._find_cosmic_database()
            if not cosmic_file:
                logger.warning("No COSMIC database file found")
                return pd.DataFrame()

            # Process the single dataset
            result = self._process_cosmic_dataset(cosmic_file, gene_name)

            if result is not None and not result.empty:
                return self._process_cosmic_variants(result)
            else:
                logger.warning(f"No COSMIC data found for {gene_name}")
                return pd.DataFrame()

        except Exception as e:
            logger.error(f"COSMIC data fetch failed: {str(e)}")
            return pd.DataFrame()

    def _process_cosmic_dataset(self, file_path: str, gene_name: str) -> pd.DataFrame:
        """Process the COSMIC variants dataset."""
        try:
            logger.debug(f"Querying COSMIC database for {gene_name}")

            result = self._query_parquet_by_gene(file_path, gene_name, "GENE_SYMBOL")

            if result is not None and not result.empty:
                return result

        except Exception as e:
            logger.warning(f"COSMIC database query failed: {str(e)}")

        return pd.DataFrame()

    def _process_cosmic_variants(self, cosmic_df: pd.DataFrame) -> pd.DataFrame:
        """Process raw COSMIC data into standardized format."""
        if cosmic_df.empty:
            return pd.DataFrame()

        processed_variants = []

        for _, variant in cosmic_df.iterrows():
            # Extract position from GENOME_START
            position = variant.get("GENOME_START")
            if pd.isna(position):
                continue

            variant_info = {
                "position": int(position),
                "cosmic_id": variant.get("GENOMIC_MUTATION_ID", ""),
                "gene_name": variant.get("GENE_SYMBOL", ""),
                "transcript_id": variant.get("TRANSCRIPT_ACCESSION", ""),
                "chromosome": variant.get("CHROMOSOME", ""),
                "reference": variant.get("GENOMIC_WT_ALLELE", ""),
                "alternate": variant.get("GENOMIC_MUT_ALLELE", ""),
                # HGVS notations - these are already properly formatted
                "hgvs_genomic": variant.get("HGVSG", ""),
                "hgvs_coding": variant.get("HGVSC", ""),
                "hgvs_protein": variant.get("HGVSP", ""),
                # Mutation descriptions
                "mutation_cds": variant.get("MUTATION_CDS", ""),
                "mutation_aa": variant.get("MUTATION_AA", ""),
                "consequence": variant.get("MUTATION_DESCRIPTION", ""),
                # Sample and study information
                "sample_name": variant.get("SAMPLE_NAME", ""),
                "cosmic_sample_id": variant.get("COSMIC_SAMPLE_ID", ""),
                "cosmic_study_id": variant.get("COSMIC_STUDY_ID", ""),
                "sample_count": variant.get("GENOME_SCREEN_SAMPLE_COUNT", 0),
                "pubmed_pmid": variant.get("PUBMED_PMID", ""),
                "somatic_status": variant.get("MUTATION_SOMATIC_STATUS", ""),
                "zygosity": variant.get("MUTATION_ZYGOSITY", ""),
            }

            processed_variants.append(variant_info)

        result_df = pd.DataFrame(processed_variants)

        # Deduplicate based on unique genomic variant (chromosome, position, ref, alt)
        # Keep first occurrence to maintain consistency
        if not result_df.empty:
            original_count = len(result_df)
            result_df = result_df.drop_duplicates(
                subset=["chromosome", "position", "reference", "alternate"],
                keep="first",
            )
            deduplicated_count = len(result_df)
            if original_count > deduplicated_count:
                logger.debug(
                    f"Deduplicated COSMIC variants: {original_count} → {deduplicated_count} "
                    f"({original_count - deduplicated_count} duplicates removed)"
                )

        return result_df

    def _query_parquet_by_gene(
        self, file_path: str, gene_name: str, gene_column: str
    ) -> pd.DataFrame:
        """Query Parquet file by gene name using efficient filtering.

        Uses PyArrow filters to read only relevant rows, avoiding loading entire file.
        """
        import pyarrow.parquet as pq
        import pyarrow.dataset as ds

        # Use PyArrow dataset with filter pushdown for efficient querying
        # This only reads rows matching the filter, not the entire file
        dataset = ds.dataset(file_path, format="parquet")

        # Filter for exact gene name match (case-insensitive)
        # Note: PyArrow filters are pushed down to the Parquet reader for efficiency
        table = dataset.to_table(
            filter=(ds.field(gene_column) == gene_name)
            | (ds.field(gene_column) == gene_name.upper())
            | (ds.field(gene_column) == gene_name.lower())
        )

        return table.to_pandas()

    async def get_cosmic_summary(self, gene_name: str) -> Dict:
        """Get summary statistics for COSMIC variants.

        Args:
            gene_name (str): Gene symbol.

        Returns:
            Dict: Dictionary of summary statistics.
        """
        cosmic_df = await self.get_cosmic_variants(gene_name)

        if cosmic_df is None or cosmic_df.empty:
            return self._get_empty_cosmic_summary()

        summary = {
            "total_variants": len(cosmic_df),
            "mutation_descriptions": cosmic_df["consequence"]
            .value_counts(dropna=False)
            .to_dict(),
            "mutation_cds": cosmic_df["mutation_cds"]
            .value_counts(dropna=False)
            .head(10)
            .to_dict(),
            "mutation_aa": cosmic_df["mutation_aa"]
            .value_counts(dropna=False)
            .head(10)
            .to_dict(),
            "somatic_status": cosmic_df["somatic_status"]
            .value_counts(dropna=False)
            .to_dict(),
            "zygosity": cosmic_df["zygosity"].value_counts(dropna=False).to_dict(),
            "sample_count": cosmic_df["cosmic_sample_id"].nunique(),
            "study_count": cosmic_df["cosmic_study_id"].nunique(),
        }
        return summary

    def _get_empty_cosmic_summary(self) -> Dict:
        """Return empty COSMIC summary structure."""
        return {
            "total_variants": 0,
            "mutation_descriptions": {},
            "mutation_cds": {},
            "mutation_aa": {},
            "somatic_status": {},
            "zygosity": {},
            "sample_count": 0,
            "study_count": 0,
        }

    async def get_custom_variants(
        self, parquet_path: str, gene_name: str = None
    ) -> pd.DataFrame:
        """Get variants from a custom user-provided parquet file.

        Args:
            parquet_path (str): Path to custom mutation parquet file.
            gene_name (str, optional): Optional gene name to filter by.

        Returns:
            pd.DataFrame: DataFrame of custom variants (already in standardized format).
        """
        cache_key = f"custom_{parquet_path}_{gene_name}"
        if cache_key in self.cached_data:
            return self.cached_data[cache_key]

        try:
            # Load the parquet file
            df = pd.read_parquet(parquet_path)

            # Validate that it has the required columns
            required_columns = [
                "position",
                "variant_id",
                "reference",
                "alternate",
                "source",
                "impact",
                "hgvsc",
                "hgvsp",
                "allele_frequency",
                "allele_count",
                "allele_count_hom",
                "clinical_significance",
            ]

            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                logger.warning(
                    f"Custom parquet file missing required columns: {missing_columns}"
                )
                # Add missing columns with default values
                for col in missing_columns:
                    if col in ["allele_frequency", "allele_count", "allele_count_hom"]:
                        df[col] = None
                    elif col == "source":
                        df[col] = "Custom"
                    else:
                        df[col] = ""

            # Filter by gene if specified
            if gene_name and "gene_name" in df.columns:
                df = df[df["gene_name"].str.upper() == gene_name.upper()]

            # Ensure position is numeric
            df["position"] = pd.to_numeric(df["position"], errors="coerce")

            # Cache and return
            self.cached_data[cache_key] = df
            logger.info(f"Loaded {len(df)} custom variants from {parquet_path}")
            return df

        except Exception as e:
            logger.error(f"Error loading custom parquet file {parquet_path}: {str(e)}")
            return pd.DataFrame()

    def _standardize_hgvsc(self, hgvsc_str: str) -> str:
        """Standardize HGVS coding sequence notation.

        Args:
            hgvsc_str: Raw HGVS coding sequence string

        Returns:
            Standardized HGVS coding sequence string
        """
        if pd.isna(hgvsc_str) or not isinstance(hgvsc_str, str):
            return ""

        hgvsc_str = hgvsc_str.strip()
        if not hgvsc_str:
            return ""

        # Handle ClinVar format: "NM_005228.5(EGFR):c.-216G>T" -> "c.-216G>T"
        if ":" in hgvsc_str:
            parts = hgvsc_str.split(":")
            if len(parts) >= 2:
                # Take the part after the colon (the actual HGVS description)
                hgvsc_str = parts[-1].strip()

        return hgvsc_str

    def _extract_ref_alt_from_hgvs(self, hgvs_str: str) -> tuple:
        """Extract reference and alternate alleles from HGVS notation.

        Args:
            hgvs_str: HGVS string (coding sequence)

        Returns:
            Tuple of (reference, alternate) or ("", "") if cannot parse
        """
        import re

        if pd.isna(hgvs_str) or not isinstance(hgvs_str, str):
            return "", ""

        hgvs_str = hgvs_str.strip()

        # Clean up ClinVar format first
        if ":" in hgvs_str:
            parts = hgvs_str.split(":")
            if len(parts) >= 2:
                hgvs_str = parts[-1].strip()

        # Handle substitutions: c.14G>T, c.*65A>G, c.-75A>G
        subst_match = re.search(r"([ATCG])>([ATCG])", hgvs_str)
        if subst_match:
            ref, alt = subst_match.groups()
            return ref, alt

        # Handle deletions: c.*71del, c.-72del
        del_match = re.search(r"([ATCG]+)?del", hgvs_str)
        if del_match:
            ref = del_match.group(1) if del_match.group(1) else "N"
            return ref, ""

        # Handle duplications: c.*294dup
        dup_match = re.search(r"([ATCG]+)?dup", hgvs_str)
        if dup_match:
            dup_seq = dup_match.group(1) if dup_match.group(1) else "N"
            return "", dup_seq

        # Handle insertions: c.123_124insA
        ins_match = re.search(r"ins([ATCG]+)", hgvs_str)
        if ins_match:
            ins_seq = ins_match.group(1)
            return "", ins_seq

        return "", ""

    def _standardize_hgvsp(self, hgvsp_str: str, impact: str = "") -> str:
        """Standardize HGVS protein notation.

        For synonymous variants, returns empty string since there's no protein change.

        Args:
            hgvsp_str: Raw HGVS protein string
            impact: Variant impact type to check for synonymous variants

        Returns:
            Standardized HGVS protein string, or empty string for synonymous variants
        """
        if pd.isna(hgvsp_str) or not isinstance(hgvsp_str, str):
            return ""

        # If this is a synonymous variant, there's no protein change
        if impact and "synonymous" in impact.lower():
            return ""

        hgvsp_str = hgvsp_str.strip()
        if not hgvsp_str:
            return ""

        # Handle ClinVar format: remove transcript prefix
        if ":" in hgvsp_str:
            parts = hgvsp_str.split(":")
            if len(parts) >= 2:
                hgvsp_str = parts[-1].strip()

        # Convert protein HGVS to short format if it starts with 'p.'
        if hgvsp_str.startswith("p."):
            return self._convert_protein_to_short_format(hgvsp_str)

        return hgvsp_str

    def _convert_protein_to_short_format(self, protein_hgvs: str) -> str:
        """Convert protein HGVS descriptions to short format.

        Converts 'p.Pro3His' to 'P3H' format.

        Args:
            protein_hgvs: Protein HGVS string starting with 'p.'

        Returns:
            Short format protein change description
        """
        if not protein_hgvs.startswith("p."):
            return protein_hgvs

        # Remove 'p.' prefix
        change = protein_hgvs[2:].strip()

        # Handle complex changes (frameshift, etc.) - return as-is
        if any(
            term in change.lower() for term in ["fs", "ter", "ext", "del", "dup", "ins"]
        ):
            return change

        # Convert to short format
        return self._convert_to_short_aa_format(change)

    def _convert_to_short_aa_format(self, change: str) -> str:
        """Convert amino acid change to short format."""
        import re

        # Amino acid three-letter to one-letter mapping
        aa_map = {
            "Ala": "A",
            "Arg": "R",
            "Asn": "N",
            "Asp": "D",
            "Cys": "C",
            "Gln": "Q",
            "Glu": "E",
            "Gly": "G",
            "His": "H",
            "Ile": "I",
            "Leu": "L",
            "Lys": "K",
            "Met": "M",
            "Phe": "F",
            "Pro": "P",
            "Ser": "S",
            "Thr": "T",
            "Trp": "W",
            "Tyr": "Y",
            "Val": "V",
            "Ter": "*",
            "Stop": "*",
        }

        # Try three-letter format: Pro3His
        match = re.match(r"([A-Za-z]{3})(\d+)([A-Za-z]{3})", change)
        if match:
            aa1, pos, aa2 = match.groups()
            aa1_short = aa_map.get(aa1, aa1)
            aa2_short = aa_map.get(aa2, aa2)
            return f"{aa1_short}{pos}{aa2_short}"

        # If already in single-letter format or doesn't match, return as-is
        return change

    def standardize_mutation_data(
        self, variants_df: pd.DataFrame, source: str, gene_name: str = None
    ) -> pd.DataFrame:
        """Standardize mutation data from different sources into a consistent format.

        Args:
            variants_df (pd.DataFrame): Raw variants DataFrame from any source.
            source (str): Source of the data ('gnomad', 'clinvar', or 'cosmic').
            gene_name (str, optional): Gene name to add to variants missing this field.

        Returns:
            pd.DataFrame: Standardized mutation data.
        """
        if variants_df.empty:
            return self._get_empty_standardized_df()

        source_lower = source.lower()

        if source_lower == "gnomad":
            standardized_df = self._standardize_gnomad_data(variants_df, gene_name)
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
            "position",
            "variant_id",
            "reference",
            "alternate",
            "source",
            "impact",
            "hgvsc",
            "hgvsp",
            "allele_frequency",
            "allele_count",
            "allele_count_hom",
            "sample_count",
            "clinical_significance",
        ]
        return pd.DataFrame(columns=standard_columns)

    def _standardize_gnomad_data(
        self, variants_df: pd.DataFrame, gene_name: str = None
    ) -> pd.DataFrame:
        """Standardize gnomAD data."""
        return pd.DataFrame(
            {
                "position": variants_df["position"],
                "variant_id": variants_df["variant_id"],
                "reference": variants_df["reference"],
                "alternate": variants_df["alternate"],
                "source": "gnomAD",
                "impact": variants_df["consequence"].apply(
                    self._standardize_impact_category
                ),
                "hgvsc": variants_df["hgvsc"].apply(self._standardize_hgvsc),
                "hgvsp": variants_df["hgvsp"].apply(self._standardize_hgvsp),
                "allele_frequency": variants_df["allele_frequency"],
                "allele_count": variants_df["allele_count"],
                "allele_count_hom": variants_df["allele_count_hom"],
                "sample_count": None,
                "clinical_significance": None,
                "chromosome": variants_df["chromosome"]
                if "chromosome" in variants_df.columns
                else "",
                "gene_name": gene_name if gene_name else "",
            }
        )

    def _infer_clinvar_alleles(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Infer missing reference and alternate alleles for ClinVar variants.

        When ClinVar doesn't provide SPDI format (canonical_spdi is empty),
        we need to infer the alleles from the genomic coordinates and the genome sequence.
        This commonly happens for deletions and insertions.

        Args:
            variants_df: DataFrame with ClinVar variant data

        Returns:
            DataFrame with inferred ref_allele and alt_allele filled in
        """
        if self.genome_handler is None:
            logger.debug("No genome handler available for allele inference")
            return variants_df

        variants_df = variants_df.copy()
        inferred_count = 0

        for idx, row in variants_df.iterrows():
            # Skip if alleles are already present
            ref_allele = str(row.get("ref_allele", ""))
            alt_allele = str(row.get("alt_allele", ""))

            if ref_allele and alt_allele:
                continue

            # Need chromosome, start, and stop to infer alleles
            chrom = row.get("chromosome", "")
            start = row.get("start", "")
            stop = row.get("stop", "")

            if not chrom or not start or not stop:
                continue

            try:
                start = int(start)
                stop = int(stop)
            except (ValueError, TypeError):
                continue

            # Infer variant type based on coordinates
            if start == stop:
                # Single nucleotide variant - extract 1bp reference
                try:
                    ref_seq = str(self.genome_handler.get_sequence(chrom, start, start))
                    if ref_seq and not alt_allele:
                        # This is likely an SNV but we don't know the alt
                        # Keep the reference, alt will stay empty
                        variants_df.at[idx, "ref_allele"] = ref_seq
                except Exception as e:
                    logger.debug(f"Could not extract sequence for {chrom}:{start}: {e}")

            elif stop > start:
                # Deletion or complex variant
                # For deletions in VCF format, we need the base BEFORE the deletion
                # plus the deleted sequence as reference, and just the base before as alternate
                try:
                    # Get the base before the deletion
                    anchor_base = str(
                        self.genome_handler.get_sequence(chrom, start - 1, start - 1)
                    )
                    # Get the deleted sequence
                    deleted_seq = str(
                        self.genome_handler.get_sequence(chrom, start, stop)
                    )

                    if anchor_base and deleted_seq:
                        # VCF-style deletion format:
                        # Reference = anchor base + deleted sequence
                        # Alternate = anchor base only
                        ref_vcf = anchor_base + deleted_seq
                        alt_vcf = anchor_base

                        variants_df.at[idx, "ref_allele"] = ref_vcf
                        variants_df.at[idx, "alt_allele"] = alt_vcf
                        # Adjust position to the anchor base
                        variants_df.at[idx, "start"] = start - 1
                        inferred_count += 1
                except Exception as e:
                    logger.debug(
                        f"Could not infer deletion alleles for {row.get('accession', 'unknown')}: {e}"
                    )

        # Log summary if any alleles were inferred
        if inferred_count > 0:
            logger.debug(f"Inferred alleles for {inferred_count} ClinVar variant(s)")

        return variants_df

    def _standardize_clinvar_data(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Standardize ClinVar data."""
        variants_df = variants_df.copy()

        # Create impact column first
        variants_df["standardized_impact"] = variants_df[
            "molecular_consequences"
        ].apply(self._standardize_impact_category)

        # Apply HGVS standardization
        variants_df["standardized_hgvsc"] = variants_df["cdna_change"].apply(
            self._standardize_hgvsc
        )
        variants_df["standardized_hgvsp"] = variants_df.apply(
            lambda row: self._standardize_hgvsp(
                row["protein_change"], row["standardized_impact"]
            ),
            axis=1,
        )

        # Fix missing ref/alt alleles for deletions and insertions
        # This happens when ClinVar doesn't provide SPDI format
        variants_df = self._infer_clinvar_alleles(variants_df)

        return pd.DataFrame(
            {
                "position": variants_df["start"],
                "variant_id": variants_df["accession"],
                "reference": variants_df["ref_allele"],
                "alternate": variants_df["alt_allele"],
                "source": "ClinVar",
                "impact": variants_df["standardized_impact"],
                "hgvsc": variants_df["standardized_hgvsc"],
                "hgvsp": variants_df["standardized_hgvsp"],
                "allele_frequency": None,
                "allele_count": None,
                "allele_count_hom": None,
                "sample_count": None,
                "clinical_significance": variants_df["clinical_significance"],
                "chromosome": variants_df["chromosome"],
                "gene_name": variants_df["gene_name"],
            }
        )

    def _standardize_cosmic_data(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Standardize COSMIC data."""
        return pd.DataFrame(
            {
                "position": variants_df["position"],
                "variant_id": variants_df["cosmic_id"],
                "reference": variants_df["reference"],
                "alternate": variants_df["alternate"],
                "source": "COSMIC",
                # COSMIC MUTATION_DESCRIPTION contains variant types (SNV/deletion/insertion), not functional impacts.
                # Label as "unclassified variant" so they pass through filtering to validation, where genomic
                # changes are applied to determine actual functional impact (missense/nonsense/frameshift/etc.).
                "impact": "unclassified variant",
                "hgvsc": variants_df["hgvs_coding"].apply(self._standardize_hgvsc),
                "hgvsp": variants_df["hgvs_protein"].apply(self._standardize_hgvsp),
                "allele_frequency": None,  # COSMIC doesn't have population frequency
                "allele_count": None,  # COSMIC doesn't have this field
                "allele_count_hom": None,  # COSMIC doesn't have this field
                "sample_count": variants_df["sample_count"]
                if "sample_count" in variants_df.columns
                else None,
                "clinical_significance": None,  # COSMIC doesn't have clinical significance
                "chromosome": variants_df["chromosome"],
                "gene_name": variants_df["gene_name"],
            }
        )

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
            "unclassified": "unclassified variant",
        }
        return mapping.get(category, "other variant")

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

    def _filter_out_low_impact_variants(
        self, mutations_df: pd.DataFrame
    ) -> pd.DataFrame:
        """Filter out low impact variants (splice, synonymous, intronic) by default.

        Args:
            mutations_df: DataFrame with mutations

        Returns:
            Filtered DataFrame excluding low impact variants
        """
        if mutations_df.empty:
            return mutations_df

        # Get all low impact terms
        low_impact_terms = []
        for category in self.low_impact_categories:
            low_impact_terms.extend(self.mutation_categories[category])

        # Add explicit terms to filter
        low_impact_terms.extend(["other variant", "unknown"])

        # Filter out mutations with low impact
        filtered_df = mutations_df.copy()
        for term in low_impact_terms:
            mask = ~filtered_df["impact"].str.contains(term, case=False, na=False)
            filtered_df = filtered_df[mask]

        return filtered_df

    def _bulk_filter_intronic_variants(
        self,
        mutations_df: pd.DataFrame,
        transcript_id: str,
        protein_generator,
        current_feature: Optional[pd.Series] = None,
    ) -> Tuple[pd.DataFrame, List[Dict]]:
        """Bulk filter out intronic variants using transcript region checks.

        Args:
            mutations_df: DataFrame with mutations
            transcript_id: Transcript ID to check regions against
            protein_generator: AlternativeProteinGenerator instance for region checks
            current_feature: Optional feature information

        Returns:
            Tuple of (filtered DataFrame, list of intronic variant info dicts)
        """
        if mutations_df.empty:
            return mutations_df, []

        intronic_variants = []
        non_intronic_indices = []

        # Check if debug mode is enabled
        debug_mode = getattr(protein_generator, "debug", False)

        for idx, mutation in mutations_df.iterrows():
            try:
                genomic_pos = int(mutation["position"])
                ref_allele = str(mutation.get("reference", "")).upper()
                variant_id = mutation.get("variant_id", f"pos_{genomic_pos}")
                source = mutation.get("source", "unknown")

                # Format position display (range for multi-bp variants)
                ref_len = len(ref_allele) if ref_allele else 1
                if ref_len > 1:
                    pos_display = f"{genomic_pos}-{genomic_pos + ref_len - 1}"
                else:
                    pos_display = str(genomic_pos)

                # Use the same position check as predict_consequence_fast
                is_in_region = protein_generator._is_position_in_transcript_regions(
                    transcript_id, genomic_pos, current_feature, ref_allele
                )

                if is_in_region:
                    non_intronic_indices.append(idx)
                    if debug_mode:
                        logger.debug(
                            f"✓ Position {pos_display} ({variant_id}, {source}) in coding/UTR"
                        )
                else:
                    # Track intronic variant for logging
                    intronic_variants.append(
                        {
                            "position": genomic_pos,
                            "variant_id": variant_id,
                            "source": source,
                        }
                    )
                    if debug_mode:
                        logger.debug(
                            f"🚫 Position {pos_display} ({variant_id}, {source}) intronic (filtered out)"
                        )

            except Exception as e:
                # If check fails, keep the variant (safer to validate than skip)
                logger.debug(f"Error checking position {mutation.get('position')}: {e}")
                non_intronic_indices.append(idx)

        # Filter to non-intronic mutations
        filtered_df = mutations_df.loc[non_intronic_indices].copy()

        return filtered_df, intronic_variants

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

    async def get_mutations_for_gene(
        self,
        gene_name: str,
        alt_features: Optional[pd.DataFrame] = None,
        sources: Optional[List[str]] = None,
        refseq_id: Optional[str] = None,
        custom_parquet_path: Optional[str] = None,
        verbose: bool = False,
    ) -> pd.DataFrame:
        """Get mutation data from specified sources for a gene.

        Fetches, standardizes, and combines mutation data from multiple sources
        (gnomAD, ClinVar, COSMIC, custom parquet files).

        Args:
            gene_name (str): Name of the gene.
            alt_features (Optional[pd.DataFrame]): DataFrame containing alternative isoform features.
            sources (Optional[List[str]]): List of mutation sources to query.
            refseq_id (Optional[str]): RefSeq transcript ID to filter mutations (used for ClinVar and other sources supporting transcript-level queries).
            custom_parquet_path (Optional[str]): Path to custom parquet file (if using 'custom' source).
            verbose (bool): Whether to print detailed information.

        Returns:
            pd.DataFrame: DataFrame containing standardized mutations from all requested sources.
        """
        if sources is None:
            sources = ["gnomad", "clinvar", "cosmic"]
        sources = [s.lower() for s in sources]

        # Fetch and combine data from sources
        dfs_to_combine = []
        if verbose:
            logger.info(f"Fetching mutations from sources: {', '.join(sources)}...")

        for source in sources:
            try:
                df = await self._fetch_source_data(source, gene_name, refseq_id)
                standardized_df = self.standardize_mutation_data(df, source, gene_name)

                if not standardized_df.empty:
                    dfs_to_combine.append(standardized_df)
                    if verbose:
                        logger.info(
                            f"Found {len(standardized_df)} variants from {source}"
                        )

            except Exception as e:
                logger.error(f"Error fetching {source} data: {str(e)}")

        # Add custom data if specified
        if "custom" in sources and custom_parquet_path:
            try:
                custom_df = await self.get_custom_variants(
                    custom_parquet_path, gene_name
                )
                if not custom_df.empty:
                    # Apply standardization to custom data
                    custom_df["impact"] = custom_df["impact"].apply(
                        self._standardize_impact_category
                    )
                    custom_df["hgvsc"] = custom_df["hgvsc"].apply(
                        self._standardize_hgvsc
                    )
                    custom_df["hgvsp"] = custom_df.apply(
                        lambda row: self._standardize_hgvsp(
                            row["hgvsp"], row["impact"]
                        ),
                        axis=1,
                    )
                    dfs_to_combine.append(custom_df)
                    logger.info(f"Found {len(custom_df)} variants from custom")
            except Exception as e:
                logger.error(f"Error fetching custom data: {str(e)}")

        if not dfs_to_combine:
            if verbose:
                logger.info("No mutations found")
            return pd.DataFrame()

        # Process combined data
        combined_df = pd.concat(dfs_to_combine, ignore_index=True)
        combined_df["position"] = pd.to_numeric(
            combined_df["position"], errors="coerce"
        )

        # Filter by alt_features if provided
        if alt_features is not None and not alt_features.empty:
            return self._filter_by_alt_features(combined_df, alt_features, verbose)
        else:
            # No filtering, just add relative positions
            if len(combined_df) > 0:
                min_pos = combined_df["position"].min()
                combined_df["relative_position"] = combined_df["position"] - min_pos
            return combined_df

    async def _fetch_source_data(
        self, source: str, gene_name: str, refseq_id: str = None
    ) -> pd.DataFrame:
        """Fetch data from a specific source."""
        if source == "gnomad":
            return await self.get_gnomad_variants(gene_name)
        elif source == "clinvar":
            return await self.get_clinvar_variants(gene_name, refseq_id)
        elif source == "cosmic":
            return await self.get_cosmic_variants(gene_name)
        else:
            return pd.DataFrame()

    def _filter_by_alt_features(
        self, combined_df: pd.DataFrame, alt_features: pd.DataFrame, verbose: bool
    ) -> pd.DataFrame:
        """Filter mutations by alternative features."""
        alt_features["start"] = pd.to_numeric(alt_features["start"])
        alt_features["end"] = pd.to_numeric(alt_features["end"])

        # Also convert mutation coordinates if they exist
        if "mutation_start" in alt_features.columns:
            alt_features["mutation_start"] = pd.to_numeric(
                alt_features["mutation_start"]
            )
        if "mutation_end" in alt_features.columns:
            alt_features["mutation_end"] = pd.to_numeric(alt_features["mutation_end"])

        if verbose:
            logger.info("Analyzing mutations in alternative features:")

        filtered_dfs = []

        for idx, feature in alt_features.iterrows():
            # Use mutation coordinates if available, otherwise fall back to original coordinates
            start_pos = feature.get("mutation_start", feature["start"])
            end_pos = feature.get("mutation_end", feature["end"])
            name = feature.get("name", f"Feature {idx + 1}")

            # Show both original and mutation coordinates if they differ
            if verbose and ("mutation_start" in feature or "mutation_end" in feature):
                orig_start = feature["start"]
                orig_end = feature["end"]
                if start_pos != orig_start or end_pos != orig_end:
                    logger.info(
                        f"Feature {name}: {orig_start}-{orig_end} extended to {start_pos}-{end_pos} for mutation detection"
                    )

            # Filter mutations for this feature using mutation coordinates
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
                logger.info(
                    f"Feature {name}({start_pos}-{end_pos}): {len(feature_mutations)} mutations"
                )

        if filtered_dfs:
            filtered_df = pd.concat(filtered_dfs, ignore_index=True)
            if verbose:
                logger.info(f"Total mutations in all features: {len(filtered_df)}")

            # Calculate relative positions
            if len(filtered_df) > 0:
                min_pos = filtered_df["position"].min()
                filtered_df["relative_position"] = filtered_df["position"] - min_pos
            return filtered_df
        else:
            if verbose:
                logger.info("No mutations found in any feature")
            return pd.DataFrame()

    async def _validate_mutation_consequences(
        self,
        mutations_df: pd.DataFrame,
        transcript_id: str,
        protein_generator,
        current_feature: Optional[pd.Series] = None,
    ) -> pd.DataFrame:
        """Validate mutation consequences using translation-based prediction."""
        validated_mutations = []

        # Get gene name from first mutation if available
        gene_name = (
            mutations_df.iloc[0].get("gene_name", "unknown")
            if not mutations_df.empty
            else "unknown"
        )

        logger.info(
            f"Validating {len(mutations_df)} mutations for {gene_name} ({transcript_id})"
        )

        if current_feature is not None:
            feature_type = current_feature.get("region_type", "unknown")
            feature_range = (
                f"{current_feature.get('start', '?')}-{current_feature.get('end', '?')}"
            )
            logger.info(f"  Region: {feature_type} at {feature_range}")

        validation_stats = {
            "total_processed": 0,
            "successful_validations": 0,
            "validation_failures": 0,
            "impact_agreements": 0,
            "impact_disagreements": 0,
            "alt_start_site_mutations": 0,
            "disagreement_patterns": {},
        }

        for idx, mutation in mutations_df.iterrows():
            validation_stats["total_processed"] += 1

            try:
                # Get mutation details for logging and validation
                genomic_pos = int(mutation["position"])
                ref_allele = str(mutation["reference"]).upper()
                alt_allele = str(mutation["alternate"]).upper()
                original_impact = mutation.get("impact", "unknown")
                variant_id = mutation.get("variant_id", f"pos_{genomic_pos}")
                source = mutation.get("source", "unknown")
                chrom = mutation.get("chromosome", "?")

                # Create searchable variant identifier for databases (truncate long indels for display)
                def truncate_allele(allele: str, max_len: int = 10) -> str:
                    """Truncate long alleles for readable logging."""
                    if len(allele) <= max_len:
                        return allele
                    return f"{allele[:max_len]}...[{len(allele)}bp]"

                # Full variant for downstream use
                variant_lookup = f"{chrom}:{genomic_pos}:{ref_allele}>{alt_allele}"

                # Truncated variant for logging
                ref_display = truncate_allele(ref_allele)
                alt_display = truncate_allele(alt_allele)
                variant_display = f"{chrom}:{genomic_pos}:{ref_display}>{alt_display}"

                # Log mutation being evaluated at the start
                logger.info(
                    f"  Evaluating: {source}: {variant_id} | {variant_display} | Database impact: {original_impact}"
                )

                # Check if mutation is in alternative start site
                in_alt_start_site = False
                is_alt_start_loss = False
                alt_start_note = ""

                if current_feature is not None:
                    # Extract info directly from current_feature
                    alt_start_pos = current_feature.get("alternative_start_pos", None)
                    start_codon = current_feature.get(
                        "start_codon", "AUG"
                    )  # RNA notation
                    strand = current_feature.get("strand", "+")  # Get strand info

                    if alt_start_pos is not None:
                        try:
                            alt_start_pos = int(alt_start_pos)

                            # Handle strand-specific codon coordinates
                            if strand == "+":
                                codon_start = alt_start_pos
                                codon_end = alt_start_pos + 2
                            else:  # strand == "-"
                                codon_start = alt_start_pos - 2
                                codon_end = alt_start_pos

                            if codon_start <= genomic_pos <= codon_end:
                                in_alt_start_site = True
                                validation_stats["alt_start_site_mutations"] += 1

                                # Check for start codon loss (SNVs only for now)
                                if (
                                    len(ref_allele) == 1
                                    and len(alt_allele) == 1
                                    and self.genome_handler
                                ):
                                    try:
                                        # Get reference codon from genome (3 bases)
                                        if strand == "+":
                                            ref_codon_seq = (
                                                self.genome_handler.get_sequence(
                                                    chrom, codon_start, codon_end
                                                )
                                            )
                                        else:
                                            ref_codon_seq = (
                                                self.genome_handler.get_sequence(
                                                    chrom, codon_start, codon_end
                                                )
                                            )
                                            # Reverse complement for minus strand
                                            complement = {
                                                "A": "T",
                                                "T": "A",
                                                "G": "C",
                                                "C": "G",
                                            }
                                            ref_codon_seq = "".join(
                                                complement.get(b, b)
                                                for b in reversed(ref_codon_seq)
                                            )

                                        # Convert DNA to RNA
                                        ref_codon_rna = ref_codon_seq.replace("T", "U")

                                        # Compute variant position within codon
                                        if strand == "+":
                                            variant_pos_in_codon = (
                                                genomic_pos - codon_start
                                            )
                                        else:
                                            variant_pos_in_codon = (
                                                codon_end - genomic_pos
                                            )

                                        # Get mutant codon
                                        mutant_codon_rna = get_mutant_codon(
                                            ref_codon_seq,
                                            variant_pos_in_codon,
                                            ref_allele,
                                            alt_allele,
                                            strand,
                                        )

                                        # Check if this is a start loss
                                        is_alt_start_loss = is_start_codon_loss(
                                            ref_codon_rna, mutant_codon_rna
                                        )

                                        if is_alt_start_loss:
                                            logger.info(
                                                f"  Start loss detected: {ref_codon_rna} → {mutant_codon_rna}"
                                            )
                                    except Exception as e:
                                        logger.debug(
                                            f"Could not compute start loss for variant: {e}"
                                        )

                                alt_start_note = (
                                    f" (in alternative start site {start_codon})"
                                )
                                logger.info(
                                    f"  {source}: {variant_id} | {variant_display} | Alt start site at {start_codon} codon"
                                )

                        except (ValueError, TypeError):
                            logger.debug(
                                f"Warning: Invalid alt_start_pos '{alt_start_pos}'"
                            )

                # Helper function to compute mutant codon
                def get_mutant_codon(
                    ref_codon_seq: str,
                    variant_pos_in_codon: int,
                    ref_base: str,
                    alt_base: str,
                    strand: str,
                ) -> str:
                    """Compute the mutant codon sequence.

                    Args:
                        ref_codon_seq: Reference codon sequence (DNA, 3 bases)
                        variant_pos_in_codon: Position within codon (0, 1, or 2)
                        ref_base: Reference base from variant
                        alt_base: Alternate base from variant
                        strand: Strand ('+' or '-')

                    Returns:
                        Mutant codon sequence (RNA notation)
                    """
                    if len(ref_codon_seq) != 3:
                        return ""

                    # Convert to list for mutation
                    mutant_codon = list(ref_codon_seq.upper())

                    # Apply mutation
                    mutant_codon[variant_pos_in_codon] = alt_base.upper()
                    mutant_codon_seq = "".join(mutant_codon)

                    # Convert DNA to RNA (T→U)
                    mutant_codon_rna = mutant_codon_seq.replace("T", "U")

                    return mutant_codon_rna

                # Use fast prediction
                validated_impact = await protein_generator.predict_consequence_fast(
                    transcript_id, genomic_pos, ref_allele, alt_allele, current_feature
                )

                # Add annotation for alternative start site mutations
                if in_alt_start_site:
                    validated_impact_with_note = f"{validated_impact}{alt_start_note}"
                else:
                    validated_impact_with_note = validated_impact

                # Create validated mutation record - PRESERVE ORIGINAL IMPACT
                validated_mutation = mutation.copy()

                # Add new columns
                validated_mutation["impact_validated"] = validated_impact
                validated_mutation["impact_validated_with_note"] = (
                    validated_impact_with_note
                )
                validated_mutation["in_alt_start_site"] = in_alt_start_site
                validated_mutation["is_alt_start_loss"] = is_alt_start_loss
                validated_mutation["alt_start_note"] = alt_start_note

                # Check if mutation is in canonical start site
                in_canonical_start_site = False
                is_canonical_start_loss = False
                canonical_start_codon_extracted = ""
                canonical_start_note = ""

                if current_feature is not None:
                    # Extract canonical start position from feature
                    canonical_start_pos = current_feature.get(
                        "canonical_start_pos", None
                    )
                    strand = current_feature.get("strand", "+")

                    if canonical_start_pos is not None:
                        try:
                            canonical_start_pos = int(canonical_start_pos)

                            # Handle strand-specific codon coordinates
                            if strand == "+":
                                codon_start = canonical_start_pos
                                codon_end = canonical_start_pos + 2
                            else:  # strand == "-"
                                codon_start = canonical_start_pos - 2
                                codon_end = canonical_start_pos

                            if codon_start <= genomic_pos <= codon_end:
                                in_canonical_start_site = True
                                validation_stats["canonical_start_site_mutations"] = (
                                    validation_stats.get(
                                        "canonical_start_site_mutations", 0
                                    )
                                    + 1
                                )

                                # Check for start codon loss (SNVs only for now)
                                if (
                                    len(ref_allele) == 1
                                    and len(alt_allele) == 1
                                    and self.genome_handler
                                ):
                                    try:
                                        # Get reference codon from genome (3 bases)
                                        if strand == "+":
                                            ref_codon_seq = (
                                                self.genome_handler.get_sequence(
                                                    chrom, codon_start, codon_end
                                                )
                                            )
                                        else:
                                            ref_codon_seq = (
                                                self.genome_handler.get_sequence(
                                                    chrom, codon_start, codon_end
                                                )
                                            )
                                            # Reverse complement for minus strand
                                            complement = {
                                                "A": "T",
                                                "T": "A",
                                                "G": "C",
                                                "C": "G",
                                            }
                                            ref_codon_seq = "".join(
                                                complement.get(b, b)
                                                for b in reversed(ref_codon_seq)
                                            )

                                        # Convert DNA to RNA and store
                                        ref_codon_rna = ref_codon_seq.replace("T", "U")
                                        canonical_start_codon_extracted = ref_codon_rna

                                        # Compute variant position within codon
                                        if strand == "+":
                                            variant_pos_in_codon = (
                                                genomic_pos - codon_start
                                            )
                                        else:
                                            variant_pos_in_codon = (
                                                codon_end - genomic_pos
                                            )

                                        # Get mutant codon
                                        mutant_codon_rna = get_mutant_codon(
                                            ref_codon_seq,
                                            variant_pos_in_codon,
                                            ref_allele,
                                            alt_allele,
                                            strand,
                                        )

                                        # Check if this is a start loss
                                        is_canonical_start_loss = is_start_codon_loss(
                                            ref_codon_rna, mutant_codon_rna
                                        )

                                        if is_canonical_start_loss:
                                            logger.info(
                                                f"  Canonical start loss detected: {ref_codon_rna} → {mutant_codon_rna}"
                                            )
                                    except Exception as e:
                                        logger.debug(
                                            f"Could not compute canonical start loss for variant: {e}"
                                        )

                                canonical_start_note = " (in canonical start site)"
                                logger.info(
                                    f"  {source}: {variant_id} | {variant_display} | Canonical start site"
                                )

                        except (ValueError, TypeError):
                            logger.debug(
                                f"Warning: Invalid canonical_start_pos '{canonical_start_pos}'"
                            )

                validated_mutation["in_canonical_start_site"] = in_canonical_start_site
                validated_mutation["is_canonical_start_loss"] = is_canonical_start_loss
                validated_mutation["canonical_start_codon_extracted"] = (
                    canonical_start_codon_extracted
                )
                validated_mutation["canonical_start_note"] = canonical_start_note
                validated_mutation["validation_method"] = "fast_codon_analysis"
                validated_mutation["validation_successful"] = True

                # Log intronic variants detected during validation (BUG - shouldn't reach here)
                if validated_impact == "intronic":
                    logger.warning(
                        f"⚠️ VALIDATION BUG: {source}: {variant_id} classified as intronic during validation - should have been filtered earlier!"
                    )

                # Track agreement/disagreement
                if original_impact == validated_impact:
                    validated_mutation["impact_agreement"] = True
                    validation_stats["impact_agreements"] += 1
                    # Log all validated changes at INFO level
                    alt_site_tag = " [alt start]" if in_alt_start_site else ""
                    logger.info(
                        f"  {source}: {variant_id} | {variant_display} | Validated: {validated_impact}{alt_site_tag}"
                    )
                else:
                    validated_mutation["impact_agreement"] = False
                    validation_stats["impact_disagreements"] += 1

                    # Track disagreement patterns
                    disagreement_key = f"{original_impact} → {validated_impact}"
                    if (
                        disagreement_key
                        not in validation_stats["disagreement_patterns"]
                    ):
                        validation_stats["disagreement_patterns"][disagreement_key] = 0
                    validation_stats["disagreement_patterns"][disagreement_key] += 1

                    # Log disagreements at INFO level with full database lookup info
                    alt_site_tag = " [alt start]" if in_alt_start_site else ""
                    logger.info(
                        f"  {source}: {variant_id} | {variant_display} | {original_impact} → {validated_impact}{alt_site_tag}"
                    )

                validated_mutations.append(validated_mutation)
                validation_stats["successful_validations"] += 1

            except Exception as e:
                logger.warning(
                    f"  {source}: {variant_id} | {variant_display} | Validation failed: {str(e)}"
                )

                # Keep original mutation with error info
                validated_mutation = mutation.copy()
                validated_mutation["impact_validated"] = "validation_failed"
                validated_mutation["impact_validated_with_note"] = "validation_failed"
                validated_mutation["in_alt_start_site"] = False
                validated_mutation["alt_start_note"] = ""
                validated_mutation["in_canonical_start_site"] = False
                validated_mutation["canonical_start_note"] = ""
                validated_mutation["validation_method"] = "fast_codon_analysis"
                validated_mutation["validation_successful"] = False
                validated_mutation["validation_error"] = str(e)
                validated_mutation["impact_agreement"] = False

                validated_mutations.append(validated_mutation)
                validation_stats["validation_failures"] += 1

        # Create final DataFrame
        result_df = pd.DataFrame(validated_mutations)

        # Print concise summary at INFO level
        total = validation_stats["total_processed"]
        agreements = validation_stats["impact_agreements"]
        disagreements = validation_stats["impact_disagreements"]
        failures = validation_stats["validation_failures"]
        alt_start = validation_stats["alt_start_site_mutations"]

        agreement_pct = (agreements / total * 100) if total > 0 else 0
        logger.info(
            f"Validation complete: {total} total, {agreements} agree ({agreement_pct:.1f}%), "
            f"{disagreements} disagree, {failures} failed, {alt_start} in alt start sites"
        )

        # Show disagreement patterns if any exist
        if validation_stats["disagreement_patterns"]:
            logger.info(f"Disagreement breakdown:")
            for pattern, count in sorted(
                validation_stats["disagreement_patterns"].items(),
                key=lambda x: x[1],
                reverse=True,
            ):
                logger.info(f"  {pattern}: {count}")

        # Show final distributions
        if not result_df.empty:
            logger.debug(f"Final validated impact distribution:")
            if "impact_validated" in result_df.columns:
                for impact, count in (
                    result_df["impact_validated"].value_counts().items()
                ):
                    logger.debug(f"{impact}: {count}")

            # Show alternative start site breakdown
            if "in_alt_start_site" in result_df.columns:
                alt_start_count = result_df["in_alt_start_site"].sum()
                total_count = len(result_df)
                percentage = (
                    (alt_start_count / total_count * 100) if total_count > 0 else 0
                )
                logger.debug(
                    f"Alternative start site mutations: {alt_start_count}/{total_count}({percentage:.1f}%)"
                )

        return result_df

    async def analyze_gene_mutations_comprehensive(
        self,
        gene_name: str,
        genome_handler,
        alt_isoform_handler,
        output_dir: str,
        visualize: bool = False,
        impact_types: Optional[Dict[str, List[str]]] = None,
        custom_parquet_path: Optional[str] = None,
        sources: Optional[List[str]] = None,
        top_n_per_type_per_transcript: Optional[int] = None,
        validate_consequences: bool = True,
    ) -> Dict:
        """Comprehensive mutation analysis for a gene with transcript-alternative feature pairs.

        Args:
            gene_name: Name of the gene to analyze
            genome_handler: Initialized GenomeHandler instance
            alt_isoform_handler: Initialized AlternativeIsoform instance
            output_dir: Directory to save output files
            visualize: Whether to generate visualizations
            impact_types: Dict mapping sources to impact types to filter by
            custom_parquet_path: Path to custom parquet file (if using 'custom' source)
            sources: Data sources to pull mutations from
            top_n_per_type_per_transcript: Top N features per type per transcript
            validate_consequences: Whether to validate consequences via translation
        Returns:
            Dictionary containing comprehensive analysis results
        """
        try:
            # Get alternative isoform features
            logger.debug(f"  Getting alternative features...")
            alt_features = alt_isoform_handler.get_mutation_features(
                gene_name, top_n_per_type_per_transcript
            )
            if alt_features.empty:
                logger.debug(f"  No alternative features found")
                return {"gene_name": gene_name, "status": "no_features", "error": None}
            logger.debug(f"  Found {len(alt_features)} alternative features")

            # Get unique transcript IDs from the BED file (alt_features)
            logger.debug(f"  Using transcript IDs from BED file...")
            transcript_ids_from_bed = alt_features["transcript_id"].unique()
            # Filter out any "NA" assignments
            transcript_ids_from_bed = [
                tid for tid in transcript_ids_from_bed if tid != "NA"
            ]
            if not transcript_ids_from_bed:
                logger.debug(f"  No valid transcript assignments in BED file")
                return {
                    "gene_name": gene_name,
                    "status": "no_transcript_assignments",
                    "error": None,
                }

            # Get transcript info for only the transcripts specified in BED file
            transcript_info = genome_handler.get_transcript_ids(gene_name)
            if transcript_info.empty:
                logger.debug(f"  No transcript info found in genome")
                return {
                    "gene_name": gene_name,
                    "status": "no_transcripts",
                    "error": None,
                }

            # Filter to only transcripts that are in the BED file
            transcript_info = transcript_info[
                transcript_info["transcript_id"].isin(transcript_ids_from_bed)
            ]
            if transcript_info.empty:
                logger.debug(
                    f"  None of the BED transcript IDs found in genome annotation"
                )
                return {
                    "gene_name": gene_name,
                    "status": "no_matching_transcripts",
                    "error": None,
                }
            logger.debug(
                f"  Using {len(transcript_info)} transcripts from BED file assignments"
            )

            # Create transcript-alternative feature pairs directly from BED file assignments
            logger.debug(f"  Creating transcript-feature pairs from BED assignments...")
            transcript_feature_pairs = []
            # Process each alternative feature (which already has transcript assignment)
            for idx, feature in alt_features.iterrows():
                feature_transcript_id = feature["transcript_id"]
                # Skip if no transcript assignment
                if feature_transcript_id == "NA":
                    continue
                # Get transcript info for this specific transcript
                transcript_match = transcript_info[
                    transcript_info["transcript_id"] == feature_transcript_id
                ]
                if transcript_match.empty:
                    # This transcript ID from BED wasn't found in genome annotation
                    continue
                # Get refseq_id directly from feature if present
                refseq_id = feature.get("refseq_transcript_id", None)
                if refseq_id == "NA":
                    refseq_id = None
                transcript = transcript_match.iloc[0]
                transcript_start = transcript["start"]
                transcript_end = transcript["end"]
                transcript_chromosome = transcript["chromosome"]
                transcript_strand = transcript["strand"]
                feature_start = feature.get("mutation_start", feature["start"])
                feature_end = feature.get("mutation_end", feature["end"])
                feature_chrom = feature["chromosome"]
                feature_type = feature.get("region_type", "unknown")
                # Create feature ID
                feature_id = f"{feature_type}_{idx}"
                if "start_codon" in feature and not pd.isna(feature["start_codon"]):
                    feature_id = f"{feature_type}_{feature['start_codon']}_{feature_start}_{feature_end}"
                # Create the transcript-feature pair
                transcript_feature_pairs.append(
                    {
                        "transcript_id": feature_transcript_id,
                        "refseq_id": refseq_id,
                        "feature_id": feature_id,
                        "bed_name": feature.get(
                            "name", feature_id
                        ),  # Add BED name for pre-validated lookup
                        "feature_idx": idx,
                        "feature_type": feature_type,
                        "transcript_start": transcript_start,
                        "transcript_end": transcript_end,
                        "transcript_strand": transcript_strand,
                        "feature_start": feature_start,
                        "feature_end": feature_end,
                    }
                )

            if not transcript_feature_pairs:
                logger.debug(
                    f"  No valid transcript-feature pairs from BED assignments"
                )
                return {
                    "gene_name": gene_name,
                    "status": "no_valid_pairs",
                    "error": None,
                }
            logger.debug(
                f"  Created {len(transcript_feature_pairs)} transcript-feature pairs from BED assignments"
            )

            # Determine sources to use
            if sources is None:
                sources = ["clinvar"]

            # Get the list of desired impact types for each source
            desired_impact_types = []
            if impact_types:
                for source, impacts in impact_types.items():
                    desired_impact_types.extend(impacts)

            # Initialize validation components once if needed
            protein_generator = None
            if validate_consequences:
                logger.info(f"Initializing mutation validation for {gene_name}")
                from swissisoform.translation import AlternativeProteinGenerator

                protein_generator = AlternativeProteinGenerator(
                    genome_handler=genome_handler,
                    alt_isoform_handler=alt_isoform_handler,
                    output_dir=output_dir,
                    debug=True,
                )

            # Container for all transcript-feature analysis results
            pair_results = []
            all_mutations_for_viz = pd.DataFrame()

            # Track errors and warnings for this gene
            gene_errors = []
            gene_warnings = []
            skipped_pairs = 0

            # Process each transcript-feature pair individually
            logger.debug(
                f"Processing {len(transcript_feature_pairs)} transcript-feature pairs..."
            )
            for pair_idx, pair in enumerate(transcript_feature_pairs, 1):
                transcript_id = pair["transcript_id"]
                feature_id = pair["feature_id"]
                feature_idx = pair["feature_idx"]
                feature_type = pair["feature_type"]
                refseq_id = pair.get("refseq_id", None)
                feature_start = pair["feature_start"]
                feature_end = pair["feature_end"]

                logger.debug(
                    f"Processing pair {pair_idx}/{len(transcript_feature_pairs)}: {transcript_id} × {feature_id}"
                )

                # STEP 1: Get mutations specifically for this transcript-feature pair
                logger.debug(f"Fetching mutations for {refseq_id or 'gene-level'}...")
                pair_mutations = await self.get_mutations_for_gene(
                    gene_name=gene_name,
                    alt_features=alt_features.loc[[feature_idx]],  # Only this feature
                    sources=sources,
                    refseq_id=refseq_id,
                    custom_parquet_path=custom_parquet_path,
                    verbose=True,
                )

                if pair_mutations is None or pair_mutations.empty:
                    logger.debug(f"No mutations found for this pair")
                    gene_warnings.append(f"no_mutations:{transcript_id}")
                    # Still record the pair with zero counts
                    mutation_categories = {}
                    if impact_types:
                        for impact_type in desired_impact_types:
                            category_key = (
                                f"mutations_{impact_type.replace(' ', '_').lower()}"
                            )
                            mutation_categories[category_key] = 0

                    # Initialize per-source counts
                    source_categories = {}
                    if sources:
                        for source in sources:
                            source_key = f"mutations_{source.lower()}"
                            source_categories[source_key] = 0

                    pair_results.append(
                        {
                            "transcript_id": transcript_id,
                            "feature_id": feature_id,
                            "feature_type": feature_type,
                            "feature_start": feature_start,
                            "feature_end": feature_end,
                            "mutation_count_total": 0,
                            "mutation_sources": "",  # NEW
                            "variant_ids": "",
                            "mutations_in_alt_start_site": 0,
                            "alt_start_site_variant_ids": "",
                            "mutations_in_canonical_start_site": 0,
                            "canonical_start_site_variant_ids": "",
                            "alternative_start_codon": "",
                            "canonical_start_codon": "AUG",
                            "alternative_start_loss_count": 0,
                            "alternative_start_loss_variant_ids": "",
                            "canonical_start_loss_count": 0,
                            "canonical_start_loss_variant_ids": "",
                            "gnomad_allele_count": 0,
                            "cosmic_sample_count": 0,
                            **mutation_categories,
                            **source_categories,  # NEW
                        }
                    )
                    continue

                logger.debug(f"Found {len(pair_mutations)} mutations")

                # STEP 2: Filter to mutations in the specific alternative region
                # For multi-bp variants, ensure the ENTIRE variant is within the feature
                logger.debug(
                    f"Filtering to feature region ({feature_start}-{feature_end})..."
                )

                # Calculate end position for each variant based on reference allele length
                def get_variant_end_pos(row):
                    ref_len = (
                        len(str(row.get("reference", "")))
                        if pd.notna(row.get("reference"))
                        else 1
                    )
                    return row["position"] + ref_len - 1

                pair_mutations_copy = pair_mutations.copy()
                pair_mutations_copy["variant_end_pos"] = pair_mutations_copy.apply(
                    get_variant_end_pos, axis=1
                )

                # Filter: variant start must be >= feature_start AND variant end must be <= feature_end
                region_mutations = pair_mutations_copy[
                    (pair_mutations_copy["position"] >= feature_start)
                    & (pair_mutations_copy["variant_end_pos"] <= feature_end)
                ].copy()

                # Log boundary-spanning variants that were filtered out
                boundary_filtered = pair_mutations_copy[
                    ~(
                        (pair_mutations_copy["position"] >= feature_start)
                        & (pair_mutations_copy["variant_end_pos"] <= feature_end)
                    )
                ]
                if not boundary_filtered.empty:
                    logger.debug(
                        f"Filtered {len(boundary_filtered)} boundary-spanning variant(s):"
                    )
                    for _, var in boundary_filtered.iterrows():
                        var_id = var.get("variant_id", "unknown")
                        var_pos = var.get("position", "?")
                        var_end = var.get("variant_end_pos", "?")
                        var_ref = var.get("reference", "?")
                        var_source = var.get("source", "?")
                        logger.debug(
                            f"  🚫 {var_id} ({var_source}): pos {var_pos}-{var_end} extends beyond feature boundary ({feature_start}-{feature_end}), ref={var_ref}"
                        )

                # Drop the temporary column
                if "variant_end_pos" in region_mutations.columns:
                    region_mutations = region_mutations.drop(
                        columns=["variant_end_pos"]
                    )

                if region_mutations.empty:
                    logger.debug(f"No mutations in feature region")
                    mutation_categories = {}
                    if impact_types:
                        for impact_type in desired_impact_types:
                            category_key = (
                                f"mutations_{impact_type.replace(' ', '_').lower()}"
                            )
                            mutation_categories[category_key] = 0

                    # Initialize per-source counts
                    source_categories = {}
                    if sources:
                        for source in sources:
                            source_key = f"mutations_{source.lower()}"
                            source_categories[source_key] = 0

                    pair_results.append(
                        {
                            "transcript_id": transcript_id,
                            "feature_id": feature_id,
                            "feature_type": feature_type,
                            "feature_start": feature_start,
                            "feature_end": feature_end,
                            "mutation_count_total": 0,
                            "mutation_sources": "",  # NEW
                            "variant_ids": "",
                            "mutations_in_alt_start_site": 0,
                            "alt_start_site_variant_ids": "",
                            "mutations_in_canonical_start_site": 0,
                            "canonical_start_site_variant_ids": "",
                            "alternative_start_codon": "",
                            "canonical_start_codon": "AUG",
                            "alternative_start_loss_count": 0,
                            "alternative_start_loss_variant_ids": "",
                            "canonical_start_loss_count": 0,
                            "canonical_start_loss_variant_ids": "",
                            "gnomad_allele_count": 0,
                            "cosmic_sample_count": 0,
                            **mutation_categories,
                            **source_categories,  # NEW
                        }
                    )
                    continue

                logger.debug(f"{len(region_mutations)} mutations in feature region")

                # STEP 3: Bulk filter out intronic variants
                logger.debug(f"Filtering out intronic variants...")

                # Get current feature for context (needed for position checks)
                current_feature = (
                    alt_features.loc[feature_idx]
                    if feature_idx in alt_features.index
                    else None
                )

                # Get CDS and UTR regions for debug logging
                if protein_generator is not None:
                    try:
                        features = protein_generator.genome.get_transcript_features(
                            transcript_id
                        )
                        transcript_data = protein_generator.genome.get_transcript_features_with_sequence(
                            transcript_id
                        )
                        strand = (
                            transcript_data["sequence"]["strand"]
                            if transcript_data
                            else "+"
                        )
                        chrom = (
                            transcript_data["sequence"]["chromosome"]
                            if transcript_data
                            else "?"
                        )

                        cds_regions = features[features["feature_type"] == "CDS"]
                        utr_regions = features[
                            features["feature_type"].str.contains("UTR", na=False)
                        ]

                        # Combine and sort all regions by position (strand-aware)
                        all_regions = []
                        for _, cds in cds_regions.iterrows():
                            all_regions.append(
                                {
                                    "start": int(cds["start"]),
                                    "end": int(cds["end"]),
                                    "type": "CDS",
                                }
                            )
                        for _, utr in utr_regions.iterrows():
                            all_regions.append(
                                {
                                    "start": int(utr["start"]),
                                    "end": int(utr["end"]),
                                    "type": "UTR",
                                }
                            )

                        # Sort by position (ascending for +, descending for -)
                        all_regions.sort(
                            key=lambda x: x["start"], reverse=(strand == "-")
                        )

                        # Print regions once
                        logger.debug(
                            f"Transcript regions ({strand} strand, {len(all_regions)} total):"
                        )
                        for region in all_regions:
                            logger.debug(
                                f"  {chrom}:{region['start']}-{region['end']} ({region['type']})"
                            )

                        # Sort mutations by position (same order as regions)
                        region_mutations_sorted = region_mutations.sort_values(
                            "position", ascending=(strand != "-")
                        ).copy()

                        # Perform bulk intronic filtering
                        high_impact_mutations, intronic_variants = (
                            self._bulk_filter_intronic_variants(
                                region_mutations_sorted,
                                transcript_id,
                                protein_generator,
                                current_feature,
                            )
                        )

                        # Log breakdown by source
                        if intronic_variants:
                            source_counts = {}
                            for v in intronic_variants:
                                source = v["source"]
                                source_counts[source] = source_counts.get(source, 0) + 1
                            source_summary = ", ".join(
                                [f"{src}: {cnt}" for src, cnt in source_counts.items()]
                            )
                            logger.debug(
                                f"Intronic variants by source (filtered out): {source_summary}"
                            )

                        logger.debug(
                            f"{len(high_impact_mutations)} variants remain after intronic filtering"
                        )

                    except Exception as e:
                        logger.debug(
                            f"Error during intronic filtering: {e}, keeping all variants"
                        )
                        high_impact_mutations = region_mutations.copy()
                else:
                    # No protein generator available, skip filtering
                    high_impact_mutations = region_mutations.copy()

                if high_impact_mutations.empty:
                    logger.debug(f"No high-impact mutations remain")
                    gene_warnings.append(f"all_mutations_filtered:{transcript_id}")
                    mutation_categories = {}
                    if impact_types:
                        for impact_type in desired_impact_types:
                            category_key = (
                                f"mutations_{impact_type.replace(' ', '_').lower()}"
                            )
                            mutation_categories[category_key] = 0

                    # Initialize per-source counts
                    source_categories = {}
                    if sources:
                        for source in sources:
                            source_key = f"mutations_{source.lower()}"
                            source_categories[source_key] = 0

                    pair_results.append(
                        {
                            "transcript_id": transcript_id,
                            "feature_id": feature_id,
                            "feature_type": feature_type,
                            "feature_start": feature_start,
                            "feature_end": feature_end,
                            "mutation_count_total": 0,
                            "mutation_sources": "",
                            "variant_ids": "",
                            **mutation_categories,
                            **source_categories,
                        }
                    )
                    continue

                # STEP 4: Validate consequences for this specific subset
                validated_mutations = high_impact_mutations.copy()
                if validate_consequences and protein_generator is not None:
                    logger.debug(
                        f"Validating consequences for {len(high_impact_mutations)} mutations..."
                    )

                    # Simple protein validation check
                    if current_feature is not None:
                        test_result = protein_generator.extract_alternative_protein(
                            transcript_id, current_feature
                        )
                        if not test_result:
                            logger.debug(
                                f"🛑 Protein extraction failed - skipping this pair"
                            )
                            gene_errors.append(
                                f"protein_extraction_failed:{transcript_id}"
                            )
                            skipped_pairs += 1
                            continue

                    validated_mutations = await self._validate_mutation_consequences(
                        high_impact_mutations,
                        transcript_id,
                        protein_generator,
                        current_feature,
                    )
                    logger.debug(f"Validation complete")

                    # Filter out intronic variants detected during validation (shouldn't happen but track it)
                    if "impact_validated" in validated_mutations.columns:
                        intronic_mask = (
                            validated_mutations["impact_validated"] == "intronic"
                        )
                        intronic_count = intronic_mask.sum()
                        if intronic_count > 0:
                            # Track the variant IDs that were intronic
                            intronic_variants = validated_mutations[intronic_mask]
                            for _, row in intronic_variants.iterrows():
                                variant_id = row.get("variant_id", "unknown")
                                gene_errors.append(
                                    f"intronic_in_validation:{variant_id}"
                                )

                            logger.warning(
                                f"Filtering {intronic_count} intronic variants detected during validation"
                            )
                            validated_mutations = validated_mutations[~intronic_mask]

                    # Apply user-specified impact filtering on VALIDATED impacts
                    if (
                        impact_types
                        and "impact_validated" in validated_mutations.columns
                    ):
                        pre_filter_count = len(validated_mutations)
                        validated_mutations = validated_mutations[
                            validated_mutations["impact_validated"].isin(
                                desired_impact_types
                            )
                        ]
                        logger.debug(
                            f"User impact filtering: {pre_filter_count}→{len(validated_mutations)} mutations"
                        )

                # STEP 5: Collect results for this pair
                pair_mutation_count = len(validated_mutations)

                # Initialize default counts for all desired impact types
                mutation_categories = {}
                source_categories = {}  # Track per-source counts
                source_impact_categories = {}  # Track per-source per-impact counts (e.g., clinvar_missense)

                if impact_types:
                    for impact_type in desired_impact_types:
                        category_key = (
                            f"mutations_{impact_type.replace(' ', '_').lower()}"
                        )
                        mutation_categories[category_key] = 0

                # Initialize per-source counts for each requested source
                if sources:
                    for source in sources:
                        source_key = f"mutations_{source.lower()}"
                        source_categories[source_key] = 0

                        # Initialize per-source per-impact counts (source×impact matrix)
                        if impact_types:
                            for impact_type in desired_impact_types:
                                impact_normalized = impact_type.replace(
                                    " ", "_"
                                ).lower()
                                source_impact_key = (
                                    f"mutations_{source.lower()}_{impact_normalized}"
                                )
                                source_impact_categories[source_impact_key] = 0
                                # Also track IDs for source×impact
                                source_impact_id_key = (
                                    f"variant_ids_{source.lower()}_{impact_normalized}"
                                )
                                source_impact_categories[source_impact_id_key] = ""

                variant_ids = []
                mutation_sources = set()  # Track which sources have mutations
                if pair_mutation_count > 0:
                    # Count by impact category
                    impact_column = (
                        "impact_validated"
                        if "impact_validated" in validated_mutations.columns
                        else "impact"
                    )

                    for impact in validated_mutations[impact_column].unique():
                        if pd.isna(impact):
                            continue

                        impact_mutations = validated_mutations[
                            validated_mutations[impact_column] == impact
                        ]
                        category_count = len(impact_mutations)

                        # Store count for this impact type
                        category_key = f"mutations_{impact.replace(' ', '_').lower()}"
                        if category_key in mutation_categories:
                            mutation_categories[category_key] = category_count

                    # Count by source and by source×impact
                    # Also sum allele counts (gnomAD) and sample counts (COSMIC)
                    gnomad_allele_count_sum = 0
                    cosmic_sample_count_sum = 0

                    if "source" in validated_mutations.columns:
                        for source in validated_mutations["source"].unique():
                            if pd.notna(source):
                                mutation_sources.add(str(source))
                                source_mutations = validated_mutations[
                                    validated_mutations["source"] == source
                                ]
                                source_key = f"mutations_{str(source).lower()}"
                                source_categories[source_key] = len(source_mutations)

                                # Sum allele counts for gnomAD
                                if (
                                    str(source).lower() == "gnomad"
                                    and "allele_count" in source_mutations.columns
                                ):
                                    gnomad_allele_count_sum = (
                                        source_mutations["allele_count"].fillna(0).sum()
                                    )

                                # Sum sample counts for COSMIC
                                if (
                                    str(source).lower() == "cosmic"
                                    and "sample_count" in source_mutations.columns
                                ):
                                    cosmic_sample_count_sum = (
                                        source_mutations["sample_count"].fillna(0).sum()
                                    )

                                # Count by source×impact (e.g., clinvar_missense_variant)
                                for impact in source_mutations[impact_column].unique():
                                    if pd.notna(impact):
                                        source_impact_mutations = source_mutations[
                                            source_mutations[impact_column] == impact
                                        ]
                                        impact_normalized = impact.replace(
                                            " ", "_"
                                        ).lower()
                                        source_impact_key = f"mutations_{str(source).lower()}_{impact_normalized}"

                                        if (
                                            source_impact_key
                                            in source_impact_categories
                                        ):
                                            source_impact_categories[
                                                source_impact_key
                                            ] = len(source_impact_mutations)

                                            # Also collect IDs for source×impact
                                            if (
                                                "variant_id"
                                                in source_impact_mutations.columns
                                            ):
                                                source_impact_ids = (
                                                    source_impact_mutations[
                                                        "variant_id"
                                                    ]
                                                    .dropna()
                                                    .unique()
                                                    .tolist()
                                                )
                                                source_impact_ids = [
                                                    str(id).strip()
                                                    for id in source_impact_ids
                                                    if str(id).strip()
                                                ]
                                                if source_impact_ids:
                                                    source_impact_id_key = f"variant_ids_{str(source).lower()}_{impact_normalized}"
                                                    source_impact_categories[
                                                        source_impact_id_key
                                                    ] = ",".join(source_impact_ids)

                    # Collect all variant IDs
                    if "variant_id" in validated_mutations.columns:
                        variant_ids = (
                            validated_mutations["variant_id"].dropna().unique().tolist()
                        )
                        variant_ids = [
                            str(id).strip() for id in variant_ids if str(id).strip()
                        ]

                    # Add to collection for final visualization
                    validated_mutations = validated_mutations.copy()
                    validated_mutations["pair_transcript_id"] = transcript_id
                    validated_mutations["pair_feature_id"] = feature_id
                    all_mutations_for_viz = pd.concat(
                        [all_mutations_for_viz, validated_mutations], ignore_index=True
                    )

                # Create result record for this pair
                pair_result = {
                    "transcript_id": transcript_id,
                    "feature_id": feature_id,
                    "bed_name": pair.get("bed_name", ""),  # BED feature name for lookup
                    "feature_type": feature_type,
                    "feature_start": feature_start,
                    "feature_end": feature_end,
                    "mutation_count_total": pair_mutation_count,
                    "mutation_sources": ",".join(sorted(mutation_sources))
                    if mutation_sources
                    else "",  # List of sources
                    "variant_ids": ",".join(variant_ids) if variant_ids else "",
                    "mutations_in_alt_start_site": 0,
                    "alt_start_site_variant_ids": "",
                    "mutations_in_canonical_start_site": 0,
                    "canonical_start_site_variant_ids": "",
                    "alternative_start_codon": "",
                    "canonical_start_codon": "AUG",
                    "alternative_start_loss_count": 0,
                    "alternative_start_loss_variant_ids": "",
                    "canonical_start_loss_count": 0,
                    "canonical_start_loss_variant_ids": "",
                    "gnomad_allele_count": int(gnomad_allele_count_sum),
                    "cosmic_sample_count": int(cosmic_sample_count_sum),
                    **mutation_categories,
                    **source_categories,  # Add per-source counts
                    **source_impact_categories,  # Add source×impact matrix (e.g., clinvar_missense_variant)
                }

                # Add validation metrics if available
                if validate_consequences and len(validated_mutations) > 0:
                    # Check for alternative start site mutations
                    if "in_alt_start_site" in validated_mutations.columns:
                        alt_start_mutations = validated_mutations[
                            validated_mutations["in_alt_start_site"] == True
                        ]
                        pair_result["mutations_in_alt_start_site"] = len(
                            alt_start_mutations
                        )

                        if (
                            not alt_start_mutations.empty
                            and "variant_id" in alt_start_mutations.columns
                        ):
                            alt_start_ids = (
                                alt_start_mutations["variant_id"]
                                .dropna()
                                .unique()
                                .tolist()
                            )
                            alt_start_ids = [
                                str(id).strip()
                                for id in alt_start_ids
                                if str(id).strip()
                            ]
                            pair_result["alt_start_site_variant_ids"] = (
                                ",".join(alt_start_ids) if alt_start_ids else ""
                            )

                        # Extract alternative start codon
                        if (
                            current_feature is not None
                            and not alt_start_mutations.empty
                        ):
                            alt_codon = current_feature.get("start_codon", "")
                            pair_result["alternative_start_codon"] = alt_codon

                        # Check for alternative start losses
                        if "is_alt_start_loss" in validated_mutations.columns:
                            alt_start_loss_mutations = validated_mutations[
                                (validated_mutations["in_alt_start_site"] == True)
                                & (validated_mutations["is_alt_start_loss"] == True)
                            ]
                            pair_result["alternative_start_loss_count"] = len(
                                alt_start_loss_mutations
                            )

                            if (
                                not alt_start_loss_mutations.empty
                                and "variant_id" in alt_start_loss_mutations.columns
                            ):
                                alt_start_loss_ids = (
                                    alt_start_loss_mutations["variant_id"]
                                    .dropna()
                                    .unique()
                                    .tolist()
                                )
                                alt_start_loss_ids = [
                                    str(id).strip()
                                    for id in alt_start_loss_ids
                                    if str(id).strip()
                                ]
                                pair_result["alternative_start_loss_variant_ids"] = (
                                    ",".join(alt_start_loss_ids)
                                    if alt_start_loss_ids
                                    else ""
                                )

                                # Log start losses
                                if alt_start_loss_ids:
                                    logger.info(
                                        f"Alternative start codon loss: {alt_codon} → variant(s), count={len(alt_start_loss_ids)}, IDs={','.join(alt_start_loss_ids[:3])}{'...' if len(alt_start_loss_ids) > 3 else ''}"
                                    )

                    # Check for canonical start site mutations
                    if "in_canonical_start_site" in validated_mutations.columns:
                        canonical_start_mutations = validated_mutations[
                            validated_mutations["in_canonical_start_site"] == True
                        ]
                        pair_result["mutations_in_canonical_start_site"] = len(
                            canonical_start_mutations
                        )

                        if (
                            not canonical_start_mutations.empty
                            and "variant_id" in canonical_start_mutations.columns
                        ):
                            canonical_start_ids = (
                                canonical_start_mutations["variant_id"]
                                .dropna()
                                .unique()
                                .tolist()
                            )
                            canonical_start_ids = [
                                str(id).strip()
                                for id in canonical_start_ids
                                if str(id).strip()
                            ]
                            pair_result["canonical_start_site_variant_ids"] = (
                                ",".join(canonical_start_ids)
                                if canonical_start_ids
                                else ""
                            )

                        # Extract canonical start codon if available
                        if (
                            not canonical_start_mutations.empty
                            and "canonical_start_codon_extracted"
                            in canonical_start_mutations.columns
                        ):
                            # Get first non-empty extracted codon
                            extracted_codons = canonical_start_mutations[
                                "canonical_start_codon_extracted"
                            ].dropna()
                            extracted_codons = extracted_codons[extracted_codons != ""]
                            if not extracted_codons.empty:
                                pair_result["canonical_start_codon"] = (
                                    extracted_codons.iloc[0]
                                )

                        # Check for canonical start losses
                        if "is_canonical_start_loss" in validated_mutations.columns:
                            canonical_start_loss_mutations = validated_mutations[
                                (validated_mutations["in_canonical_start_site"] == True)
                                & (
                                    validated_mutations["is_canonical_start_loss"]
                                    == True
                                )
                            ]
                            pair_result["canonical_start_loss_count"] = len(
                                canonical_start_loss_mutations
                            )

                            if (
                                not canonical_start_loss_mutations.empty
                                and "variant_id"
                                in canonical_start_loss_mutations.columns
                            ):
                                canonical_start_loss_ids = (
                                    canonical_start_loss_mutations["variant_id"]
                                    .dropna()
                                    .unique()
                                    .tolist()
                                )
                                canonical_start_loss_ids = [
                                    str(id).strip()
                                    for id in canonical_start_loss_ids
                                    if str(id).strip()
                                ]
                                pair_result["canonical_start_loss_variant_ids"] = (
                                    ",".join(canonical_start_loss_ids)
                                    if canonical_start_loss_ids
                                    else ""
                                )

                                # Log start losses
                                if canonical_start_loss_ids:
                                    logger.info(
                                        f"Canonical start codon loss: AUG → variant(s), count={len(canonical_start_loss_ids)}, IDs={','.join(canonical_start_loss_ids[:3])}{'...' if len(canonical_start_loss_ids) > 3 else ''}"
                                    )

                    # Original validation metrics (if they exist)
                    if "consequence_matches" in validated_mutations.columns:
                        total_validated = (
                            validated_mutations["consequence_matches"].notna().sum()
                        )
                        matches = validated_mutations["consequence_matches"].sum()
                        pair_result["consequences_validated"] = total_validated
                        pair_result["consequences_match_reported"] = matches
                        pair_result["validation_accuracy"] = (
                            matches / total_validated if total_validated > 0 else 0
                        )

                pair_results.append(pair_result)

                # Enhanced completion message
                completion_msg = (
                    f"  Pair complete: {pair_mutation_count} final mutations"
                )
                if (
                    validate_consequences
                    and len(validated_mutations) > 0
                    and "in_alt_start_site" in validated_mutations.columns
                ):
                    alt_start_count = validated_mutations["in_alt_start_site"].sum()
                    if alt_start_count > 0:
                        completion_msg += f" (🎯 {alt_start_count} in alt start site)"
                logger.debug(completion_msg)

            logger.debug(
                f"All pairs processed. Total mutations for visualization: {len(all_mutations_for_viz)}"
            )

            # Generate visualizations if requested
            if visualize:
                from swissisoform.visualize import GenomeVisualizer

                visualizer = GenomeVisualizer(genome_handler)
                gene_dir = Path(output_dir) / gene_name
                gene_dir.mkdir(parents=True, exist_ok=True)
                logger.debug(
                    f"Generating visualizations for each transcript-feature pair:"
                )

                for pair_idx, pair in enumerate(transcript_feature_pairs, 1):
                    transcript_id = pair["transcript_id"]
                    feature_id = pair["feature_id"]
                    feature_idx = pair["feature_idx"]

                    # Get the specific feature for this pair
                    if feature_idx in alt_features.index:
                        feature_for_viz = alt_features.loc[[feature_idx]].copy()
                    else:
                        logger.debug(
                            f"Warning: Invalid feature index {feature_idx}, skipping visualization"
                        )
                        continue

                    logger.debug(
                        f"Visualizing pair {pair_idx}/{len(transcript_feature_pairs)}: {transcript_id} × {feature_id}"
                    )

                    # Create directories organized by transcript and feature
                    transcript_dir = gene_dir / transcript_id
                    transcript_dir.mkdir(exist_ok=True)

                    # Prepare output paths
                    pair_base_filename = f"{transcript_id}_{feature_id}"

                    # Always create visualization, even with empty mutations
                    pair_viz_mutations = pd.DataFrame()  # Start with empty
                    if not all_mutations_for_viz.empty:
                        pair_viz_mutations = all_mutations_for_viz[
                            (
                                all_mutations_for_viz["pair_transcript_id"]
                                == transcript_id
                            )
                            & (all_mutations_for_viz["pair_feature_id"] == feature_id)
                        ].copy()

                    logger.debug(
                        f"Creating view with {len(pair_viz_mutations)} validated mutations"
                    )

                    if "impact_validated" in pair_viz_mutations.columns:
                        logger.debug(
                            f"Validated impacts: {pair_viz_mutations['impact_validated'].unique()}"
                        )

                    # Always call visualizer, even with empty mutations DataFrame
                    visualizer.visualize_transcript(
                        gene_name=gene_name,
                        transcript_id=transcript_id,
                        alt_features=feature_for_viz,
                        mutations_df=pair_viz_mutations,  # Can be empty
                        output_file=str(
                            transcript_dir / f"{pair_base_filename}_filtered.pdf"
                        ),
                    )

                    logger.debug(f"Creating zoomed view")
                    visualizer.visualize_transcript_zoomed(
                        gene_name=gene_name,
                        transcript_id=transcript_id,
                        alt_features=feature_for_viz,
                        mutations_df=pair_viz_mutations,  # Can be empty
                        output_file=str(
                            transcript_dir / f"{pair_base_filename}_filtered_zoom.pdf"
                        ),
                        padding=100,
                    )

            # Calculate total mutations across all pairs
            total_mutations = (
                sum(pair["mutation_count_total"] for pair in pair_results)
                if pair_results
                else 0
            )

            logger.debug("Processing complete")

            # Determine status based on errors and results
            status = "success"
            error_msg = None
            if gene_errors:
                status = "success_with_errors"
                error_msg = "; ".join(gene_errors[:3])  # Limit to first 3 errors
                if len(gene_errors) > 3:
                    error_msg += f" (+{len(gene_errors) - 3} more)"
            elif skipped_pairs > 0 and len(pair_results) == 0:
                status = "all_pairs_failed"
                error_msg = f"{skipped_pairs} pairs skipped"
            elif gene_warnings and len(pair_results) > 0 and total_mutations == 0:
                # Has pairs but no mutations - might be legitimate or API failure
                if any("no_mutations" in w for w in gene_warnings):
                    error_msg = "no_mutations_found"

            return {
                "gene_name": gene_name,
                "status": status,
                "total_transcripts": len(transcript_info),
                "alternative_features": len(alt_features),
                "transcript_feature_pairs": len(transcript_feature_pairs),
                "pairs_processed": len(pair_results),
                "pairs_skipped": skipped_pairs,
                "mutations_filtered": total_mutations,
                "mutations_in_alt_start_sites": sum(
                    pair.get("mutations_in_alt_start_site", 0) for pair in pair_results
                ),
                "gnomad_allele_count": sum(
                    pair.get("gnomad_allele_count", 0) for pair in pair_results
                ),
                "cosmic_sample_count": sum(
                    pair.get("cosmic_sample_count", 0) for pair in pair_results
                ),
                "pair_results": pair_results,
                "consequences_validated": validate_consequences,
                "error": error_msg,
                "warnings": "; ".join(gene_warnings[:5]) if gene_warnings else None,
            }

        except Exception as e:
            logger.error(f"Error processing gene {gene_name}: {str(e)}")
            logger.debug(f"Error: {str(e)}")
            return {"gene_name": gene_name, "status": "error", "error": str(e)}
