"""Mutation data handling and analysis module.

This module provides the MutationHandler class for fetching, processing, and
analyzing mutation data from various sources including gnomAD and ClinVar.
"""

from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport
import pandas as pd
from typing import Dict, Optional, List, Any
import requests
import json
import asyncio
import logging
from datetime import datetime

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
            "retmax": 500,  # Adjust as needed
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
            sources = ["gnomad", "clinvar"]
        sources = [s.lower() for s in sources]

        # Fetch and combine data from sources
        dfs_to_combine = []

        if verbose:
            print(f"Fetching mutations from sources: {', '.join(sources)}...")

        for source in sources:
            try:
                if source == "gnomad":
                    df = await self.get_gnomad_variants(gene_name)
                elif source == "clinvar":
                    df = await self.get_clinvar_variants(gene_name)
                else:
                    continue

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

        # Remove duplicates
        combined_df = combined_df.drop_duplicates(
            subset=["position", "reference", "alternate"], keep="first"
        )

        if verbose:
            print(f"Found {len(combined_df)} unique variants after deduplication")

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


async def analyze_mutations(
    gene_name: str,
    mutation_handler: "MutationHandler",
    alt_features: pd.DataFrame,
    sources: Optional[List[str]] = None,
    impact_types: Optional[Dict[str, List[str]]] = None,
    aggregator_csv_path: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """Analyze mutations in alternative isoform regions.

    Args:
        gene_name: Name of the gene to analyze
        mutation_handler: Initialized MutationHandler instance
        alt_features: DataFrame containing alternative isoform features
        sources: List of mutation sources to query ('clinvar', 'gnomad', 'aggregator')
        impact_types: Dictionary mapping sources to impact types to filter by
        aggregator_csv_path: Path to aggregator CSV file (only needed if 'aggregator' in sources)

    Returns:
        DataFrame containing mutations in alternative isoform regions or None if no mutations found
    """
    if sources is None:
        sources = ["clinvar"]

    if impact_types is None:
        impact_types = {}  # Default to no filtering by impact type

    print(f"Fetching mutations from sources: {', '.join(sources)}...")

    mutations = await mutation_handler.get_visualization_ready_mutations(
        gene_name=gene_name,
        alt_features=alt_features,
        sources=sources,
        aggregator_csv_path=aggregator_csv_path,
    )

    # Filter mutations by impact type for each source if specified
    if not mutations.empty and impact_types:
        print(f"Filtering for impact types by source:")
        filtered_mutations = pd.DataFrame()

        for source, impacts in impact_types.items():
            print(f"  - {source}: {', '.join(impacts)}")
            source_mutations = mutations[
                mutations["source"].str.lower() == source.lower()
            ]

            if not source_mutations.empty:
                filtered_source = source_mutations[
                    source_mutations["impact"].isin(impacts)
                ]
                filtered_mutations = pd.concat([filtered_mutations, filtered_source])

            # Keep mutations from sources that don't have filters specified
            other_sources = [s for s in sources if s.lower() != source.lower()]
            for other_source in other_sources:
                other_mutations = mutations[
                    mutations["source"].str.lower() == other_source.lower()
                ]
                filtered_mutations = pd.concat([filtered_mutations, other_mutations])

        mutations = filtered_mutations

    if mutations.empty:
        print("No matching mutations found")
        return None

    print(f"Found {len(mutations)} mutations in truncation regions")

    # Get mutation statistics
    mutation_impacts = mutations["impact"].value_counts().to_dict()
    clinical_sig = mutations["clinical_significance"].value_counts().to_dict()

    # Create truncation regions string
    truncation_regions = alt_features.apply(
        lambda x: f"{x['start']}-{x['end']}", axis=1
    ).tolist()

    print("\nMutation Analysis:")
    print(f"Impact types: {mutation_impacts}")
    print(f"Clinical significance: {clinical_sig}")
    print(f"Truncation regions: {truncation_regions}")

    return mutations
