import pandas as pd
from typing import List, Dict, Optional, Union
import asyncio
import logging
from pathlib import Path

# Configure logger
logger = logging.getLogger(__name__)


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


def save_analysis_results(
    results: List[Dict],
    output_dir: Union[str, Path],
    filename: str = "analysis_results.csv",
) -> None:
    """Save analysis results to CSV file.

    Args:
        results: List of result dictionaries to save
        output_dir: Directory to save the results file
        filename: Name of the results file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results_df = pd.DataFrame(results)
    output_path = output_dir / filename

    results_df.to_csv(output_path, index=False)
    print(f"Results saved to {output_path}")


def validate_gene_list(
    gene_list_path: Union[str, Path], genome_handler: "GenomeHandler"
) -> List[str]:
    """Validate a list of genes against available genomic data.

    Args:
        gene_list_path: Path to file containing gene names
        genome_handler: Initialized GenomeHandler instance

    Returns:
        List of valid gene names
    """
    gene_list_path = Path(gene_list_path)

    if not gene_list_path.exists():
        raise FileNotFoundError(f"Gene list file not found: {gene_list_path}")

    # Read gene list
    with open(gene_list_path, "r") as f:
        gene_names = [line.strip() for line in f if line.strip()]

    # Validate genes
    valid_genes = []
    invalid_genes = []

    for gene_name in gene_names:
        features = genome_handler.find_gene_features(gene_name)
        if features.empty:
            invalid_genes.append(gene_name)
        else:
            valid_genes.append(gene_name)

    if invalid_genes:
        print(f"Warning: {len(invalid_genes)} genes not found in genome annotations:")
        print(
            ", ".join(invalid_genes[:10]) + (", ..." if len(invalid_genes) > 10 else "")
        )

    print(f"Found {len(valid_genes)} valid genes")
    return valid_genes


def parse_gene_list(gene_list_path: Union[str, Path]) -> List[str]:
    """Parse a file containing a list of gene names.

    Args:
        gene_list_path: Path to file containing gene names

    Returns:
        List of gene names
    """
    gene_list_path = Path(gene_list_path)

    if not gene_list_path.exists():
        raise FileNotFoundError(f"Gene list file not found: {gene_list_path}")

    with open(gene_list_path, "r") as f:
        gene_names = [line.strip() for line in f if line.strip()]

    print(f"Read {len(gene_names)} genes from {gene_list_path}")
    return gene_names
