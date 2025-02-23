import pandas as pd
from typing import List, Dict, Optional
from swissisoform.mutations import MutationHandler  # For MutationHandler
from swissisoform.isoform import AlternativeIsoform  # If alt_features is of this type
import asyncio  # Since analyze_mutations is an async function

async def analyze_mutations(gene_name: str, 
                            mutation_handler: MutationHandler,
                            alt_features: pd.DataFrame,
                            sources: Optional[List[str]] = None,
                            impact_types: Optional[Dict[str, List[str]]] = None,
                            aggregator_csv_path: Optional[str] = None):
    
    if sources is None:
        sources = ['clinvar']
    
    if impact_types is None:
        impact_types = {}  # Default to an empty dict, meaning no filtering by impact type
    
    print(f"Fetching mutations from sources: {', '.join(sources)}...")
    
    mutations = await mutation_handler.get_visualization_ready_mutations(
        gene_name=gene_name,
        alt_features=alt_features,
        sources=sources,
        aggregator_csv_path=aggregator_csv_path
    )
    
    # Filter mutations by impact type for each source if specified
    if not mutations.empty and impact_types:
        print(f"Filtering for impact types by source:")
        for source, impacts in impact_types.items():
            print(f"  - {source}: {', '.join(impacts)}")
            mutations = mutations[mutations['impact'].isin(impacts)]
    
    if mutations.empty:
        print("No matching mutations found")
        return None
    
    print(f"Found {len(mutations)} mutations in truncation regions")
    
    # Get mutation statistics (as in script)
    mutation_impacts = mutations['impact'].value_counts().to_dict()
    clinical_sig = mutations['clinical_significance'].value_counts().to_dict()
    
    # Create truncation regions string without nested f-string
    truncation_regions = alt_features.apply(
        lambda x: f"{x.start}-{x.end}", 
        axis=1
    ).tolist()
    
    print("\nMutation Analysis:")
    print(f"Impact types: {mutation_impacts}")
    print(f"Clinical significance: {clinical_sig}")
    print(f"Truncation regions: {truncation_regions}")
    
    return mutations