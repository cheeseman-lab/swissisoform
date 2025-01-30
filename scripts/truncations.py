import pandas as pd
import asyncio
from typing import List, Dict
from swissisoform.genome import GenomeHandler
from swissisoform.isoform import AlternativeIsoform
from swissisoform.mutations import MutationHandler
from datetime import datetime

async def analyze_gene_truncation_mutations(
    gene_names: List[str],
    genome: GenomeHandler,
    alt_isoforms: AlternativeIsoform,
    mutation_handler: MutationHandler
) -> pd.DataFrame:
    """
    Analyze mutations in truncation regions for a list of genes.
    
    Args:
        gene_names: List of gene names to analyze
        genome: Initialized GenomeHandler instance
        alt_isoforms: Initialized AlternativeIsoform instance
        mutation_handler: Initialized MutationHandler instance
        
    Returns:
        DataFrame containing analysis results for each gene
    """
    results = []
    total_genes = len(gene_names)
    
    print(f"\nStarting analysis of {total_genes} genes at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    for idx, gene_name in enumerate(gene_names, 1):
        try:
            print(f"\nProcessing gene {idx}/{total_genes}: {gene_name}")
            
            # Get alternative isoform features
            print(f"  ├─ Getting alternative features...", end='', flush=True)
            alt_features = alt_isoforms.get_visualization_features(gene_name)
            
            if alt_features.empty:
                print(f"\r  ├─ No alternative features found")
                continue
            print(f"\r  ├─ Found {len(alt_features)} alternative features")
                
            # Get transcript info
            print(f"  ├─ Getting transcript information...", end='', flush=True)
            transcript_info = genome.get_transcript_ids(gene_name)
            if transcript_info.empty:
                print(f"\r  ├─ No transcript info found")
                continue
            print(f"\r  ├─ Found {len(transcript_info)} transcripts")
                
            # Get ClinVar mutations
            print(f"  ├─ Fetching ClinVar mutations...", end='', flush=True)
            mutations = await mutation_handler.get_visualization_ready_mutations(
                gene_name=gene_name,
                alt_features=alt_features,
                sources=['clinvar']
            )
            
            if mutations.empty:
                print(f"\r  └─ No ClinVar mutations found")
                results.append({
                    'gene_name': gene_name,
                    'has_truncation_region': True,
                    'truncation_regions': len(alt_features),
                    'mutations_in_region': 0,
                    'mutation_impacts': None,
                    'clinical_significance': None,
                    'truncation_coordinates': alt_features.apply(
                        lambda x: f"{x['start']}-{x['end']}", axis=1).tolist()
                })
                continue
            
            # Analyze mutations
            print(f"\r  └─ Found {len(mutations)} mutations in truncation regions")
            mutation_impacts = mutations['impact'].value_counts().to_dict()
            clinical_sig = mutations['clinical_significance'].value_counts().to_dict()
            
            results.append({
                'gene_name': gene_name,
                'has_truncation_region': True,
                'truncation_regions': len(alt_features),
                'mutations_in_region': len(mutations),
                'mutation_impacts': mutation_impacts,
                'clinical_significance': clinical_sig,
                'truncation_coordinates': alt_features.apply(
                    lambda x: f"{x['start']}-{x['end']}", axis=1).tolist()
            })
            
        except Exception as e:
            print(f"\r  └─ Error: {str(e)}")
            continue
    
    print(f"\nAnalysis completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    return pd.DataFrame(results)

async def main(gene_list_path: str):
    print(f"Initializing analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Initialize handlers
    print("Loading genome data...")
    genome = GenomeHandler(
        '../data/genome_data/hg38.fa',
        '../data/genome_data/hg38.ncbiRefSeq.gtf'
    )
    
    print("Loading alternative isoform data...")
    alt_isoforms = AlternativeIsoform()
    alt_isoforms.load_bed('../data/ribosome_profiling/RiboTISHV6_Ly2024_AnnoToTruncation_exonintersect.bed')
    
    print("Initializing mutation handler...")
    mutation_handler = MutationHandler()
    
    # Read gene list
    print(f"Reading gene list from {gene_list_path}")
    with open(gene_list_path, 'r') as f:
        gene_names = [line.strip() for line in f if line.strip()]
    
    # Run analysis
    results_df = await analyze_gene_truncation_mutations(
        gene_names,
        genome,
        alt_isoforms,
        mutation_handler
    )
    
    # Basic summary statistics
    print("\nAnalysis Summary:")
    print(f"Total genes analyzed: {len(results_df)}")
    
    genes_with_mutations = results_df[results_df['mutations_in_region'] > 0]
    print(f"Genes with mutations in truncation regions: {len(genes_with_mutations)}")
    
    if not genes_with_mutations.empty:
        print("\nBreakdown of genes with mutations:")
        for _, row in genes_with_mutations.iterrows():
            print(f"\n{row['gene_name']}:")
            print(f"  Number of mutations: {row['mutations_in_region']}")
            print(f"  Impact types: {row['mutation_impacts']}")
            print(f"  Clinical significance: {row['clinical_significance']}")
            print(f"  Truncation regions: {row['truncation_coordinates']}")
    
    # Save results
    output_file = 'truncation_mutation_analysis.csv'
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")
    
    return results_df

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python script.py path_to_gene_list.txt")
        sys.exit(1)
    
    asyncio.run(main(sys.argv[1]))