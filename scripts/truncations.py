#!/usr/bin/env python3

import asyncio
import pandas as pd
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Dict
import os
import warnings

# Suppress pandas FutureWarning
warnings.filterwarnings('ignore', category=FutureWarning)

from swissisoform.genome import GenomeHandler
from swissisoform.visualize import GenomeVisualizer
from swissisoform.isoform import AlternativeIsoform
from swissisoform.mutations import MutationHandler
from swissisoform.utils import analyze_mutations

async def process_gene(
    gene_name: str,
    genome: GenomeHandler,
    alt_isoforms: AlternativeIsoform,
    mutation_handler: MutationHandler,
    output_dir: str,
    visualize: bool = False,
    include_unfiltered: bool = False,
    impact_types: Optional[Dict[str, List[str]]] = None
) -> dict:
    """
    Process a single gene with visualizations and mutation analysis.
    
    Args:
        gene_name: Name of the gene to process
        genome: Initialized GenomeHandler instance
        alt_isoforms: Initialized AlternativeIsoform instance
        mutation_handler: Initialized MutationHandler instance
        output_dir: Directory to save output files
        visualize: Whether to generate visualizations
        include_unfiltered: Whether to include unfiltered mutation analysis
        impact_types: Optional dict of mutation impact types to filter by source
        
    Returns:
        Dictionary containing analysis results
    """
    try:
        # Get alternative isoform features
        print(f"  ├─ Getting alternative features...", end='', flush=True)
        alt_features = alt_isoforms.get_visualization_features(gene_name)
        
        if alt_features.empty:
            print(f"\r  ├─ No alternative features found")
            return {
                'gene_name': gene_name,
                'status': 'no_features',
                'error': None
            }
        
        print(f"\r  ├─ Found {len(alt_features)} alternative features")
        
        # Get transcript information
        print(f"  ├─ Getting transcript information...", end='', flush=True)
        transcript_info = genome.get_transcript_ids(gene_name)
        
        if transcript_info.empty:
            print(f"\r  ├─ No transcript info found")
            return {
                'gene_name': gene_name,
                'status': 'no_transcripts',
                'error': None
            }
        
        print(f"\r  ├─ Found {len(transcript_info)} transcripts")
        
        # Get mutations (unfiltered) if requested
        mutations_unfiltered = None
        unfiltered_count = 0
        if include_unfiltered:
            print(f"  ├─ Fetching unfiltered mutations...", end='', flush=True)
            mutations_unfiltered = await analyze_mutations(
                gene_name=gene_name,
                mutation_handler=mutation_handler,
                alt_features=alt_features,
                sources=["clinvar"]
            )
            unfiltered_count = len(mutations_unfiltered) if mutations_unfiltered is not None else 0
            print(f"\r  ├─ Found {unfiltered_count} unfiltered mutations")
        
        # Get mutations (filtered)
        print(f"  ├─ Analyzing mutations...", end='', flush=True)
        mutations_filtered = await analyze_mutations(
            gene_name=gene_name,
            mutation_handler=mutation_handler,
            alt_features=alt_features,
            sources=["clinvar"],
            impact_types=impact_types
        )
        
        filtered_count = len(mutations_filtered) if mutations_filtered is not None else 0
        print(f"\r  ├─ Found {filtered_count} mutations after filtering")
        
        if visualize:
            visualizer = GenomeVisualizer(genome)
            gene_dir = Path(output_dir) / gene_name
            gene_dir.mkdir(parents=True, exist_ok=True)
            
            print(f"  ├─ Generating visualizations:")
            for _, transcript in transcript_info.iterrows():
                transcript_id = transcript['transcript_id']
                print(f"  │  ├─ Processing {transcript_id}")
                
                if include_unfiltered and mutations_unfiltered is not None:
                    print(f"  │  │  ├─ Creating unfiltered view")
                    visualizer.visualize_transcript(
                        gene_name=gene_name,
                        transcript_id=transcript_id,
                        alt_features=alt_features,
                        mutations_df=mutations_unfiltered,
                        output_file=str(gene_dir / f'{transcript_id}_unfiltered.png')
                    )
                
                if mutations_filtered is not None:
                    print(f"  │  │  ├─ Creating filtered view")
                    visualizer.visualize_transcript(
                        gene_name=gene_name,
                        transcript_id=transcript_id,
                        alt_features=alt_features,
                        mutations_df=mutations_filtered,
                        output_file=str(gene_dir / f'{transcript_id}_filtered.png')
                    )
                    
                    print(f"  │  │  └─ Creating zoomed view")
                    visualizer.visualize_transcript_zoomed(
                        gene_name=gene_name,
                        transcript_id=transcript_id,
                        alt_features=alt_features,
                        mutations_df=mutations_filtered,
                        output_file=str(gene_dir / f'{transcript_id}_filtered_zoom.png'),
                        padding=100
                    )
                print(f"  │  └─ Completed {transcript_id}")
        
        print("  └─ Processing complete")
        return {
            'gene_name': gene_name,
            'status': 'success',
            'transcripts': len(transcript_info),
            'alt_features': len(alt_features),
            'mutations_unfiltered': unfiltered_count if include_unfiltered else None,
            'mutations_filtered': filtered_count,
            'error': None
        }
        
    except Exception as e:
        print(f"  └─ Error: {str(e)}")
        return {
            'gene_name': gene_name,
            'status': 'error',
            'error': str(e)
        }

async def main(
    gene_list_path: str, 
    output_dir: str,
    visualize: bool = False,
    include_unfiltered: bool = False
):
    print(f"Starting analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Initialize handlers
    print("\nInitializing components:")
    print("  ├─ Loading genome data...")
    genome = GenomeHandler(
        '../data/genome_data/hg38.fa',
        '../data/genome_data/hg38.ncbiRefSeq.gtf'
    )
    
    print("  ├─ Loading alternative isoform data...")
    alt_isoforms = AlternativeIsoform()
    alt_isoforms.load_bed('../data/ribosome_profiling/RiboTISHV6_MD2025_AnnoToTruncation_exonintersect.bed')
    
    print("  └─ Initializing mutation handler...")
    mutation_handler = MutationHandler()
    
    # Define impact types for filtering
    impact_types = {
        "clinvar": ["missense variant", "nonsense variant", "frameshift variant"]
    }
    
    # Read gene list
    print(f"\nReading gene list from {gene_list_path}")
    with open(gene_list_path, 'r') as f:
        gene_names = [line.strip() for line in f if line.strip()]
    
    total_genes = len(gene_names)
    print(f"\nStarting analysis of {total_genes} genes")
    print(f"Impact types filter: {impact_types['clinvar']}")
    
    # Process all genes
    results = []
    for idx, gene_name in enumerate(gene_names, 1):
        print(f"\nProcessing gene {idx}/{total_genes}: {gene_name}")
        result = await process_gene(
            gene_name=gene_name,
            genome=genome,
            alt_isoforms=alt_isoforms,
            mutation_handler=mutation_handler,
            output_dir=output_dir,
            visualize=visualize,
            include_unfiltered=include_unfiltered,
            impact_types=impact_types
        )
        results.append(result)
        
        # Save intermediate results
        results_df = pd.DataFrame(results)
        results_df.to_csv(Path(output_dir) / 'analysis_results.csv', index=False)
    
    # Final summary
    print("\nAnalysis Summary:")
    results_df = pd.DataFrame(results)
    print(f"  ├─ Total genes processed: {len(results_df)}")
    print("\n  ├─ Status breakdown:")
    for status, count in results_df['status'].value_counts().items():
        print(f"  │  ├─ {status}: {count}")
    
    # Genes with errors
    error_genes = results_df[results_df['status'] == 'error']
    if not error_genes.empty:
        print("\n  ├─ Genes with errors:")
        for _, row in error_genes.iterrows():
            print(f"  │  ├─ {row['gene_name']}: {row['error']}")
    
    print(f"\n  ├─ Results saved to: {output_dir}")
    print(f"  └─ Analysis completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Batch process genes for mutation analysis')
    parser.add_argument('gene_list', help='Path to file containing gene names')
    parser.add_argument('output_dir', help='Directory to save output files')
    parser.add_argument('--visualize', action='store_true', help='Generate visualizations')
    parser.add_argument('--include-unfiltered', action='store_true', help='Include unfiltered mutation analysis')
    
    args = parser.parse_args()
    
    asyncio.run(main(
        gene_list_path=args.gene_list,
        output_dir=args.output_dir,
        visualize=args.visualize,
        include_unfiltered=args.include_unfiltered
    ))