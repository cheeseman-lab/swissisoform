#!/usr/bin/env python3

"""Analyze truncation sites and mutations in gene transcripts.

This script batch processes a list of genes to identify alternative isoform
truncation sites and analyze mutations that occur within these regions.
It can generate visualizations of transcript features and mutations.
"""

import asyncio
import pandas as pd
from datetime import datetime
from pathlib import Path
import argparse
import logging
import warnings
from typing import Optional

# Suppress pandas FutureWarning
warnings.filterwarnings("ignore", category=FutureWarning)

from swissisoform.genome import GenomeHandler
from swissisoform.visualize import GenomeVisualizer
from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.mutations import MutationHandler
from swissisoform.utils import (
    parse_gene_list,
    save_gene_level_results,
    save_truncation_level_results,
    print_analysis_summary,
    load_preferred_transcripts,  # New utility function to add
)

# Configure logger
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


async def process_gene(
    gene_name: str,
    genome: GenomeHandler,
    alt_isoforms: AlternativeIsoform,
    mutation_handler: MutationHandler,
    output_dir: str,
    visualize: bool = False,
    include_unfiltered: bool = False,
    impact_types: dict = None,
    preferred_transcripts: set = None,  # New parameter
) -> dict:
    """Process a single gene with visualizations and mutation analysis.

    Args:
        gene_name: Name of the gene to process
        genome: Initialized GenomeHandler instance
        alt_isoforms: Initialized AlternativeIsoform instance
        mutation_handler: Initialized MutationHandler instance
        output_dir: Directory to save output files
        visualize: Whether to generate visualizations
        include_unfiltered: Whether to include unfiltered mutation analysis
        impact_types: Optional dict of mutation impact types to filter by source
        preferred_transcripts: Optional set of transcript IDs to prioritize
    Returns:
        Dictionary containing analysis results
    """
    try:
        # Get alternative isoform features
        print(f"  ├─ Getting alternative features...", end="", flush=True)
        alt_features = alt_isoforms.get_visualization_features(gene_name)
        if alt_features.empty:
            print(f"\r  ├─ No alternative features found")
            return {"gene_name": gene_name, "status": "no_features", "error": None}
        print(f"\r  ├─ Found {len(alt_features)} alternative features")

        # Get transcript information
        print(f"  ├─ Getting transcript information...", end="", flush=True)
        transcript_info = genome.get_transcript_ids(gene_name)
        if transcript_info.empty:
            print(f"\r  ├─ No transcript info found")
            return {"gene_name": gene_name, "status": "no_transcripts", "error": None}

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
                if not (transcript_end < trunc_start or transcript_start > trunc_end):
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
        all_mutations = await mutation_handler.get_visualization_ready_mutations(
            gene_name=gene_name,
            alt_features=None,  # Don't filter by alt_features yet
            sources=["clinvar"],
            aggregator_csv_path=None,
        )

        if all_mutations is None or all_mutations.empty:
            print(f"\r  ├─ No mutations found for this gene")
            all_mutations = pd.DataFrame()
        else:
            print(f"\r  ├─ Found {len(all_mutations)} total mutations for this gene")

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
                        all_mutations = pd.concat([all_mutations, filtered_mutations])

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
                        category_key = f"mutations_{impact.replace(' ', '_').lower()}"
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
                                str(id).strip() for id in impact_ids if str(id).strip()
                            ]

                            # Add to the category-specific IDs
                            if impact_ids:
                                impact_key = (
                                    f"clinvar_ids_{impact.replace(' ', '_').lower()}"
                                )
                                mutation_categories[impact_key] = ",".join(impact_ids)

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
            visualizer = GenomeVisualizer(genome)
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

async def main(
    gene_list_path: str,
    output_dir: str,
    genome_path: str,
    annotation_path: str,
    bed_path: str,
    preferred_transcripts_path: Optional[str] = None,
    visualize: bool = False,
    include_unfiltered: bool = False,
):
    """Main function to process genes for mutation analysis.

    Args:
        gene_list_path: Path to file containing gene names
        output_dir: Directory to save output files
        genome_path: Path to the genome FASTA file
        annotation_path: Path to the genome annotation GTF file
        bed_path: Path to the alternative isoform BED file
        preferred_transcripts_path: Optional path to file with preferred transcript IDs
        visualize: Whether to generate visualizations
        include_unfiltered: Whether to include unfiltered mutation analysis
    """
    start_time = datetime.now()
    print(f"Starting analysis at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    # Create output directory
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # Initialize handlers
    print("\nInitializing components:")
    print("  ├─ Loading genome data...")
    genome = GenomeHandler(genome_path, annotation_path)
    print("  ├─ Loading alternative isoform data...")
    alt_isoforms = AlternativeIsoform()
    alt_isoforms.load_bed(bed_path)
    print("  └─ Initializing mutation handler...")
    mutation_handler = MutationHandler()

    # Load preferred transcripts if provided
    preferred_transcripts = None
    if preferred_transcripts_path:
        print(f"\nLoading preferred transcripts from {preferred_transcripts_path}")
        preferred_transcripts = load_preferred_transcripts(preferred_transcripts_path)
        print(f"Loaded {len(preferred_transcripts)} preferred transcript IDs")

    # Define impact types for filtering
    impact_types = {
        "clinvar": ["missense variant", "nonsense variant", "frameshift variant"]
    }

    # Read gene list
    print(f"\nReading gene list from {gene_list_path}")
    gene_names = parse_gene_list(gene_list_path)

    total_genes = len(gene_names)
    print(f"\nStarting analysis of {total_genes} genes")
    print(f"Impact types filter: {impact_types['clinvar']}")
    print(f"Analyzing transcript-truncation pairs with mutations in truncation regions")
    if preferred_transcripts:
        print(
            f"Using {len(preferred_transcripts)} preferred transcripts when available"
        )

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
            impact_types=impact_types,
            preferred_transcripts=preferred_transcripts,
        )
        results.append(result)

        # Save both levels of results
        save_gene_level_results(results, output_dir)
        save_truncation_level_results(results, output_dir)

    # Final summary
    end_time = datetime.now()
    duration = end_time - start_time

    # Create and print summary
    results_df = pd.DataFrame(results)
    print_analysis_summary(results_df, output_dir)

    print(
        f"  └─ Analysis completed in {duration} at {end_time.strftime('%Y-%m-%d %H:%M:%S')}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Batch process genes for mutation analysis"
    )
    parser.add_argument("gene_list", help="Path to file containing gene names")
    parser.add_argument("output_dir", help="Directory to save output files")
    parser.add_argument(
        "--genome",
        default="../data/genome_data/GRCh38.p7.genome.fa",
        help="Path to genome FASTA",
    )
    parser.add_argument(
        "--annotation",
        default="../data/genome_data/gencode.v25.annotation.ensembl_cleaned.gtf",
        help="Path to genome annotation",
    )
    parser.add_argument(
        "--bed",
        default="../data/ribosome_profiling/full_truncations_JL_cleaned.bed",
        help="Path to alternative isoform BED file",
    )
    parser.add_argument(
        "--preferred-transcripts",
        default="../data/genome_data/hela_top_transcript.txt",
        help="Path to file containing preferred transcript IDs",
    )
    parser.add_argument(
        "--visualize", action="store_true", help="Generate visualizations"
    )
    parser.add_argument(
        "--include-unfiltered",
        action="store_true",
        help="Include unfiltered mutation analysis",
    )

    args = parser.parse_args()

    asyncio.run(
        main(
            gene_list_path=args.gene_list,
            output_dir=args.output_dir,
            genome_path=args.genome,
            annotation_path=args.annotation,
            bed_path=args.bed,
            preferred_transcripts_path=args.preferred_transcripts,
            visualize=args.visualize,
            include_unfiltered=args.include_unfiltered,
        )
    )
