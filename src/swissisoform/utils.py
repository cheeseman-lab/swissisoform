"""Utility functions for the SwissIsoform package.

This module provides various helper functions for mutation analysis,
file handling, and data processing used throughout the package.
"""
import pandas as pd
from typing import List, Dict, Optional, Union
import asyncio
import logging
from pathlib import Path
from collections import defaultdict
from biomart import BiomartServer


# Configure logger
logger = logging.getLogger(__name__)

def save_gene_level_results(
    results: List[Dict],
    output_dir: Union[str, Path],
    filename: str = "gene_level_results.csv",
) -> None:
    """Save gene-level analysis results to CSV file.

    Args:
        results: List of result dictionaries to save
        output_dir: Directory to save the results file
        filename: Name of the results file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Extract gene-level results, removing pair-specific data
    gene_results = []
    for result in results:
        # Skip pair-specific details to keep the summary concise
        if "pair_results" in result:
            result = {k: v for k, v in result.items() if k != "pair_results"}
        gene_results.append(result)
    # Save to CSV
    gene_df = pd.DataFrame(gene_results)
    output_path = output_dir / filename
    gene_df.to_csv(output_path, index=False)
    print(f"Gene-level results saved to {output_path}")

def save_truncation_level_results(
    results: List[Dict],
    output_dir: Union[str, Path],
    filename: str = "truncation_level_results.csv",
) -> None:
    """Save detailed transcript-truncation pair analysis results to CSV file.

    Args:
        results: List of result dictionaries containing pair_results
        output_dir: Directory to save the results file
        filename: Name of the results file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Extract all pair results
    all_pairs = []
    for result in results:
        if "pair_results" in result and result["pair_results"]:
            gene_name = result["gene_name"]
            # Add gene name to each pair
            for pair in result["pair_results"]:
                pair_with_gene = {"gene_name": gene_name, **pair}
                all_pairs.append(pair_with_gene)
    # If no pairs found, return early
    if not all_pairs:
        print(f"No transcript-truncation pairs found to save")
        return
    # Convert to DataFrame
    pairs_df = pd.DataFrame(all_pairs)
    # Organize columns in preferred order
    column_order = [
        "gene_name",              # Gene name first
        "transcript_id",          # Transcript ID second
        "truncation_id",          # Alternative feature ID third
        "truncation_start", 
        "truncation_end",
        "mutation_count_total"    # Total mutation count
    ]
    
    # Add mutation category columns to the order (they start with "mutations_")
    mutation_category_columns = [col for col in pairs_df.columns if col.startswith("mutations_")]
    column_order.extend(mutation_category_columns)
    
    # Add the ClinVar IDs column 
    clinvar_columns = [col for col in pairs_df.columns if col.startswith("clinvar")]
    column_order.extend(clinvar_columns)
    
    # Reorder columns that exist in the dataframe
    available_columns = [col for col in column_order if col in pairs_df.columns]
    
    # Add any remaining columns that weren't explicitly ordered
    remaining_columns = [col for col in pairs_df.columns if col not in available_columns]
    final_column_order = available_columns + remaining_columns
    
    # Apply the column order and save
    if final_column_order:  # Only reorder if we have columns to reorder
        pairs_df = pairs_df[final_column_order]
        
    # Save to CSV
    output_path = output_dir / filename
    pairs_df.to_csv(output_path, index=False)
    print(f"Transcript-truncation pair analysis saved to {output_path}")

def print_analysis_summary(results_df, output_dir):
    """Print summary statistics from the analysis results.
    
    Args:
        results_df: DataFrame containing analysis results
        output_dir: Directory where output files are saved
    """
    print("\nAnalysis Summary:")
    print(f"  ├─ Total genes processed: {len(results_df)}")
    print("\n  ├─ Status breakdown:")
    for status, count in results_df["status"].value_counts().items():
        print(f"  │  ├─ {status}: {count}")
    
    # Transcript-truncation statistics
    successful_genes = results_df[results_df["status"] == "success"]
    if not successful_genes.empty:
        total_transcripts = successful_genes["total_transcripts"].sum()
        total_truncations = successful_genes["truncation_features"].sum()
        total_pairs = successful_genes["transcript_truncation_pairs"].sum()
        
        # Calculate average pairs per gene
        avg_pairs_per_gene = total_pairs / len(successful_genes) if len(successful_genes) > 0 else 0
        
        print(f"\n  ├─ Transcript-Truncation Analysis:")
        print(f"  │  ├─ Total transcripts across all genes: {total_transcripts}")
        print(f"  │  ├─ Total truncation features: {total_truncations}")
        print(f"  │  ├─ Total transcript-truncation pairs: {total_pairs}")
        print(f"  │  └─ Average pairs per gene: {avg_pairs_per_gene:.2f}")
        
        # Mutation statistics if available
        if "mutations_filtered" in successful_genes.columns:
            total_mutations = successful_genes["mutations_filtered"].sum()
            print(f"\n  ├─ Mutation Analysis:")
            print(f"  │  ├─ Total mutations in truncation regions: {total_mutations}")
            
            # Try to load mutation analysis results to report statistics
            mutation_analysis_path = Path(output_dir) / "mutation_analysis_by_pair.csv"
            if mutation_analysis_path.exists():
                try:
                    pairs_df = pd.read_csv(mutation_analysis_path)
                    
                    # Get total mutations across all pairs
                    total_pair_mutations = pairs_df["mutation_count_total"].sum() if "mutation_count_total" in pairs_df.columns else 0
                    avg_mutations_per_pair = total_pair_mutations / len(pairs_df) if len(pairs_df) > 0 else 0
                    print(f"  │  ├─ Total mutations across all pairs: {total_pair_mutations}")
                    print(f"  │  ├─ Average mutations per transcript-truncation pair: {avg_mutations_per_pair:.2f}")
                    
                    # Print statistics for each mutation category
                    mutation_categories = [col for col in pairs_df.columns if col.startswith("mutations_")]
                    if mutation_categories:
                        print(f"  │  │")
                        print(f"  │  ├─ Breakdown by mutation category:")
                        
                        for category in mutation_categories:
                            # Skip if the column doesn't exist in the dataframe
                            if category not in pairs_df.columns:
                                continue
                                
                            # Convert category name from mutations_missense_variant to "Missense Variant"
                            category_name = category.replace("mutations_", "").replace("_", " ")
                            category_name = category_name.title()
                            
                            # Calculate statistics for this category
                            category_total = pairs_df[category].sum()
                            category_percent = (category_total / total_pair_mutations * 100) if total_pair_mutations > 0 else 0
                            
                            print(f"  │  │  ├─ {category_name}: {category_total} ({category_percent:.1f}%)")
                    
                    print(f"  │  └─ Detailed results available in mutation_analysis_by_pair.csv")
                except Exception as e:
                    logger.error(f"Error reading mutation analysis: {str(e)}")
                    print(f"  │  └─ Error reading detailed mutation analysis: {str(e)}")
    # Genes with errors
    error_genes = results_df[results_df["status"] == "error"]
    if not error_genes.empty:
        print("\n  ├─ Genes with errors:")
        for _, row in error_genes.iterrows():
            print(f"  │  ├─ {row['gene_name']}: {row['error']}")
    print(f"\n  ├─ Results saved to: {output_dir}")
    print(f"  ├─ Gene-level results saved to: {Path(output_dir) / 'gene_level_results.csv'}")
    print(f"  ├─ Detailed mutation analysis by pair saved to: {Path(output_dir) / 'mutation_analysis_by_pair.csv'}")

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


def get_ensembl_reference():
    """Fetch current ENSEMBL ID to gene name mappings using biomart."""
    server = BiomartServer("http://www.ensembl.org/biomart")
    database = server.databases["ENSEMBL_MART_ENSEMBL"]
    dataset = database.datasets["hsapiens_gene_ensembl"]
    attributes = ["ensembl_gene_id", "external_gene_name"]
    response = dataset.search({"attributes": attributes})
    reference_data = {}
    for line in response.iter_lines():
        ensembl_id, gene_name = line.decode("utf-8").split("\t")
        if ensembl_id and gene_name:
            reference_data[ensembl_id] = gene_name
    return reference_data


def parse_bed_line(line):
    """Parse a single BED file line to extract ENSEMBL ID and gene name.
    
    Args:
        line: BED file line to parse
    Returns:
        Dictionary containing parsed data or None if invalid
    """
    fields = line.strip().split("\t")
    if len(fields) < 4: return None
    name_parts = fields[3].split("_")
    if len(name_parts) < 3: return None
    full_ensembl_id = name_parts[0]
    ensembl_id = full_ensembl_id.split(".")[0]
    gene_name = name_parts[1]
    return {
        "ensembl_id": ensembl_id,
        "full_ensembl_id": full_ensembl_id,
        "gene_name": gene_name,
        "fields": fields,
        "line": line.strip(),
    }

def cleanup_bed(input_bed_path, output_bed_path, verbose=True):
    """Clean BED file: fix gene names, remove invalid entries and duplicates.
    
    Args:
        input_bed_path: Path to input BED file
        output_bed_path: Path to save cleaned BED file
        verbose: Print detailed cleanup summary
    Returns:
        Dictionary containing cleanup statistics
    """
    input_bed_path, output_bed_path = Path(input_bed_path), Path(output_bed_path)
    
    # Get reference data
    if verbose: print("Fetching Ensembl reference data...")
    reference_data = get_ensembl_reference()
    if verbose: print(f"Retrieved {len(reference_data)} mappings")
    
    # Process file
    valid_entries, seen_positions = [], set()
    stats = {"total": 0, "invalid_format": 0, "invalid_ensembl": 0, 
             "duplicates": 0, "updated": 0, "valid": 0}
    
    with open(input_bed_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            if not line.strip(): continue
            stats["total"] += 1
            
            # Validate format and Ensembl ID
            parsed = parse_bed_line(line)
            if not parsed or not parsed["ensembl_id"].startswith("ENSG"):
                stats["invalid_format"] += 1
                continue
            if parsed["ensembl_id"] not in reference_data:
                stats["invalid_ensembl"] += 1
                continue
                
            # Update gene name if needed
            ref_gene_name = reference_data[parsed["ensembl_id"]]
            fields = parsed["fields"].copy()
            name_parts = fields[3].split("_")
            
            if name_parts[1].upper() != ref_gene_name.upper():
                name_parts[1] = ref_gene_name
                fields[3] = "_".join(name_parts)
                stats["updated"] += 1
                
            # Check for duplicates
            position_key = f"{parsed['ensembl_id']}_{fields[0]}_{fields[1]}_{fields[2]}"
            if position_key in seen_positions:
                stats["duplicates"] += 1
                continue
                
            # Add valid entry
            seen_positions.add(position_key)
            valid_entries.append("\t".join(fields))
            stats["valid"] += 1
    
    # Write output
    output_bed_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_bed_path, "w") as f_out:
        f_out.write("\n".join(valid_entries))
    
    if verbose:
        print(f"\nCleanup Summary:")
        print(f"  Total entries: {stats['total']}")
        print(f"  Invalid entries removed: {stats['invalid_format'] + stats['invalid_ensembl']}")
        print(f"  Duplicates removed: {stats['duplicates']}")
        print(f"  Gene names updated: {stats['updated']}")
        print(f"  Valid entries in final file: {stats['valid']}")
    
    return stats
