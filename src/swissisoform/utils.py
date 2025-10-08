"""Refactored BED cleanup utilities with improved transcript selection logic.

This module provides cleaner, more systematic functions for:
1. Building comprehensive transcript databases from multiple sources
2. Intelligent transcript selection using preferred lists and MANE Select
3. Enhanced BED file processing with better transcript mapping
"""

import pandas as pd
from typing import List, Dict, Optional, Union, Tuple, Set
import logging
from pathlib import Path
from collections import defaultdict
import re
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class TranscriptInfo:
    """Data class to hold transcript information."""

    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str
    start_pos: int  # TSS position
    is_preferred: bool = False
    is_mane_select: bool = False
    refseq_id: Optional[str] = None


class TranscriptDatabase:
    """Comprehensive transcript database with intelligent selection logic."""

    def __init__(self):
        """Initialize transcript containers, mappings, and transcript sets."""
        self.transcripts: Dict[str, TranscriptInfo] = {}
        self.gene_transcripts: Dict[str, List[str]] = defaultdict(list)
        self.gene_names: Dict[str, str] = {}
        self.preferred_transcripts: Set[str] = set()
        self.mane_select_transcripts: Set[str] = set()
        self.ensembl_to_refseq: Dict[str, str] = {}

    def load_from_gtf(self, gtf_path: Union[str, Path], verbose: bool = True) -> None:
        """Load transcript information from GENCODE GTF."""
        if verbose:
            logger.info(f"Loading transcripts from GTF: {gtf_path}")

        transcript_count = 0

        with open(gtf_path, "r") as f:
            for line in f:
                if line.startswith("#") or "\tstart_codon\t" not in line:
                    continue

                fields = line.strip().split("\t")
                chrom, start, end, strand = (
                    fields[0],
                    int(fields[3]),
                    int(fields[4]),
                    fields[6],
                )

                gene_match = re.search(r'gene_id "([^"]+)"', fields[8])
                transcript_match = re.search(r'transcript_id "([^"]+)"', fields[8])
                gene_name_match = re.search(r'gene_name "([^"]+)"', fields[8])

                if gene_match and transcript_match and gene_name_match:
                    gene_id = gene_match.group(1).split(".")[0]
                    transcript_id = transcript_match.group(1)
                    gene_name = gene_name_match.group(1)
                    start_pos = start if strand == "+" else end

                    transcript_info = TranscriptInfo(
                        transcript_id=transcript_id,
                        gene_id=gene_id,
                        gene_name=gene_name,
                        chrom=chrom,
                        start=start,
                        end=end,
                        strand=strand,
                        start_pos=start_pos,
                    )

                    self.transcripts[transcript_id] = transcript_info
                    self.gene_transcripts[gene_id].append(transcript_id)
                    self.gene_names[gene_id] = gene_name
                    transcript_count += 1

        if verbose:
            logger.info(
                f"Loaded {transcript_count} transcripts for {len(self.gene_transcripts)} genes"
            )

    def load_preferred_transcripts(
        self, preferred_path: Union[str, Path], verbose: bool = True
    ) -> None:
        """Load preferred transcript list."""
        if verbose:
            logger.info(f"Loading preferred transcripts: {preferred_path}")

        with open(preferred_path, "r") as f:
            self.preferred_transcripts = {line.strip() for line in f if line.strip()}

        # Mark preferred transcripts
        preferred_count = 0
        for transcript_id in self.preferred_transcripts:
            if transcript_id in self.transcripts:
                self.transcripts[transcript_id].is_preferred = True
                preferred_count += 1

        if verbose:
            logger.info(
                f"Loaded {len(self.preferred_transcripts)} preferred transcripts"
            )
            logger.info(f"Found {preferred_count} in current dataset")

    def load_refseq_mapping(
        self, refseq_gtf_path: Union[str, Path], verbose: bool = True
    ) -> None:
        """Load RefSeq mapping and MANE Select information from RefSeq GTF."""
        if verbose:
            logger.info(f"Loading RefSeq mappings: {refseq_gtf_path}")

        mane_count = 0
        mapping_count = 0

        with open(refseq_gtf_path, "r") as f:
            for line in f:
                if line.startswith("#") or "\ttranscript\t" not in line:
                    continue

                attributes = line.strip().split("\t")[8]

                # Extract RefSeq transcript ID and Ensembl cross-reference
                refseq_match = re.search(r'transcript_id "([^"]+)"', attributes)
                ensembl_match = re.search(r'db_xref "Ensembl:([^"]+)"', attributes)

                if refseq_match and ensembl_match:
                    refseq_id = refseq_match.group(1)
                    ensembl_id = ensembl_match.group(1)

                    # Only keep RefSeq transcript IDs
                    if refseq_id.startswith(("NM_", "NR_", "XM_", "XR_")):
                        # Store mapping
                        self.ensembl_to_refseq[ensembl_id] = refseq_id
                        base_ensembl = ensembl_id.split(".")[0]
                        self.ensembl_to_refseq[base_ensembl] = refseq_id
                        mapping_count += 1

                        # Check for MANE Select tag
                        if 'tag "MANE Select"' in attributes:
                            self.mane_select_transcripts.add(ensembl_id)
                            self.mane_select_transcripts.add(base_ensembl)

                            # Mark in our transcript database
                            if ensembl_id in self.transcripts:
                                self.transcripts[ensembl_id].is_mane_select = True
                                self.transcripts[ensembl_id].refseq_id = refseq_id
                            if base_ensembl in self.transcripts:
                                self.transcripts[base_ensembl].is_mane_select = True
                                self.transcripts[base_ensembl].refseq_id = refseq_id

                            mane_count += 1

        if verbose:
            logger.info(f"Found {mapping_count} Ensembl-RefSeq mappings")
            logger.info(f"Found {mane_count} MANE Select transcripts")

    def select_best_transcript_for_position(
        self, gene_id: str, target_pos: int, strand: str, max_distance: int = 100
    ) -> Optional[str]:
        """Select the best transcript for a given genomic position with intelligent ranking."""
        if gene_id not in self.gene_transcripts:
            return None

        candidate_transcripts = []

        for transcript_id in self.gene_transcripts[gene_id]:
            transcript_info = self.transcripts[transcript_id]
            distance = abs(transcript_info.start_pos - target_pos)

            # Only consider transcripts within reasonable distance
            if distance <= max_distance:
                # Calculate priority score (lower is better)
                priority_score = 0

                # Exact match gets highest priority
                if distance == 0:
                    priority_score = 0
                else:
                    priority_score = distance

                # Boost preferred transcripts
                if transcript_info.is_preferred:
                    priority_score -= 1000

                # Boost MANE Select transcripts
                if transcript_info.is_mane_select:
                    priority_score -= 2000

                candidate_transcripts.append((transcript_id, priority_score, distance))

        if not candidate_transcripts:
            # Fallback: find closest transcript regardless of distance
            all_distances = []
            for transcript_id in self.gene_transcripts[gene_id]:
                transcript_info = self.transcripts[transcript_id]
                distance = abs(transcript_info.start_pos - target_pos)

                priority_score = distance
                if transcript_info.is_preferred:
                    priority_score -= 1000
                if transcript_info.is_mane_select:
                    priority_score -= 2000

                all_distances.append((transcript_id, priority_score, distance))

            candidate_transcripts = all_distances

        if candidate_transcripts:
            # Sort by priority score, then by distance
            best_transcript = min(candidate_transcripts, key=lambda x: (x[1], x[2]))
            return best_transcript[0]

        return None

    def get_refseq_id(self, transcript_id: str) -> str:
        """Get RefSeq ID for a given Ensembl transcript with comprehensive fallback logic."""
        # Method 1: Check if it's already stored in transcript info
        if (
            transcript_id in self.transcripts
            and self.transcripts[transcript_id].refseq_id
        ):
            return self.transcripts[transcript_id].refseq_id

        # Method 2: Direct mapping (versioned)
        if transcript_id in self.ensembl_to_refseq:
            return self.ensembl_to_refseq[transcript_id]

        # Method 3: Direct mapping (base, without version)
        base_id = transcript_id.split(".")[0]
        if base_id in self.ensembl_to_refseq:
            return self.ensembl_to_refseq[base_id]

        # Method 4: Fallback - check reverse mapping for unique matches
        # Look for RefSeq transcripts that have only one Ensembl match
        for refseq_id, ensembl_id in self.ensembl_to_refseq.items():
            # Check if this RefSeq maps uniquely to our transcript
            if transcript_id == ensembl_id or base_id == ensembl_id.split(".")[0]:
                # Verify it's a unique mapping by checking reverse
                matching_ensembls = [
                    ens
                    for ens, ref in self.ensembl_to_refseq.items()
                    if ref == refseq_id
                ]
                if len(matching_ensembls) == 1:
                    return refseq_id

        # Method 5: Gene-level fallback
        # If transcript has a gene, check if gene has only one RefSeq transcript
        if transcript_id in self.transcripts:
            gene_id = self.transcripts[transcript_id].gene_id
            gene_refseq_ids = set()

            # Collect all RefSeq IDs for this gene
            for other_transcript_id in self.gene_transcripts[gene_id]:
                if other_transcript_id in self.ensembl_to_refseq:
                    gene_refseq_ids.add(self.ensembl_to_refseq[other_transcript_id])

                other_base_id = other_transcript_id.split(".")[0]
                if other_base_id in self.ensembl_to_refseq:
                    gene_refseq_ids.add(self.ensembl_to_refseq[other_base_id])

            # If gene has exactly one RefSeq ID, use it
            if len(gene_refseq_ids) == 1:
                return list(gene_refseq_ids)[0]
            elif len(gene_refseq_ids) > 1:
                # Pick the best one (prefer NM_ over XM_, then shortest)
                best_refseq = min(
                    gene_refseq_ids, key=lambda x: (not x.startswith("NM_"), len(x), x)
                )
                return best_refseq

        return "NA"


@dataclass
class BedEntry:
    """Data class for BED file entries."""

    chrom: str
    start: int
    end: int
    name: str
    score: str
    strand: str
    gene_id: str
    gene_name: str
    start_type: str  # Annotated, Truncated, Extended
    start_codon: str
    efficiency: float
    start_pos: int
    line_num: int = -1


class BedProcessor:
    """Enhanced BED file processor with improved transcript mapping."""

    def __init__(self, transcript_db: TranscriptDatabase):
        """Initialize the processor with a transcript database and statistics tracker.

        Args:
            transcript_db: TranscriptDatabase instance used for mapping BED entries.
        """
        self.transcript_db = transcript_db
        self.stats = defaultdict(int)

    def parse_bed_entry(self, line: str, line_num: int) -> Optional[BedEntry]:
        """Parse a BED line into a BedEntry object."""
        fields = line.strip().split("\t")
        if len(fields) < 6:
            self.stats["invalid_format"] += 1
            return None

        try:
            name_info = self._parse_alternative_start_name(fields[3])
            gene_id = name_info["gene_id"].split(".")[0]

            # Validate gene exists in our database
            if gene_id not in self.transcript_db.gene_names:
                self.stats["invalid_gene"] += 1
                return None

            # Filter out uORFs
            if name_info["start_type"] == "uORF":
                self.stats["uorfs_filtered"] += 1
                return None

            # Only keep target types
            if name_info["start_type"] not in ["Annotated", "Truncated", "Extended"]:
                self.stats["filtered_type"] += 1
                return None

            chrom, start, end, strand = (
                fields[0],
                int(fields[1]),
                int(fields[2]),
                fields[5],
            )
            start_pos = start if strand == "+" else end
            gene_name = self.transcript_db.gene_names[gene_id]

            # Update name with correct gene name if needed
            if name_info["gene_name"].upper() != gene_name.upper():
                updated_name = f"{gene_name}_{name_info['gene_id']}_{name_info['start_type']}_{name_info['start_codon']}_{name_info['efficiency']}"
                self.stats["updated_gene_names"] += 1
            else:
                updated_name = fields[3]

            return BedEntry(
                chrom=chrom,
                start=start,
                end=end,
                name=updated_name,
                score=fields[4],
                strand=strand,
                gene_id=gene_id,
                gene_name=gene_name,
                start_type=name_info["start_type"],
                start_codon=name_info["start_codon"],
                efficiency=name_info["efficiency"],
                start_pos=start_pos,
                line_num=line_num,
            )

        except (IndexError, ValueError, KeyError) as e:
            self.stats["parse_error"] += 1
            return None

    def select_relevant_annotated_start(
        self, gene_id: str, alternative_pos: int, alternative_type: str, strand: str
    ) -> Optional[BedEntry]:
        """Select biologically relevant annotated start for pairing with alternative."""
        gene_entries = self.gene_entries[
            gene_id
        ]  # You'll need to store this in BedProcessor
        annotated_entries = [e for e in gene_entries if e.start_type == "Annotated"]

        if len(annotated_entries) <= 1:
            return annotated_entries[0] if annotated_entries else None

        relevant_starts = []

        for ann_entry in annotated_entries:
            ann_pos = ann_entry.start_pos

            if alternative_type == "Extended":
                if strand == "+":
                    # Extension upstream, canonical downstream
                    if ann_pos > alternative_pos:
                        distance = ann_pos - alternative_pos
                        relevant_starts.append((ann_entry, distance))
                else:
                    # Extension upstream in transcript (lower genomic coord)
                    if ann_pos < alternative_pos:
                        distance = alternative_pos - ann_pos
                        relevant_starts.append((ann_entry, distance))

            elif alternative_type == "Truncated":
                if strand == "+":
                    # Truncation downstream, canonical upstream
                    if ann_pos < alternative_pos:
                        distance = alternative_pos - ann_pos
                        relevant_starts.append((ann_entry, distance))
                else:
                    # Truncation downstream in transcript (higher genomic coord)
                    if ann_pos > alternative_pos:
                        distance = ann_pos - alternative_pos
                        relevant_starts.append((ann_entry, distance))

        if relevant_starts:
            # Return the closest biologically relevant annotated start
            return min(relevant_starts, key=lambda x: x[1])[0]
        else:
            # Fallback to highest efficiency if no biologically relevant pairing
            return max(annotated_entries, key=lambda x: x.efficiency)

    def _parse_alternative_start_name(self, name: str) -> Dict[str, str]:
        """Parse alternative start site name field."""
        parts = name.split("_")

        if len(parts) >= 5:
            return {
                "gene_name": parts[0],
                "gene_id": parts[1],
                "start_type": parts[2],
                "start_codon": parts[3],
                "efficiency": float(parts[4])
                if parts[4].replace(".", "").isdigit()
                else 0.0,
            }
        else:
            return {
                "gene_id": name,
                "gene_name": name,
                "start_type": "Unknown",
                "start_codon": "UNK",
                "efficiency": 0.0,
            }

    def add_missing_annotated_starts(
        self, gene_entries: Dict[str, List[BedEntry]]
    ) -> int:
        """Add missing annotated starts for genes that have alternatives but no annotated."""
        added_count = 0

        for gene_id, entries in gene_entries.items():
            has_annotated = any(e.start_type == "Annotated" for e in entries)
            has_alternatives = any(
                e.start_type in ["Truncated", "Extended"] for e in entries
            )

            if has_alternatives and not has_annotated:
                # Find the best transcript for this gene
                if gene_id in self.transcript_db.gene_transcripts:
                    # Use the first preferred/MANE Select transcript if available
                    best_transcript_id = None
                    for transcript_id in self.transcript_db.gene_transcripts[gene_id]:
                        transcript_info = self.transcript_db.transcripts[transcript_id]
                        if (
                            transcript_info.is_mane_select
                            or transcript_info.is_preferred
                        ):
                            best_transcript_id = transcript_id
                            break

                    # Fallback to first transcript
                    if not best_transcript_id:
                        best_transcript_id = self.transcript_db.gene_transcripts[
                            gene_id
                        ][0]

                    transcript_info = self.transcript_db.transcripts[best_transcript_id]

                    # Create BED coordinates for the start codon
                    if transcript_info.strand == "+":
                        bed_start = transcript_info.start
                        bed_end = transcript_info.start + 2
                    else:
                        bed_start = transcript_info.end - 2
                        bed_end = transcript_info.end

                    # Create annotated entry
                    annotated_entry = BedEntry(
                        chrom=transcript_info.chrom,
                        start=bed_start,
                        end=bed_end,
                        name=f"{transcript_info.gene_name}_{gene_id}_Annotated_AUG_0.0",
                        score="0",
                        strand=transcript_info.strand,
                        gene_id=gene_id,
                        gene_name=transcript_info.gene_name,
                        start_type="Annotated",
                        start_codon="AUG",
                        efficiency=0.0,
                        start_pos=transcript_info.start_pos,
                        line_num=-1,  # Mark as added
                    )

                    gene_entries[gene_id].append(annotated_entry)
                    added_count += 1

        return added_count

    def process_bed_file(
        self,
        input_bed_path: Union[str, Path],
        output_bed_path: Union[str, Path],
        verbose: bool = True,
    ) -> Dict:
        """Process BED file with enhanced transcript mapping and biologically relevant pairing."""
        if verbose:
            logger.info(f"Processing BED file: {input_bed_path}")

        # Parse all entries
        valid_entries = []
        gene_entries = defaultdict(list)

        with open(input_bed_path, "r") as f:
            for line_num, line in enumerate(f, 1):
                if not line.strip():
                    continue

                self.stats["total_entries"] += 1
                entry = self.parse_bed_entry(line, line_num)

                if entry:
                    valid_entries.append(entry)
                    gene_entries[entry.gene_id].append(entry)
                    self.stats["valid_entries"] += 1

        if verbose:
            logger.info(
                f"  ├─ Valid entries: {self.stats['valid_entries']} / {self.stats['total_entries']}"
            )

        # Add missing annotated starts
        added_annotated = self.add_missing_annotated_starts(gene_entries)
        if verbose:
            logger.info(f"  ├─ Added annotated starts: {added_annotated}")

        # Store gene entries for pairing logic
        self.gene_entries = gene_entries

        # Map entries to transcripts with biologically relevant pairing
        enhanced_entries = []
        transcript_assignments = 0

        for gene_id, entries in gene_entries.items():
            # Only process complete genes (have both annotated and alternatives)
            has_annotated = any(e.start_type == "Annotated" for e in entries)
            has_alternatives = any(
                e.start_type in ["Truncated", "Extended"] for e in entries
            )

            if not (has_annotated and has_alternatives):
                continue

            for entry in entries:
                if entry.start_type == "Annotated":
                    # For annotated starts, use direct position matching
                    transcript_id = (
                        self.transcript_db.select_best_transcript_for_position(
                            gene_id=entry.gene_id,
                            target_pos=entry.start_pos,
                            strand=entry.strand,
                            max_distance=100,
                        )
                    )
                else:
                    # For alternatives, find relevant annotated start first
                    relevant_annotated = self.select_relevant_annotated_start(
                        gene_id=gene_id,
                        alternative_pos=entry.start_pos,
                        alternative_type=entry.start_type,
                        strand=entry.strand,
                    )

                    if relevant_annotated:
                        # Use the transcript from the relevant annotated start
                        transcript_id = (
                            self.transcript_db.select_best_transcript_for_position(
                                gene_id=relevant_annotated.gene_id,
                                target_pos=relevant_annotated.start_pos,
                                strand=relevant_annotated.strand,
                                max_distance=100,
                            )
                        )
                    else:
                        # Fallback to direct position matching
                        transcript_id = (
                            self.transcript_db.select_best_transcript_for_position(
                                gene_id=entry.gene_id,
                                target_pos=entry.start_pos,
                                strand=entry.strand,
                                max_distance=100,
                            )
                        )

                if transcript_id is None:
                    transcript_id = "NA"
                else:
                    transcript_assignments += 1

                # Get RefSeq ID
                refseq_id = (
                    self.transcript_db.get_refseq_id(transcript_id)
                    if transcript_id != "NA"
                    else "NA"
                )

                # Create enhanced BED line
                enhanced_line = f"{entry.chrom}\t{entry.start}\t{entry.end}\t{entry.name}\t{entry.score}\t{entry.strand}\t{transcript_id}\t{refseq_id}"
                enhanced_entries.append(enhanced_line)

        # Write output
        output_bed_path = Path(output_bed_path)
        output_bed_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_bed_path, "w") as f:
            f.write("\n".join(enhanced_entries))

        # Calculate final statistics
        final_stats = dict(self.stats)
        final_stats.update(
            {
                "added_annotated": added_annotated,
                "enhanced_entries": len(enhanced_entries),
                "transcript_assignments": transcript_assignments,
                "final_genes": len(gene_entries),
            }
        )

        if verbose:
            logger.info(f"  ├─ Enhanced entries: {len(enhanced_entries)}")
            logger.info(f"  ├─ Transcript assignments: {transcript_assignments}")
            logger.info(f"  └─ Output saved: {output_bed_path}")

        return final_stats


def comprehensive_cleanup_bed_with_intelligent_selection(
    input_bed_path: Union[str, Path],
    gtf_path: Union[str, Path],
    preferred_transcripts_path: Union[str, Path],
    output_bed_path: Union[str, Path],
    refseq_gtf_path: Optional[Union[str, Path]] = None,
    verbose: bool = True,
) -> Dict:
    """Comprehensive BED cleanup with intelligent transcript selection.

    This function provides:
    1. Improved transcript selection using preferred lists and MANE Select
    2. Clean, systematic processing pipeline
    3. Enhanced RefSeq mapping with fallback logic
    4. Better error handling and statistics

    Args:
        input_bed_path: Path to input BED file
        gtf_path: Path to GENCODE GTF annotation file
        preferred_transcripts_path: Path to preferred transcripts file
        output_bed_path: Path to save enhanced BED file
        refseq_gtf_path: Path to RefSeq GTF file (auto-detected if None)
        verbose: Print detailed progress

    Returns:
        Dictionary with comprehensive processing statistics
    """
    if verbose:
        logger.info(
            "🔧 Starting comprehensive BED cleanup with intelligent transcript selection..."
        )
        logger.info(f"  ├─ Input BED: {input_bed_path}")
        logger.info(f"  ├─ GENCODE GTF: {gtf_path}")
        logger.info(f"  ├─ Preferred transcripts: {preferred_transcripts_path}")
        logger.info(f"  └─ Output BED: {output_bed_path}")

    # Auto-detect RefSeq GTF if not provided
    if refseq_gtf_path is None:
        genome_data_dir = Path(gtf_path).parent
        refseq_gtf_path = find_refseq_gtf_file(genome_data_dir)

        if refseq_gtf_path is None:
            if verbose:
                logger.warning("No RefSeq GTF found! RefSeq mapping will be limited.")
                logger.info(
                    "💡 Run: bash 0_download_genome.sh to get RefSeq annotations"
                )

    if verbose and refseq_gtf_path:
        logger.info(f"  ├─ RefSeq GTF: {refseq_gtf_path}")

    # Build comprehensive transcript database
    if verbose:
        logger.info("\n📊 Building transcript database...")

    transcript_db = TranscriptDatabase()
    transcript_db.load_from_gtf(gtf_path, verbose=verbose)
    transcript_db.load_preferred_transcripts(
        preferred_transcripts_path, verbose=verbose
    )

    if refseq_gtf_path:
        transcript_db.load_refseq_mapping(refseq_gtf_path, verbose=verbose)

    # Process BED file
    if verbose:
        logger.info("\n🔄 Processing BED file...")

    processor = BedProcessor(transcript_db)
    final_stats = processor.process_bed_file(
        input_bed_path, output_bed_path, verbose=verbose
    )

    if verbose:
        logger.info(f"\n✅ BED cleanup complete!")
        logger.info(f"  ├─ Enhanced entries: {final_stats['enhanced_entries']}")
        logger.info(
            f"  ├─ Transcript assignments: {final_stats['transcript_assignments']}"
        )
        logger.info(f"  ├─ Final genes: {final_stats['final_genes']}")
        logger.info(
            f"  └─ Success rate: {final_stats['transcript_assignments'] / final_stats['enhanced_entries'] * 100:.1f}%"
        )

    return final_stats


# Helper functions that you can add to your existing utils.py


def find_refseq_gtf_file(genome_data_dir: Union[str, Path]) -> Optional[Path]:
    """Find available RefSeq GTF file in the genome data directory."""
    genome_data_dir = Path(genome_data_dir)

    possible_files = [
        "GRCh38_latest_genomic.gtf",
        "ncbiRefSeq.gtf",
        "refseq_genomic.gtf",
        "GRCh38_refseq.gtf",
    ]

    for filename in possible_files:
        gtf_path = genome_data_dir / filename
        if gtf_path.exists():
            return gtf_path

    return None


def print_translation_summary(
    dataset: pd.DataFrame,
    successful_genes: int,
    total_genes: int,
    failed_genes: List[str],
    mutations_mode: bool,
    output_dir: str,
) -> None:
    """Print comprehensive summary of protein sequence generation results.

    Args:
        dataset: Generated protein sequence dataset
        successful_genes: Number of successfully processed genes
        total_genes: Total number of genes attempted
        failed_genes: List of gene names that failed processing
        mutations_mode: Whether mutations were included
        output_dir: Output directory path
    """
    logger.info(f"\nProtein Sequence Generation Summary:")
    logger.info(f"  ├─ Total genes processed: {total_genes}")

    # Status breakdown
    logger.info(f"\n  ├─ Status breakdown:")
    logger.info(f"  │  ├─ Success: {successful_genes}")
    if failed_genes:
        logger.info(f"  │  └─ Failed: {len(failed_genes)}")
    else:
        logger.info(f"  │  └─ Failed: 0")

    # Dataset statistics
    if not dataset.empty:
        logger.info(f"\n  ├─ Sequence Generation:")
        logger.info(f"  │  ├─ Total sequences generated: {len(dataset):,}")

        # Calculate transcript-isoform pairs
        genes_with_data = dataset["gene"].nunique()
        if mutations_mode and "variant_type" in dataset.columns:
            # Count unique transcript-isoform pairs (canonical + alternative base pairs)
            base_sequences = dataset[
                dataset["variant_type"].isin(
                    ["canonical", "alternative"]
                )  # Updated terminology
            ]
            unique_pairs = (
                len(base_sequences) // 2
                if len(base_sequences) % 2 == 0
                else (len(base_sequences) + 1) // 2
            )
        else:
            # For pairs mode, total sequences / 2 = pairs
            unique_pairs = (
                len(dataset) // 2 if len(dataset) % 2 == 0 else (len(dataset) + 1) // 2
            )

        logger.info(
            f"  │  ├─ Transcript-isoform pairs: {unique_pairs}"
        )  # Updated terminology
        logger.info(
            f"  │  └─ Average sequences per gene: {len(dataset) / genes_with_data:.1f}"
        )

        # Mode-specific breakdown
        if mutations_mode and "variant_type" in dataset.columns:
            logger.info(f"\n  ├─ Sequence breakdown:")
            type_counts = dataset["variant_type"].value_counts()
            for variant_type, count in type_counts.items():
                percentage = (count / len(dataset)) * 100
                logger.info(f"  │  ├─ {variant_type}: {count:,} ({percentage:.1f}%)")

            # Mutation-specific statistics
            if "canonical_mutated" in type_counts:
                mutation_data = dataset[dataset["variant_type"] == "canonical_mutated"]
                if not mutation_data.empty:
                    logger.info(f"\n  ├─ Mutation Analysis:")
                    logger.info(
                        f"  │  ├─ Total mutations integrated: {len(mutation_data)}"
                    )
                    logger.info(
                        f"  │  ├─ Unique mutation positions: {mutation_data['mutation_position'].nunique()}"
                    )

                    if "aa_change" in mutation_data.columns:
                        aa_changes = mutation_data["aa_change"].dropna()
                        silent_mutations = len(mutation_data) - len(aa_changes)
                        logger.info(
                            f"  │  ├─ Mutations with AA changes: {len(aa_changes)}"
                        )
                        if silent_mutations > 0:
                            logger.info(f"  │  ├─ Silent mutations: {silent_mutations}")

                    # Breakdown by impact type if available
                    if "mutation_impact" in mutation_data.columns:
                        impact_counts = mutation_data["mutation_impact"].value_counts()
                        logger.info(f"  │  │")
                        logger.info(f"  │  ├─ Breakdown by impact type:")
                        for impact_type, count in impact_counts.items():
                            percentage = (count / len(mutation_data)) * 100
                            logger.info(
                                f"  │  │  ├─ {impact_type}: {count} ({percentage:.1f}%)"
                            )

                    logger.info(
                        f"  │  └─ Average mutations per gene: {len(mutation_data) / genes_with_data:.1f}"
                    )
        else:
            # Pairs mode breakdown
            if "is_alternative" in dataset.columns:  # Updated column name
                canonical_count = len(dataset[dataset["is_alternative"] == 0])
                alternative_count = len(dataset[dataset["is_alternative"] == 1])
                logger.info(f"\n  ├─ Sequence breakdown:")
                logger.info(f"  │  ├─ Canonical: {canonical_count:,}")
                logger.info(
                    f"  │  └─ Alternative isoform: {alternative_count:,}"
                )  # Updated terminology

        # Length statistics
        logger.info(f"\n  ├─ Length statistics:")
        logger.info(f"  │  ├─ Average: {dataset['length'].mean():.1f} amino acids")
        logger.info(
            f"  │  ├─ Range: {dataset['length'].min()}-{dataset['length'].max()}"
        )
        logger.info(f"  │  └─ Median: {dataset['length'].median():.1f}")

    else:
        logger.info(f"\n  ├─ No sequences generated")

    # Failed genes details
    if failed_genes:
        logger.info(f"\n  ├─ Genes with errors:")
        for gene in failed_genes[:5]:  # Show first 5
            logger.info(
                f"  │  ├─ {gene}: No transcript-isoform pairs found"
            )  # Updated terminology
        if len(failed_genes) > 5:
            logger.info(f"  │  └─ ... and {len(failed_genes) - 5} more")

    # Output files
    output_path = Path(output_dir)
    logger.info(f"\n  ├─ Results saved to: {output_dir}")

    if mutations_mode:
        fasta_file = output_path / "protein_sequences_with_mutations.fasta"
        csv_file = output_path / "protein_sequences_with_mutations.csv"
        logger.info(f"  ├─ Protein sequences with mutations saved to: {csv_file.name}")
        if fasta_file.exists():
            logger.info(f"  ├─ FASTA format saved to: {fasta_file.name}")
    else:
        fasta_file = output_path / "protein_sequences_pairs.fasta"
        csv_file = output_path / "protein_sequences_pairs.csv"
        logger.info(f"  ├─ Protein sequence pairs saved to: {csv_file.name}")
        if fasta_file.exists():
            logger.info(f"  ├─ FASTA format saved to: {fasta_file.name}")


def update_gencode_gene_names(
    input_gtf_path: Union[str, Path],
    output_gtf_path: Union[str, Path],
    reference_gtf_path: Union[str, Path],
    verbose: bool = True,
) -> Dict:
    """Update gene names in GENCODE GTF using an updated GENCODE reference GTF.

    Args:
        input_gtf_path: Path to input GENCODE GTF file to be updated
        output_gtf_path: Path to save updated GTF file
        reference_gtf_path: Path to reference GTF file with preferred gene names
        verbose: Print detailed update summary

    Returns:
        Dictionary containing update statistics
    """
    input_gtf_path, output_gtf_path = Path(input_gtf_path), Path(output_gtf_path)

    # Extract Ensembl IDs to gene name mappings from both files
    if verbose:
        logger.info(
            f"Creating gene ID to name mappings from reference GTF: {reference_gtf_path}"
        )

    # Dictionary to store updated gene names: key=Ensembl ID (without version), value=updated gene name
    gene_name_updates = {}

    # First, get the current gene names from GENCODE
    gencode_gene_names = {}
    with open(input_gtf_path, "r") as f:
        for line in f:
            if line.startswith("#") or len(line.strip()) == 0:
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            # Extract Ensembl ID and gene name
            attributes = fields[8]
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)

            if gene_id_match and gene_name_match:
                full_gene_id = gene_id_match.group(1)
                gene_name = gene_name_match.group(1)

                # Remove version from Ensembl ID (e.g., ENSG00000264520.1 -> ENSG00000264520)
                base_gene_id = full_gene_id.split(".")[0]
                gencode_gene_names[base_gene_id] = gene_name

    if verbose:
        logger.info(f"Extracted {len(gencode_gene_names)} gene names from GENCODE GTF")

    # Next, get the gene names from the reference GTF
    reference_gene_names = {}
    with open(reference_gtf_path, "r") as f:
        for line in f:
            if line.startswith("#") or len(line.strip()) == 0:
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            # Extract Ensembl ID and gene name
            attributes = fields[8]

            # Look for Ensembl ID pattern (with or without version)
            ensembl_id_match = re.search(r'gene_id "((ENSG\d+)(?:\.\d+)?)"', attributes)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)

            if ensembl_id_match and gene_name_match:
                base_gene_id = ensembl_id_match.group(
                    2
                )  # This gets the ID without version
                gene_name = gene_name_match.group(1)

                reference_gene_names[base_gene_id] = gene_name

    if verbose:
        logger.info(
            f"Extracted {len(reference_gene_names)} gene names from reference GTF"
        )

    # Create update mapping
    for ensembl_id, gencode_name in gencode_gene_names.items():
        if ensembl_id in reference_gene_names:
            reference_name = reference_gene_names[ensembl_id]

            # Only update if names differ
            if gencode_name != reference_name:
                gene_name_updates[ensembl_id] = reference_name

    if verbose:
        logger.info(f"Created {len(gene_name_updates)} gene name updates")

    # Process GENCODE GTF file to update gene names
    stats = {
        "total_lines": 0,
        "updated_lines": 0,
        "genes_processed": 0,
        "genes_updated": 0,
    }

    # Create output directory if it doesn't exist
    output_dir = output_gtf_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(input_gtf_path, "r") as f_in, open(output_gtf_path, "w") as f_out:
        for line in f_in:
            stats["total_lines"] += 1

            # Pass through comment lines
            if line.startswith("#"):
                f_out.write(line)
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                f_out.write(line)
                continue

            # Parse attributes
            attributes = fields[8]
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)

            if gene_id_match:
                full_gene_id = gene_id_match.group(1)
                base_gene_id = full_gene_id.split(".")[0]

                if fields[2] == "gene":
                    stats["genes_processed"] += 1

                # Check if we have an updated name for this gene
                if base_gene_id in gene_name_updates:
                    new_name = gene_name_updates[base_gene_id]

                    # Replace the gene_name attribute
                    new_attributes = re.sub(
                        r'gene_name "([^"]+)"', f'gene_name "{new_name}"', attributes
                    )

                    if new_attributes != attributes:
                        fields[8] = new_attributes
                        stats["updated_lines"] += 1

                        if fields[2] == "gene":
                            stats["genes_updated"] += 1

                        line = "\t".join(fields) + "\n"

            f_out.write(line)

    if verbose:
        logger.info(f"\nGTF Update Summary:")
        logger.info(f"  Total lines processed: {stats['total_lines']}")
        logger.info(f"  Genes processed: {stats['genes_processed']}")
        logger.info(f"  Genes with updated names: {stats['genes_updated']}")
        logger.info(f"  Total lines updated: {stats['updated_lines']}")
        logger.info(f"  Output saved to: {output_gtf_path}")

    return stats


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
    logger.info(f"Read {len(gene_names)} genes from {gene_list_path}")
    return gene_names


def subset_gene_list(
    gene_list_path: Union[str, Path],
    subset_genes: List[str] = [
        # Truncations of interest
        "AKR7A2",
        "ALDH9A1",
        "C15orf40",
        "CHCHD1",
        "CMPK1",
        "ERGIC3",
        "FAAH2",
        "FH",
        "GADD45GIP1",
        "GARS1",
        "GLRX2",
        "GSR",
        "LAGE3",
        "MYG1",
        "NAXE",
        "NTHL1",
        "PCBD2",
        "PNPO",
        "REXO2",
        "TRNT1",
        "UXS1",
        # Extensions of interest
        "PTEN",
        "FTL",
        "MAPK14",
    ],
) -> List[str]:
    """Subset a gene list to only include specified genes.

    Args:
        gene_list_path: Path to the file containing the full gene list
        subset_genes: List of gene names to include in the subset (defaults to a predefined list)

    Returns:
        List of genes from the full list that are in the subset
    """
    full_gene_list = parse_gene_list(gene_list_path)
    subset_set = set(subset_genes)
    subsetted_genes = [gene for gene in full_gene_list if gene in subset_set]
    logger.info(
        f"Subsetted to {len(subsetted_genes)} genes from {len(full_gene_list)} total"
    )
    return subsetted_genes


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
    logger.info(f"Gene-level results saved to {output_path}")


def save_isoform_level_results(
    results: List[Dict],
    output_dir: Union[str, Path],
    filename: str = "isoform_level_results.csv",
) -> None:
    """Save detailed transcript-isoform pair analysis results to CSV file.

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
        logger.info(f"No transcript-isoform pairs found to save")
        return

    # Convert to DataFrame
    pairs_df = pd.DataFrame(all_pairs)

    # Organize columns in preferred order
    column_order = [
        "gene_name",  # Gene name first
        "transcript_id",  # Transcript ID second
        "feature_id",  # Alternative feature ID third
        "feature_type",  # truncation or extension
        "feature_start",
        "feature_end",
        "mutation_count_total",  # Total mutation count
    ]

    # Add mutation category columns to the order (they start with "mutations_")
    mutation_category_columns = [
        col for col in pairs_df.columns if col.startswith("mutations_")
    ]
    column_order.extend(mutation_category_columns)

    # Add variant IDs columns
    variant_id_columns = [
        col for col in pairs_df.columns if col.startswith("variant_ids_")
    ]
    column_order.extend(variant_id_columns)

    # Reorder columns that exist in the dataframe
    available_columns = [col for col in column_order if col in pairs_df.columns]

    # Add any remaining columns that weren't explicitly ordered
    remaining_columns = [
        col for col in pairs_df.columns if col not in available_columns
    ]
    final_column_order = available_columns + remaining_columns

    # Apply the column order and save
    if final_column_order:  # Only reorder if we have columns to reorder
        pairs_df = pairs_df[final_column_order]

    # Save to CSV
    output_path = output_dir / filename
    pairs_df.to_csv(output_path, index=False)
    logger.info(f"Transcript-isoform pair analysis saved to {output_path}")


def print_mutation_summary(results_df, output_dir):
    """Print summary statistics from the analysis results.

    Args:
        results_df: DataFrame containing analysis results
        output_dir: Directory where output files are saved
    """
    logger.info("\nAnalysis Summary:")
    logger.info(f"  ├─ Total genes processed: {len(results_df)}")
    logger.info("\n  ├─ Status breakdown:")
    for status, count in results_df["status"].value_counts().items():
        logger.info(f"  │  ├─ {status}: {count}")

    # Transcript-isoform statistics (updated for new structure)
    successful_genes = results_df[results_df["status"] == "success"]
    if not successful_genes.empty:
        total_transcripts = successful_genes["total_transcripts"].sum()
        total_features = successful_genes["alternative_features"].sum()
        total_pairs = successful_genes["transcript_feature_pairs"].sum()

        # Calculate average pairs per gene
        avg_pairs_per_gene = (
            total_pairs / len(successful_genes) if len(successful_genes) > 0 else 0
        )

        logger.info(f"\n  ├─ Transcript-Isoform Analysis:")
        logger.info(f"  │  ├─ Total transcripts across all genes: {total_transcripts}")
        logger.info(f"  │  ├─ Total alternative isoform features: {total_features}")
        logger.info(f"  │  ├─ Total transcript-isoform pairs: {total_pairs}")
        logger.info(f"  │  └─ Average pairs per gene: {avg_pairs_per_gene:.2f}")

        # Mutation statistics if available
        if "mutations_filtered" in successful_genes.columns:
            total_mutations = successful_genes["mutations_filtered"].sum()
            logger.info(f"\n  ├─ Mutation Analysis:")
            logger.info(
                f"  │  ├─ Total mutations in alternative isoform regions: {total_mutations}"
            )

            # Try to load mutation analysis results to report statistics
            mutation_analysis_path = (
                Path(output_dir) / "isoform_level_results.csv"
            )  # Updated filename
            if mutation_analysis_path.exists():
                try:
                    pairs_df = pd.read_csv(mutation_analysis_path)

                    # Get total mutations across all pairs
                    total_pair_mutations = (
                        pairs_df["mutation_count_total"].sum()
                        if "mutation_count_total" in pairs_df.columns
                        else 0
                    )
                    avg_mutations_per_pair = (
                        total_pair_mutations / len(pairs_df) if len(pairs_df) > 0 else 0
                    )
                    logger.info(
                        f"  │  ├─ Total mutations across all pairs: {total_pair_mutations}"
                    )
                    logger.info(
                        f"  │  ├─ Average mutations per transcript-isoform pair: {avg_mutations_per_pair:.2f}"
                    )

                    # Show breakdown by isoform type if available
                    if "feature_type" in pairs_df.columns:
                        logger.info(f"  │  │")
                        logger.info(f"  │  ├─ Breakdown by isoform type:")

                        type_mutations = pairs_df.groupby("feature_type")[
                            "mutation_count_total"
                        ].sum()
                        for isoform_type, count in type_mutations.items():
                            type_pairs = len(
                                pairs_df[pairs_df["feature_type"] == isoform_type]
                            )
                            avg_per_type = count / type_pairs if type_pairs > 0 else 0
                            logger.info(
                                f"  │  │  ├─ {isoform_type.capitalize()}: {count} total ({avg_per_type:.1f} avg/pair)"
                            )

                    # Print statistics for each mutation category
                    mutation_categories = [
                        col for col in pairs_df.columns if col.startswith("mutations_")
                    ]
                    if mutation_categories:
                        logger.info(f"  │  │")
                        logger.info(f"  │  ├─ Breakdown by mutation category:")

                        for category in mutation_categories:
                            # Skip if the column doesn't exist in the dataframe
                            if category not in pairs_df.columns:
                                continue

                            # Convert category name from mutations_missense_variant to "Missense Variant"
                            category_name = category.replace("mutations_", "").replace(
                                "_", " "
                            )
                            category_name = category_name.title()

                            # Calculate statistics for this category
                            category_total = pairs_df[category].sum()
                            category_percent = (
                                (category_total / total_pair_mutations * 100)
                                if total_pair_mutations > 0
                                else 0
                            )

                            logger.info(
                                f"  │  │  ├─ {category_name}: {category_total} ({category_percent:.1f}%)"
                            )

                    logger.info(
                        f"  │  └─ Detailed results available in isoform_level_results.csv"  # Updated filename
                    )
                except Exception as e:
                    logger.error(f"Error reading mutation analysis: {str(e)}")
                    logger.info(
                        f"  │  └─ Error reading detailed mutation analysis: {str(e)}"
                    )

    # Genes with errors
    error_genes = results_df[results_df["status"] == "error"]
    if not error_genes.empty:
        logger.info("\n  ├─ Genes with errors:")
        for _, row in error_genes.iterrows():
            logger.info(f"  │  ├─ {row['gene_name']}: {row['error']}")

    logger.info(f"\n  ├─ Results saved to: {output_dir}")
    logger.info(
        f"  ├─ Gene-level results saved to: {Path(output_dir) / 'gene_level_results.csv'}"
    )
    logger.info(
        f"  ├─ Detailed mutation analysis by pair saved to: {Path(output_dir) / 'isoform_level_results.csv'}"  # Updated filename
    )


def print_translation_summary(
    dataset: pd.DataFrame,
    successful_genes: int,
    total_genes: int,
    failed_genes: List[str],
    mutations_mode: bool,
    output_dir: str,
) -> None:
    """Print comprehensive summary of protein sequence generation results.

    Args:
        dataset: Generated protein sequence dataset
        successful_genes: Number of successfully processed genes
        total_genes: Total number of genes attempted
        failed_genes: List of gene names that failed processing
        mutations_mode: Whether mutations were included
        output_dir: Output directory path
    """
    logger.info(f"\nProtein Sequence Generation Summary:")
    logger.info(f"  ├─ Total genes processed: {total_genes}")

    # Status breakdown
    logger.info(f"\n  ├─ Status breakdown:")
    logger.info(f"  │  ├─ Success: {successful_genes}")
    if failed_genes:
        logger.info(f"  │  └─ Failed: {len(failed_genes)}")
    else:
        logger.info(f"  │  └─ Failed: 0")

    # Dataset statistics
    if not dataset.empty:
        logger.info(f"\n  ├─ Sequence Generation:")
        logger.info(f"  │  ├─ Total sequences generated: {len(dataset):,}")

        # Calculate transcript-isoform pairs
        genes_with_data = dataset["gene"].nunique()
        if mutations_mode and "variant_type" in dataset.columns:
            # Count unique transcript-isoform pairs (canonical + alternative base pairs)
            base_sequences = dataset[
                dataset["variant_type"].isin(
                    ["canonical", "alternative"]
                )  # Updated terminology
            ]
            unique_pairs = (
                len(base_sequences) // 2
                if len(base_sequences) % 2 == 0
                else (len(base_sequences) + 1) // 2
            )
        else:
            # For pairs mode, total sequences / 2 = pairs
            unique_pairs = (
                len(dataset) // 2 if len(dataset) % 2 == 0 else (len(dataset) + 1) // 2
            )

        logger.info(
            f"  │  ├─ Transcript-isoform pairs: {unique_pairs}"
        )  # Updated terminology
        logger.info(
            f"  │  └─ Average sequences per gene: {len(dataset) / genes_with_data:.1f}"
        )

        # Mode-specific breakdown
        if mutations_mode and "variant_type" in dataset.columns:
            logger.info(f"\n  ├─ Sequence breakdown:")
            type_counts = dataset["variant_type"].value_counts()
            for variant_type, count in type_counts.items():
                percentage = (count / len(dataset)) * 100
                logger.info(f"  │  ├─ {variant_type}: {count:,} ({percentage:.1f}%)")

            # Mutation-specific statistics
            if "canonical_mutated" in type_counts:
                mutation_data = dataset[dataset["variant_type"] == "canonical_mutated"]
                if not mutation_data.empty:
                    logger.info(f"\n  ├─ Mutation Analysis:")
                    logger.info(
                        f"  │  ├─ Total mutations integrated: {len(mutation_data)}"
                    )
                    logger.info(
                        f"  │  ├─ Unique mutation positions: {mutation_data['mutation_position'].nunique()}"
                    )

                    if "aa_change" in mutation_data.columns:
                        aa_changes = mutation_data["aa_change"].dropna()
                        silent_mutations = len(mutation_data) - len(aa_changes)
                        logger.info(
                            f"  │  ├─ Mutations with AA changes: {len(aa_changes)}"
                        )
                        if silent_mutations > 0:
                            logger.info(f"  │  ├─ Silent mutations: {silent_mutations}")

                    # Breakdown by impact type if available
                    if "mutation_impact" in mutation_data.columns:
                        impact_counts = mutation_data["mutation_impact"].value_counts()
                        logger.info(f"  │  │")
                        logger.info(f"  │  ├─ Breakdown by impact type:")
                        for impact_type, count in impact_counts.items():
                            percentage = (count / len(mutation_data)) * 100
                            logger.info(
                                f"  │  │  ├─ {impact_type}: {count} ({percentage:.1f}%)"
                            )

                    logger.info(
                        f"  │  └─ Average mutations per gene: {len(mutation_data) / genes_with_data:.1f}"
                    )
        else:
            # Pairs mode breakdown
            if "is_alternative" in dataset.columns:  # Updated column name
                canonical_count = len(dataset[dataset["is_alternative"] == 0])
                alternative_count = len(dataset[dataset["is_alternative"] == 1])
                logger.info(f"\n  ├─ Sequence breakdown:")
                logger.info(f"  │  ├─ Canonical: {canonical_count:,}")
                logger.info(
                    f"  │  └─ Alternative isoform: {alternative_count:,}"
                )  # Updated terminology

        # Length statistics
        logger.info(f"\n  ├─ Length statistics:")
        logger.info(f"  │  ├─ Average: {dataset['length'].mean():.1f} amino acids")
        logger.info(
            f"  │  ├─ Range: {dataset['length'].min()}-{dataset['length'].max()}"
        )
        logger.info(f"  │  └─ Median: {dataset['length'].median():.1f}")

    else:
        logger.info(f"\n  ├─ No sequences generated")

    # Failed genes details
    if failed_genes:
        logger.info(f"\n  ├─ Genes with errors:")
        for gene in failed_genes[:5]:  # Show first 5
            logger.info(
                f"  │  ├─ {gene}: No transcript-isoform pairs found"
            )  # Updated terminology
        if len(failed_genes) > 5:
            logger.info(f"  │  └─ ... and {len(failed_genes) - 5} more")

    # Output files
    output_path = Path(output_dir)
    logger.info(f"\n  ├─ Results saved to: {output_dir}")

    if mutations_mode:
        fasta_file = output_path / "protein_sequences_with_mutations.fasta"
        csv_file = output_path / "protein_sequences_with_mutations.csv"
        logger.info(f"  ├─ Protein sequences with mutations saved to: {csv_file.name}")
        if fasta_file.exists():
            logger.info(f"  ├─ FASTA format saved to: {fasta_file.name}")
    else:
        fasta_file = output_path / "protein_sequences_pairs.fasta"
        csv_file = output_path / "protein_sequences_pairs.csv"
        logger.info(f"  ├─ Protein sequence pairs saved to: {csv_file.name}")
        if fasta_file.exists():
            logger.info(f"  ├─ FASTA format saved to: {fasta_file.name}")


def load_pre_validated_variants(mutations_file: str) -> Dict[str, Set[str]]:
    """Load pre-validated variant IDs from step 2 results.

    Args:
        mutations_file (str): Path to isoform_level_results.csv

    Returns:
        Dict[str, Set[str]]: Dictionary mapping gene_name -> set of variant IDs
    """
    logger.info(f"Loading pre-validated variant IDs from {mutations_file}")

    # Read the mutation results
    mutations_df = pd.read_csv(mutations_file)

    if mutations_df.empty:
        logger.warning("No mutation records found in results file")
        return {}

    logger.info(f"Found {len(mutations_df)} mutation records")

    # Define the impact-specific variant ID columns
    variant_id_columns = [
        "variant_ids_missense_variant",
        "variant_ids_nonsense_variant",
        "variant_ids_frameshift_variant",
        "variant_ids_synonymous_variant",
        "variant_ids_inframe_deletion",
        "variant_ids_inframe_insertion",
    ]

    # Extract variant IDs by gene
    pre_validated_variants = {}
    total_variant_ids = 0

    for _, row in mutations_df.iterrows():
        gene_name = row.get("gene_name", "")
        if not gene_name:
            continue

        if gene_name not in pre_validated_variants:
            pre_validated_variants[gene_name] = set()

        # Collect variant IDs from all impact type columns
        for col in variant_id_columns:
            if col in row and pd.notna(row[col]) and str(row[col]).strip():
                # Parse comma-separated variant IDs
                variant_ids = [
                    vid.strip() for vid in str(row[col]).split(",") if vid.strip()
                ]
                pre_validated_variants[gene_name].update(variant_ids)
                total_variant_ids += len(variant_ids)

    genes_with_variants = len([g for g in pre_validated_variants.values() if g])
    logger.info(
        f"Loaded {total_variant_ids} pre-validated variant IDs across {genes_with_variants} genes"
    )

    # Print summary by gene (top 10)
    gene_counts = [
        (gene, len(variants))
        for gene, variants in pre_validated_variants.items()
        if variants
    ]
    gene_counts.sort(key=lambda x: x[1], reverse=True)

    if gene_counts:
        logger.info("Top genes by variant count:")
        for gene, count in gene_counts[:10]:
            logger.info(f"  ├─ {gene}: {count} variants")
        if len(gene_counts) > 10:
            logger.info(f"  └─ ... and {len(gene_counts) - 10} more genes")

    return pre_validated_variants
