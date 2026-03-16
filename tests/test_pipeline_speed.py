"""Pipeline speed profiling tests.

Times each major operation and outputs a summary.
"""

import time
import pytest
from swissisoform.genome import GenomeHandler
from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.translation import AlternativeProteinGenerator
from tests.conftest import GENOME_PATH, GTF_PATH, BED_PATH, ALL_TEST_GENES


# Collect timing results
_timings = {}


def _record(name, elapsed):
    _timings[name] = elapsed


class TestPipelineSpeed:
    """Profile major pipeline operations."""

    def test_genome_loading_time(self):
        """Time genome FASTA loading (SeqIO.to_dict)."""
        start = time.time()
        from Bio import SeqIO
        genome = SeqIO.to_dict(SeqIO.parse(GENOME_PATH, "fasta"))
        elapsed = time.time() - start
        _record("genome_fasta_loading", elapsed)
        print(f"\nGenome FASTA loading: {elapsed:.2f}s ({len(genome)} sequences)")

    def test_gtf_parsing_time(self):
        """Time GTF annotation parsing."""
        start = time.time()
        handler = GenomeHandler.__new__(GenomeHandler)
        handler.genome = {}  # Dummy
        handler.gtf_path = GTF_PATH
        handler.load_annotations(GTF_PATH)
        elapsed = time.time() - start
        _record("gtf_parsing", elapsed)
        n_features = len(handler.annotations)
        print(f"\nGTF parsing: {elapsed:.2f}s ({n_features} features)")

    def test_bed_parsing_time(self):
        """Time BED file loading."""
        start = time.time()
        handler = AlternativeIsoform()
        handler.load_bed(BED_PATH)
        elapsed = time.time() - start
        _record("bed_parsing", elapsed)
        print(f"\nBED parsing: {elapsed:.2f}s ({len(handler.start_sites)} sites)")

    def test_canonical_extraction_per_gene(self, protein_generator, genome_handler):
        """Time canonical protein extraction per gene."""
        times = {}
        for gene_name in ALL_TEST_GENES:
            gene_features = genome_handler.find_gene_features(gene_name)
            start_codons = gene_features[gene_features["feature_type"] == "start_codon"]
            if start_codons.empty:
                continue

            tid = start_codons.iloc[0]["transcript_id"]
            start = time.time()
            protein_generator.extract_canonical_protein(tid)
            elapsed = time.time() - start
            times[gene_name] = elapsed

        avg = sum(times.values()) / len(times) if times else 0
        _record("canonical_extraction_avg", avg)
        print(f"\nCanonical extraction times:")
        for gene, t in sorted(times.items(), key=lambda x: -x[1]):
            print(f"  {gene}: {t*1000:.1f}ms")
        print(f"  Average: {avg*1000:.1f}ms")

    def test_extension_extraction_per_gene(
        self, protein_generator, alt_isoform_handler
    ):
        """Time extension protein extraction per gene."""
        times = {}
        for gene_name in ALL_TEST_GENES:
            features = alt_isoform_handler.get_translation_features(gene_name)
            if features.empty:
                continue
            ext_features = features[features["region_type"] == "extension"]
            if ext_features.empty:
                continue

            feature = ext_features.iloc[0]
            tid = feature["transcript_id"]
            start = time.time()
            protein_generator.extract_alternative_protein(tid, feature)
            elapsed = time.time() - start
            times[gene_name] = elapsed

        avg = sum(times.values()) / len(times) if times else 0
        _record("extension_extraction_avg", avg)
        print(f"\nExtension extraction times:")
        for gene, t in sorted(times.items(), key=lambda x: -x[1]):
            print(f"  {gene}: {t*1000:.1f}ms")
        print(f"  Average: {avg*1000:.1f}ms")

    def test_feature_lookup_time(self, genome_handler):
        """Time per-gene feature lookup (DataFrame filtering)."""
        times = {}
        for gene_name in ALL_TEST_GENES:
            start = time.time()
            for _ in range(100):
                genome_handler.find_gene_features(gene_name)
            elapsed = (time.time() - start) / 100
            times[gene_name] = elapsed

        avg = sum(times.values()) / len(times) if times else 0
        _record("feature_lookup_avg", avg)
        print(f"\nFeature lookup times (per call, averaged over 100):")
        for gene, t in sorted(times.items(), key=lambda x: -x[1]):
            print(f"  {gene}: {t*1000:.2f}ms")
        print(f"  Average: {avg*1000:.2f}ms")

    def test_print_timing_summary(self):
        """Print final timing summary."""
        if not _timings:
            pytest.skip("No timing data collected")

        print("\n" + "=" * 60)
        print("PIPELINE SPEED SUMMARY")
        print("=" * 60)

        # Sort by time descending
        for name, elapsed in sorted(_timings.items(), key=lambda x: -x[1]):
            if elapsed >= 1.0:
                print(f"  {name}: {elapsed:.2f}s")
            else:
                print(f"  {name}: {elapsed*1000:.1f}ms")

        print("=" * 60)

        # Identify top bottlenecks
        sorted_timings = sorted(_timings.items(), key=lambda x: -x[1])
        if sorted_timings:
            print("\nTOP BOTTLENECKS:")
            for i, (name, elapsed) in enumerate(sorted_timings[:3], 1):
                print(f"  {i}. {name}: {elapsed:.2f}s")
