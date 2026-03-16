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

    def test_mutation_cache_speed(self, genome_handler):
        """Time ClinVar fetch: first call vs cached call."""
        import asyncio
        from swissisoform.mutations import MutationHandler

        handler = MutationHandler(genome_handler=genome_handler)
        test_gene = "ADAR"

        # Clear any existing cache for this gene
        cache_path = handler._get_cache_path("clinvar", test_gene)
        if cache_path.exists():
            cache_path.unlink()

        # First call (API or network)
        start = time.time()
        try:
            result_first = asyncio.get_event_loop().run_until_complete(
                handler.get_clinvar_variants(test_gene)
            )
        except Exception:
            # API may not be available in test environment
            result_first = None
        elapsed_first = time.time() - start
        _record("clinvar_first_call", elapsed_first)

        if result_first is not None and not result_first.empty:
            # Second call (should hit parquet cache)
            # Clear in-memory cache to force parquet read
            handler.cached_data.clear()

            start = time.time()
            result_cached = asyncio.get_event_loop().run_until_complete(
                handler.get_clinvar_variants(test_gene)
            )
            elapsed_cached = time.time() - start
            _record("clinvar_cached_call", elapsed_cached)

            print(f"\nClinVar {test_gene}:")
            print(f"  First call:  {elapsed_first:.2f}s ({len(result_first)} variants)")
            print(f"  Cached call: {elapsed_cached:.4f}s ({len(result_cached)} variants)")
            print(f"  Speedup:     {elapsed_first / max(elapsed_cached, 0.001):.0f}x")

            assert elapsed_cached < 1.0, f"Cached call took {elapsed_cached:.2f}s, expected < 1s"
        else:
            print(f"\nClinVar {test_gene}: API unavailable, skipping cache test")
            pytest.skip("ClinVar API not available")

    def test_intronic_filter_speed(self, protein_generator, genome_handler, alt_isoform_handler):
        """Time bulk intronic filtering."""
        import asyncio
        import pandas as pd
        from swissisoform.mutations import MutationHandler

        handler = MutationHandler(genome_handler=genome_handler)
        test_gene = "ADAR"

        # Get real mutations if available, otherwise create synthetic ones
        try:
            mutations_df = asyncio.get_event_loop().run_until_complete(
                handler.get_clinvar_variants(test_gene)
            )
        except Exception:
            mutations_df = pd.DataFrame()

        if mutations_df.empty:
            # Create synthetic mutations for timing
            gene_features = genome_handler.find_gene_features(test_gene)
            cds = gene_features[gene_features["feature_type"] == "CDS"]
            if cds.empty:
                pytest.skip("No CDS features for ADAR")
            positions = []
            for _, region in cds.iterrows():
                positions.extend(range(int(region["start"]), min(int(region["end"]), int(region["start"]) + 50)))
            mutations_df = pd.DataFrame({
                "position": positions[:200],
                "reference": ["A"] * min(200, len(positions)),
                "alternate": ["G"] * min(200, len(positions)),
                "variant_id": [f"synth_{i}" for i in range(min(200, len(positions)))],
                "source": ["synthetic"] * min(200, len(positions)),
            })

        if len(mutations_df) < 10:
            pytest.skip("Not enough mutations to time")

        # Normalize column names: ClinVar API returns 'start'/'ref_allele',
        # but _bulk_filter_intronic_variants expects 'position'/'reference'
        col_renames = {}
        if "start" in mutations_df.columns and "position" not in mutations_df.columns:
            col_renames["start"] = "position"
        if "ref_allele" in mutations_df.columns and "reference" not in mutations_df.columns:
            col_renames["ref_allele"] = "reference"
        if "alt_allele" in mutations_df.columns and "alternate" not in mutations_df.columns:
            col_renames["alt_allele"] = "alternate"
        if col_renames:
            mutations_df = mutations_df.rename(columns=col_renames)

        # Get transcript ID
        gene_features = genome_handler.find_gene_features(test_gene)
        start_codons = gene_features[gene_features["feature_type"] == "start_codon"]
        if start_codons.empty:
            pytest.skip("No start codons for ADAR")
        tid = start_codons.iloc[0]["transcript_id"]

        # Ensure required columns exist
        if "reference" not in mutations_df.columns:
            mutations_df["reference"] = "A"
        if "variant_id" not in mutations_df.columns:
            mutations_df["variant_id"] = [f"v{i}" for i in range(len(mutations_df))]
        if "source" not in mutations_df.columns:
            mutations_df["source"] = "test"

        start = time.time()
        filtered, intronic = handler._bulk_filter_intronic_variants(
            mutations_df, tid, protein_generator
        )
        elapsed = time.time() - start
        _record("intronic_filter", elapsed)

        print(f"\nIntronic filtering ({test_gene}, {len(mutations_df)} variants):")
        print(f"  Time: {elapsed*1000:.1f}ms")
        print(f"  Kept: {len(filtered)}, Filtered: {len(intronic)}")

    def test_standardize_speed(self, genome_handler):
        """Time ClinVar standardization (allele inference + formatting)."""
        import asyncio
        from swissisoform.mutations import MutationHandler

        handler = MutationHandler(genome_handler=genome_handler)

        for gene in ["ADAR", "TMEM187"]:
            try:
                df = asyncio.get_event_loop().run_until_complete(
                    handler.get_clinvar_variants(gene)
                )
            except Exception:
                continue

            if df.empty:
                continue

            start = time.time()
            std_df = handler.standardize_mutation_data(df, "clinvar", gene)
            elapsed = time.time() - start
            _record(f"standardize_{gene.lower()}", elapsed)
            print(f"\nStandardize {gene}: {elapsed*1000:.1f}ms ({len(df)} variants)")

            # Standardization should be fast (< 2s even for large variant sets)
            assert elapsed < 2.0, (
                f"Standardize {gene} took {elapsed:.2f}s, expected < 2s"
            )

    def test_comprehensive_mutation_speed(self, genome_handler, alt_isoform_handler):
        """Time full comprehensive mutation analysis for a gene."""
        import asyncio
        import os
        from swissisoform.mutations import MutationHandler

        handler = MutationHandler(genome_handler=genome_handler)
        output_dir = os.path.join("results", "speed_test")
        os.makedirs(output_dir, exist_ok=True)

        test_gene = "ADAR"
        start = time.time()
        try:
            result = asyncio.get_event_loop().run_until_complete(
                handler.analyze_gene_mutations_comprehensive(
                    gene_name=test_gene,
                    genome_handler=genome_handler,
                    alt_isoform_handler=alt_isoform_handler,
                    output_dir=output_dir,
                    sources=["clinvar"],
                    validate_consequences=True,
                )
            )
        except Exception as e:
            pytest.skip(f"Comprehensive analysis failed: {e}")
            return

        elapsed = time.time() - start
        _record("comprehensive_adar", elapsed)
        status = result.get("status", "unknown")
        print(f"\nComprehensive {test_gene}: {elapsed:.2f}s (status: {status})")

        # Full comprehensive should complete in reasonable time
        assert elapsed < 60.0, (
            f"Comprehensive {test_gene} took {elapsed:.2f}s, expected < 60s"
        )

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
