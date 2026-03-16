"""Tests for canonical protein extraction correctness."""

import warnings
import pytest
from tests.conftest import ALL_TEST_GENES


class TestCanonicalProtein:
    """Test canonical protein extraction for all test genes."""

    def _get_transcripts_for_gene(self, genome_handler, gene_name):
        """Get transcript IDs that have both start and stop codons (complete CDS)."""
        gene_features = genome_handler.find_gene_features(gene_name)
        if gene_features.empty:
            return []

        start_codons = gene_features[gene_features["feature_type"] == "start_codon"]
        stop_codons = gene_features[gene_features["feature_type"] == "stop_codon"]

        start_tids = set(start_codons["transcript_id"].unique())
        stop_tids = set(stop_codons["transcript_id"].unique())

        # Only return transcripts with complete CDS (both start and stop codons)
        return list(start_tids & stop_tids)

    @pytest.mark.parametrize("gene_name", ALL_TEST_GENES)
    def test_coding_sequence_divisible_by_3(self, protein_generator, genome_handler, gene_name):
        """Canonical coding sequence length MUST be divisible by 3."""
        transcripts = self._get_transcripts_for_gene(genome_handler, gene_name)
        if not transcripts:
            pytest.skip(f"No transcripts with start codons for {gene_name}")

        for tid in transcripts:
            result = protein_generator.extract_canonical_protein(tid)
            if result is None:
                continue

            coding_seq = result["coding_sequence"]
            assert len(coding_seq) % 3 == 0, (
                f"{gene_name} ({tid}): coding sequence length {len(coding_seq)} "
                f"is not divisible by 3 (remainder {len(coding_seq) % 3}). "
                f"This indicates broken coordinate extraction."
            )

    @pytest.mark.parametrize("gene_name", ALL_TEST_GENES)
    def test_protein_starts_with_methionine(self, protein_generator, genome_handler, gene_name):
        """Canonical protein should start with M (methionine)."""
        transcripts = self._get_transcripts_for_gene(genome_handler, gene_name)
        if not transcripts:
            pytest.skip(f"No transcripts with start codons for {gene_name}")

        for tid in transcripts:
            result = protein_generator.extract_canonical_protein(tid)
            if result is None:
                continue

            protein = result["protein"]
            if protein:
                assert protein[0] == "M", (
                    f"{gene_name} ({tid}): protein starts with '{protein[0]}', "
                    f"expected 'M'. First codon: {result['coding_sequence'][:3]}"
                )

    @pytest.mark.parametrize("gene_name", ALL_TEST_GENES)
    def test_protein_ends_with_stop(self, protein_generator, genome_handler, gene_name):
        """Canonical protein should end with * (stop codon)."""
        transcripts = self._get_transcripts_for_gene(genome_handler, gene_name)
        if not transcripts:
            pytest.skip(f"No transcripts with start codons for {gene_name}")

        for tid in transcripts:
            result = protein_generator.extract_canonical_protein(tid)
            if result is None:
                continue

            protein = result["protein"]
            if protein:
                assert protein[-1] == "*", (
                    f"{gene_name} ({tid}): protein ends with '{protein[-1]}', "
                    f"expected '*'. Last codon: {result['coding_sequence'][-3:]}"
                )

    @pytest.mark.parametrize("gene_name", ALL_TEST_GENES)
    def test_no_biopython_warnings(self, protein_generator, genome_handler, gene_name):
        """No BiopythonWarning should be raised during canonical translation."""
        transcripts = self._get_transcripts_for_gene(genome_handler, gene_name)
        if not transcripts:
            pytest.skip(f"No transcripts with start codons for {gene_name}")

        for tid in transcripts:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = protein_generator.extract_canonical_protein(tid)

                bio_warnings = [
                    x for x in w
                    if "BiopythonWarning" in str(type(x.category))
                    or "Partial codon" in str(x.message)
                ]
                assert len(bio_warnings) == 0, (
                    f"{gene_name} ({tid}): BiopythonWarning raised during canonical "
                    f"translation: {[str(x.message) for x in bio_warnings]}"
                )
