"""Tests for truncation protein extraction correctness."""

import warnings
import pytest
from tests.conftest import ALL_TEST_GENES


class TestTruncationProtein:
    """Test truncation protein extraction."""

    def _get_truncation_features(self, alt_isoform_handler, gene_name):
        """Get truncation features for a gene."""
        features = alt_isoform_handler.get_translation_features(gene_name)
        if features.empty:
            return features
        return features[features["region_type"] == "truncation"]

    @pytest.mark.parametrize("gene_name", ALL_TEST_GENES)
    def test_truncation_coding_length_divisible_by_3(
        self, protein_generator, alt_isoform_handler, gene_name
    ):
        """Truncation coding sequence length MUST be divisible by 3 — NO TRIMMING."""
        features = self._get_truncation_features(alt_isoform_handler, gene_name)
        if features.empty:
            pytest.skip(f"No truncation features for {gene_name}")

        failures = []
        for idx, feature in features.iterrows():
            tid = feature["transcript_id"]
            result = protein_generator.extract_alternative_protein(tid, feature)
            if result is None:
                continue

            coding_seq = result["coding_sequence"]
            if len(coding_seq) % 3 != 0:
                failures.append(
                    f"{gene_name} ({tid}): length {len(coding_seq)}, "
                    f"remainder {len(coding_seq) % 3}, "
                    f"alt_start={result.get('alternative_start_pos', '?')}"
                )

        assert not failures, (
            f"Truncation sequences not divisible by 3:\n"
            + "\n".join(failures)
        )

    @pytest.mark.parametrize("gene_name", ALL_TEST_GENES)
    def test_truncation_no_biopython_warnings(
        self, protein_generator, alt_isoform_handler, gene_name
    ):
        """No BiopythonWarning should be raised during truncation translation."""
        features = self._get_truncation_features(alt_isoform_handler, gene_name)
        if features.empty:
            pytest.skip(f"No truncation features for {gene_name}")

        for idx, feature in features.iterrows():
            tid = feature["transcript_id"]
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = protein_generator.extract_alternative_protein(tid, feature)

                bio_warnings = [
                    x for x in w
                    if "BiopythonWarning" in str(type(x.category))
                    or "Partial codon" in str(x.message)
                ]
                assert len(bio_warnings) == 0, (
                    f"{gene_name} ({tid}): BiopythonWarning during truncation: "
                    f"{[str(x.message) for x in bio_warnings]}"
                )

    @pytest.mark.parametrize("gene_name", ALL_TEST_GENES)
    def test_truncation_shorter_than_canonical(
        self, protein_generator, alt_isoform_handler, genome_handler, gene_name
    ):
        """Truncated protein should be shorter than canonical."""
        features = self._get_truncation_features(alt_isoform_handler, gene_name)
        if features.empty:
            pytest.skip(f"No truncation features for {gene_name}")

        for idx, feature in features.iterrows():
            tid = feature["transcript_id"]
            canonical = protein_generator.extract_canonical_protein(tid)
            if canonical is None:
                continue

            result = protein_generator.extract_alternative_protein(tid, feature)
            if result is None:
                continue

            assert len(result["protein"]) < len(canonical["protein"]), (
                f"{gene_name} ({tid}): truncated protein ({len(result['protein'])} AA) "
                f"is not shorter than canonical ({len(canonical['protein'])} AA)"
            )
