"""Tests for extension protein extraction — focuses on known-failing genes."""

import warnings
import pytest
from tests.conftest import KNOWN_FAILING_EXTENSIONS, KNOWN_GOOD_GENES


class TestExtensionProtein:
    """Test extension protein extraction for known-failing and good genes."""

    def _get_extension_features(self, alt_isoform_handler, gene_name):
        """Get extension features for a gene."""
        features = alt_isoform_handler.get_translation_features(gene_name)
        if features.empty:
            return features
        return features[features["region_type"] == "extension"]

    @pytest.mark.parametrize(
        "gene_name,info",
        list(KNOWN_FAILING_EXTENSIONS.items()),
        ids=list(KNOWN_FAILING_EXTENSIONS.keys()),
    )
    def test_known_failing_extension_coding_length(
        self, protein_generator, alt_isoform_handler, gene_name, info
    ):
        """Extension coding sequence length MUST be divisible by 3 — NO TRIMMING.

        These are the known-failing cases. If they fail, the coordinate logic
        is wrong and needs fixing.
        """
        features = self._get_extension_features(alt_isoform_handler, gene_name)
        if features.empty:
            pytest.skip(f"No extension features for {gene_name}")

        # Find the specific failing extension
        target_tid = info["transcript"]
        tid_features = features[features["transcript_id"] == target_tid]

        if tid_features.empty:
            pytest.skip(f"No extension features for {gene_name} transcript {target_tid}")

        for idx, feature in tid_features.iterrows():
            result = protein_generator.extract_alternative_protein(target_tid, feature)
            if result is None:
                pytest.skip(f"Extension extraction returned None for {gene_name}")

            coding_seq = result["coding_sequence"]
            ext_len = result.get("extension_sequence_length", 0)
            cds_len = result.get("cds_sequence_length", 0)

            assert len(coding_seq) % 3 == 0, (
                f"KNOWN FAILING: {gene_name} ({target_tid}) extension at {info['pos']}: "
                f"coding sequence length {len(coding_seq)} not divisible by 3 "
                f"(remainder {len(coding_seq) % 3}). "
                f"Extension part: {ext_len}bp, CDS part: {cds_len}bp. "
                f"Total = {ext_len + cds_len}bp. "
                f"DO NOT FIX BY TRIMMING — fix the coordinate logic."
            )

    @pytest.mark.parametrize(
        "gene_name,info",
        list(KNOWN_FAILING_EXTENSIONS.items()),
        ids=list(KNOWN_FAILING_EXTENSIONS.keys()),
    )
    def test_known_failing_no_biopython_warnings(
        self, protein_generator, alt_isoform_handler, gene_name, info
    ):
        """No BiopythonWarning should be raised for extension translations."""
        features = self._get_extension_features(alt_isoform_handler, gene_name)
        if features.empty:
            pytest.skip(f"No extension features for {gene_name}")

        target_tid = info["transcript"]
        tid_features = features[features["transcript_id"] == target_tid]
        if tid_features.empty:
            pytest.skip(f"No extension features for {target_tid}")

        for idx, feature in tid_features.iterrows():
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = protein_generator.extract_alternative_protein(target_tid, feature)

                bio_warnings = [
                    x for x in w
                    if "BiopythonWarning" in str(type(x.category))
                    or "Partial codon" in str(x.message)
                ]
                assert len(bio_warnings) == 0, (
                    f"{gene_name} ({target_tid}): BiopythonWarning during extension: "
                    f"{[str(x.message) for x in bio_warnings]}"
                )

    @pytest.mark.parametrize("gene_name", KNOWN_GOOD_GENES)
    def test_good_gene_extensions(
        self, protein_generator, alt_isoform_handler, gene_name
    ):
        """Extension proteins for known-good genes should also have valid lengths."""
        features = self._get_extension_features(alt_isoform_handler, gene_name)
        if features.empty:
            pytest.skip(f"No extension features for {gene_name}")

        for idx, feature in features.iterrows():
            tid = feature["transcript_id"]
            result = protein_generator.extract_alternative_protein(tid, feature)
            if result is None:
                continue

            coding_seq = result["coding_sequence"]
            assert len(coding_seq) % 3 == 0, (
                f"{gene_name} ({tid}): extension coding sequence length "
                f"{len(coding_seq)} not divisible by 3 (remainder {len(coding_seq) % 3})"
            )
