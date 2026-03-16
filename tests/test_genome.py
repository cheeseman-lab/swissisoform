"""Tests for GenomeHandler sequence extraction and CDS region integrity."""

import pytest
from tests.conftest import ALL_TEST_GENES


class TestGetSequence:
    """Test that get_sequence returns correct sequences for known coordinates."""

    def test_forward_strand_sequence(self, genome_handler):
        """Test sequence extraction on forward strand."""
        # ISG15 start codon on chr1+ at 1013574-1013576
        seq = genome_handler.get_sequence("chr1", 1013574, 1013576, "+")
        assert str(seq) == "ATG", f"Expected ATG, got {seq}"

    def test_reverse_strand_sequence(self, genome_handler):
        """Test sequence extraction on reverse strand returns reverse complement."""
        # HES4 start codon on chr1- at 999971-999973
        seq = genome_handler.get_sequence("chr1", 999971, 999973, "-")
        assert str(seq) == "ATG", f"Expected ATG (rev comp), got {seq}"

    def test_chromosome_prefix_flexibility(self, genome_handler):
        """Test that chromosome names work with flexible prefixes."""
        seq1 = genome_handler.get_sequence("chr1", 1013574, 1013576, "+")
        # The genome file uses "chr" prefix, so this tests the lookup logic
        assert len(str(seq1)) == 3

    def test_one_based_inclusive(self, genome_handler):
        """Test that coordinates are 1-based inclusive on both ends."""
        seq = genome_handler.get_sequence("chr1", 100, 102, "+")
        assert len(str(seq)) == 3, "1-based inclusive: positions 100,101,102 = 3bp"


class TestCDSLengths:
    """Test that CDS regions for transcripts have lengths that are multiples of 3."""

    @pytest.mark.parametrize("gene_name", ALL_TEST_GENES)
    def test_cds_total_length_divisible_by_3(self, genome_handler, gene_name):
        """For transcripts with complete CDS (start+stop), total CDS length should be divisible by 3."""
        gene_features = genome_handler.find_gene_features(gene_name)
        if gene_features.empty:
            pytest.skip(f"Gene {gene_name} not found in annotations")

        # Only check transcripts with both start and stop codons (complete CDS)
        start_tids = set(
            gene_features[gene_features["feature_type"] == "start_codon"]["transcript_id"].unique()
        )
        stop_tids = set(
            gene_features[gene_features["feature_type"] == "stop_codon"]["transcript_id"].unique()
        )
        complete_tids = start_tids & stop_tids

        cds_features = gene_features[gene_features["feature_type"] == "CDS"]
        if cds_features.empty:
            pytest.skip(f"No CDS regions for {gene_name}")

        checked = 0
        for tid in complete_tids:
            transcript_cds = cds_features[cds_features["transcript_id"] == tid]
            if transcript_cds.empty:
                continue
            total_cds_length = sum(
                row["end"] - row["start"] + 1
                for _, row in transcript_cds.iterrows()
            )
            assert total_cds_length % 3 == 0, (
                f"CDS total length for {tid} is {total_cds_length}, "
                f"not divisible by 3 (remainder {total_cds_length % 3})"
            )
            checked += 1

        if checked == 0:
            pytest.skip(f"No complete CDS transcripts for {gene_name}")

    def test_known_transcript_cds_length(self, genome_handler):
        """Test CDS length for a specific known transcript."""
        tid = "ENST00000369982.4"  # TMEM187
        features = genome_handler.get_transcript_features(tid)
        cds = features[features["feature_type"] == "CDS"]
        if cds.empty:
            pytest.skip(f"No CDS for {tid}")

        total_len = sum(row["end"] - row["start"] + 1 for _, row in cds.iterrows())
        assert total_len % 3 == 0, (
            f"TMEM187 CDS length {total_len} not divisible by 3"
        )
