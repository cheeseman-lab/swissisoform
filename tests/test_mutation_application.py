"""Tests for mutation application — uses synthetic mutations, no API calls."""

import pytest
from Bio.Seq import Seq


class TestMutationApplication:
    """Test mutation application produces valid sequences."""

    def _create_mock_mutation(self, genomic_pos, ref, alt):
        """Create a minimal mutation-like Series for testing."""
        import pandas as pd
        return pd.Series({
            "position": genomic_pos,
            "reference": ref,
            "alternate": alt,
            "impact": "missense variant",
            "hgvsc": f"c.1{ref}>{alt}",
            "variant_id": f"test_{genomic_pos}_{ref}>{alt}",
            "source": "test",
        })

    def test_mutation_preserves_coding_length_mod3(
        self, protein_generator, genome_handler
    ):
        """After applying a single-base substitution, coding length should stay mod 3."""
        # Use ADAR as a well-known gene
        # First get canonical to find a valid CDS position
        gene_features = genome_handler.find_gene_features("ADAR")
        if gene_features.empty:
            pytest.skip("ADAR not found in annotations")

        start_codons = gene_features[gene_features["feature_type"] == "start_codon"]
        if start_codons.empty:
            pytest.skip("No start codons for ADAR")

        tid = start_codons.iloc[0]["transcript_id"]
        canonical = protein_generator.extract_canonical_protein(tid)
        if canonical is None:
            pytest.skip("Could not extract canonical ADAR")

        coding_seq = canonical["coding_sequence"]
        assert len(coding_seq) % 3 == 0, "Canonical should be mod 3"

        # Find a position in the CDS to mutate
        cds = gene_features[
            (gene_features["feature_type"] == "CDS")
            & (gene_features["transcript_id"] == tid)
        ]
        if cds.empty:
            pytest.skip("No CDS for ADAR")

        # Pick a position in the middle of the first CDS
        first_cds = cds.iloc[0]
        strand = canonical["strand"]
        test_pos = int(first_cds["start"]) + 10

        # Get the actual base at that position
        ref_base = str(genome_handler.get_sequence(
            first_cds["chromosome"], test_pos, test_pos, "+"
        ))

        # Substitute with a different base
        alt_base = {"A": "T", "T": "A", "G": "C", "C": "G"}[ref_base]

        result, error = protein_generator._apply_mutation_to_sequence(
            tid, test_pos, ref_base, alt_base, "canonical"
        )

        if result is None:
            # Synonymous or other skip reason — that's okay
            pytest.skip(f"Mutation was synonymous or skipped: {error}")

        mutated_seq = result["coding_sequence"]
        assert len(mutated_seq) % 3 == 0, (
            f"Mutated coding sequence length {len(mutated_seq)} not divisible by 3"
        )

    def test_single_base_substitution_changes_one_base(
        self, protein_generator, genome_handler
    ):
        """A single-base substitution should change exactly one base in CDS."""
        gene_features = genome_handler.find_gene_features("ADAR")
        start_codons = gene_features[gene_features["feature_type"] == "start_codon"]
        if start_codons.empty:
            pytest.skip("No start codons for ADAR")

        tid = start_codons.iloc[0]["transcript_id"]
        canonical = protein_generator.extract_canonical_protein(tid)
        if canonical is None:
            pytest.skip("Could not extract canonical ADAR")

        cds = gene_features[
            (gene_features["feature_type"] == "CDS")
            & (gene_features["transcript_id"] == tid)
        ]
        first_cds = cds.iloc[0]
        test_pos = int(first_cds["start"]) + 10

        ref_base = str(genome_handler.get_sequence(
            first_cds["chromosome"], test_pos, test_pos, "+"
        ))
        alt_base = {"A": "T", "T": "A", "G": "C", "C": "G"}[ref_base]

        result, error = protein_generator._apply_mutation_to_sequence(
            tid, test_pos, ref_base, alt_base, "canonical"
        )

        if result is None:
            pytest.skip(f"Mutation skipped: {error}")

        orig_seq = canonical["coding_sequence"]
        mut_seq = result["coding_sequence"]

        # Count differences
        differences = sum(
            1 for a, b in zip(orig_seq, mut_seq) if a != b
        )
        assert differences == 1, (
            f"Expected exactly 1 base change, got {differences}"
        )
