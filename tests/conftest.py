"""Shared fixtures for SwissIsoform tests.

Loads genome and annotation data once per session since they are expensive.
"""

import pytest
import time
import logging

from swissisoform.genome import GenomeHandler
from swissisoform.alternative_isoforms import AlternativeIsoform
from swissisoform.translation import AlternativeProteinGenerator

logger = logging.getLogger(__name__)

# Data paths
GENOME_PATH = "data/genome_data/GRCh38.p7.genome.fa"
GTF_PATH = "data/genome_data/gencode.v25.annotation.v47names.gtf"
BED_PATH = "data/ribosome_profiling/hela_isoforms_with_transcripts.bed"

# Test genes
KNOWN_GOOD_GENES = ["ADAR", "DYNLRB1", "RBM10", "TCP1"]
KNOWN_FAILING_GENES = ["TMEM187", "ARL8B", "TMEM205"]

# Known failing extension cases
KNOWN_FAILING_EXTENSIONS = {
    "TMEM187": {"transcript": "ENST00000369982.4", "codon": "CTG", "pos": 153982018},
    "ARL8B": {"transcript": "ENST00000611208.4", "codon": "GTG", "pos": 5122403},
    "TMEM205": {"transcript": "ENST00000586218.5", "codon": "ATT", "pos": 11346148},
}

ALL_TEST_GENES = KNOWN_GOOD_GENES + KNOWN_FAILING_GENES


@pytest.fixture(scope="session")
def genome_handler():
    """Load genome and annotations once per session."""
    start = time.time()
    handler = GenomeHandler(GENOME_PATH, GTF_PATH)
    elapsed = time.time() - start
    logger.info(f"Genome loaded in {elapsed:.1f}s")
    return handler


@pytest.fixture(scope="session")
def alt_isoform_handler():
    """Load BED file once per session."""
    start = time.time()
    handler = AlternativeIsoform(debug=False)
    handler.load_bed(BED_PATH)
    elapsed = time.time() - start
    logger.info(f"BED loaded in {elapsed:.1f}s")
    return handler


@pytest.fixture(scope="session")
def protein_generator(genome_handler, alt_isoform_handler, tmp_path_factory):
    """Create protein generator once per session."""
    output_dir = str(tmp_path_factory.mktemp("test_output"))
    return AlternativeProteinGenerator(
        genome_handler=genome_handler,
        alt_isoform_handler=alt_isoform_handler,
        output_dir=output_dir,
        debug=True,
    )


@pytest.fixture(scope="session")
def timing_data():
    """Shared dict to accumulate timing data across tests."""
    return {}
