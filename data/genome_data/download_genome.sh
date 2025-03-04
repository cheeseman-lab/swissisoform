#!/bin/bash

#!/bin/bash
GENCODE_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25"
wget ${GENCODE_BASE}/GRCh38.p7.genome.fa.gz
wget ${GENCODE_BASE}/gencode.v25.annotation.gtf.gz

GENCODE_BASE_LATEST="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release"
wget ${GENCODE_BASE_LATEST}/gencode.v47.annotation.gtf.gz

gunzip *.gz
rm *.gz
