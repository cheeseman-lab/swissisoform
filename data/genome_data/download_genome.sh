#!/bin/bash
UCSC_BASE="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips"

wget ${UCSC_BASE}/hg38.fa.gz
wget ${UCSC_BASE}/genes/hg38.ncbiRefSeq.gtf.gz

gunzip *.gz
