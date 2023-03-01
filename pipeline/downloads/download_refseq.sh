#!/bin/bash

echo "Started download..."

url="ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz"

mkdir ../../data/pipeline/refseq

wget $url -P  ../../data/pipeline/refseq

echo "...completed downloadxw"