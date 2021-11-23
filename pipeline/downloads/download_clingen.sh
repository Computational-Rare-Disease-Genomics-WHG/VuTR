#!/bin/sh

cd $(dirname "$0")

cd ../../

# Download the ClinGen's curated gene list
curl \
 "ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv" \
 --output data/pipeline/CLINGEN/ClinGen_gene_curation_list_GRCh38.tsv