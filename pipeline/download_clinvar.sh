#!/bin/sh

# Download the ClinVar's weekly release

curl "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar.vcf.gz" \
    --output vep_data/input/clinvar.vcf.gz
