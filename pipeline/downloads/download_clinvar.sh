#!/bin/sh

cd $(dirname "$0")

cd ../../

# Download the ClinVar's weekly release
curl "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar.vcf.gz" \
    --output data/pipeline/CLINVAR/clinvar.vcf.gz
