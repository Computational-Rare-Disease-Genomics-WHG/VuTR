#!/bin/sh

curl "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar.vcf.gz" \
    --output input_vep/clinvar.vcf.gz