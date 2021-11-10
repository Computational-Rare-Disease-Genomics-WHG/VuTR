#!/bin/bash

cd $(dirname "$0")

cd ../..

# build the index file
tabix -p vcf data/pipeline/CLINVAR/clinvar.vcf.gz -f

# filter to regions
tabix data/pipeline/CLINVAR/clinvar.vcf.gz \
       	-R data/pipeline/UTR_regions.tsv \
       	> data/pipeline/clinvar_utr_filtered.vcf
