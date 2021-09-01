#!/bin/bash

# build the index file
tabix -p vcf vep_data/input/clinvar.vcf.gz -f
# filter to regions
tabix vep_data/input/clinvar.vcf.gz -R UTR_regions.tsv > vep_data/input/clinvar_utr_filtered.vcf
