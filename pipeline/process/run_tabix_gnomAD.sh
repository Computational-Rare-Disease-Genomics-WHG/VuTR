#!/bin/bash

cd $(dirname "$0")

cd ../..

GNOMAD_INPATH='/well/whiffin/projects/gnomAD'
GNOMAD_OUTPUT='data/pipeline/GNOMAD'

# Iterate through the chroms
for i in {1..22} X Y
do
    echo "Filtering gnomAD on chr${i}"

    # Tabix
    tabix -p vcf {GNOMAD_INPATH}/gnomad.genomes.v3.1.1.sites.chr${i}.vcf.bgz -f

    # Filter gnomAD to sites
    tabix ${GNOMAD_INPATH}/gnomad.genomes.v3.1.1.sites.chr${i}.vcf.bgz \
            --regions data/pipeline/UTR_regions.tsv \
            > ${GNOMAD_OUTPUT}/gnomad.genomes.v3.1.1.utr_sites.chr${i}.vcf
done