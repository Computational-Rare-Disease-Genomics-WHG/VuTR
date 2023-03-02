#!/bin/sh

# The following script subsets ClinVar and runs UTR annotator
# on 5 prime variants specifically on the BMRC

module load BCFtools

REGIONS_FILE='/well/whiffin/users/wdf915/Projects/UTR-Visualisation-App/data/pipeline'
CLINVAR_DATA_PATH='/well/whiffin/users/wdf915/Projects/UTR-Visualisation-App/data/pipeline/CLINVAR'
VEP_CACHE='/well/whiffin/projects/vep'
ASSEMBLY='GRCh38'
MANE_VERSION='1.0'

tabix -H ${CLINVAR_DATA_PATH}/clinvar.vcf.gz > 5putr.vcf
tabix ${CLINVAR_DATA_PATH}/clinvar.vcf.gz \
            --regions ${REGIONS_FILE}/UTR_regions.tsv \
            >> ${CLINVAR_DATA_PATH}/5putr.vcf

module purge 
module load VEP/103.1-GCC-10.2.0 # Specify version 

# Run VEP
vep \
    --assembly GRCh38\
    --force_overwrite\
    --species homo_sapiens\
    --cache \
    --fasta ${VEP_CACHE}/Homo_sapiens.GRCh38.dna.primary_assembly.fa\
    --mane \
    --mane_select \
    --canonical \
    --offline\
    --tab \
    --fork 10 \
    --dir_cache ${VEP_CACHE} \
    --plugin UTRannotator\
    -i ${CLINVAR_DATA_PATH}/5putr.vcf\
    -o ${CLINVAR_DATA_PATH}/5putr.annotated.txt 