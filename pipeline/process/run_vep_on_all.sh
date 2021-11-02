#!/bin/bash

# Runs VEP + UTR Annotator on all possible variants
# generated by generate_all_possible_variants.py

INPATH='../../data/pipline/vep_data/input'
OUTPATH='../../data/pipline/vep_data/output'
VEP_CACHE=''
ASSEMBLY='GRCh38'
MANE_VERSION='0.93'


for i in {1..22} X Y
do
    echo $i

    vep \
    --assembly GRCh38\
    --force_overwrite\
    --species homo_sapiens\
    --cache \
    --mane \
    --mane_select\
    --canonical\
    --offline\
    --tab \
    --dir_cache /opt/vep/.vep/Cache \
    --plugin UTRannotator, /opt/vep/.vep/Plugins/uORF_starts_ends_GRCh38_PUBLIC.txt\
    -i ${INPATH}/UTR_variants_all_possible_{ASSEMBLY}_{MANE_VERSION}_{i}.txt\
    -o ${OUTPATH}/UTR_variants_vep_all_possible_{ASSEMBLY}_{MANE_VERSION}_{i}.txt
done
