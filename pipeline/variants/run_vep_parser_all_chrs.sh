#!/bin/bash

# Flag : To be replaced by worflow management
# Runs the parser on all possible variants

INPATH='../../data/pipeline/vep_data/output'
OUTPATH='../../data/pipeline/vep_data/output'
VEP_CACHE='/well/whiffin/shared/vep'
ASSEMBLY='GRCh38'
MANE_VERSION='0.93'
HEADER_LINES=47

for i in {1..22} X Y
do
    echo "Running parser on chr ${i}"
    python3 parse_vep_outputs.py  \
    --vep_file ${INPATH}/UTR_variants_vep_all_possible_GRCh38_0.93_chr${i}.txt \
    --output_file ${OUTPATH}/UTR_variants_all_possible_GRCh38_0.93_chr_${i}_parsed_vep.txt \
    --header_lines $HEADER_LINES
done