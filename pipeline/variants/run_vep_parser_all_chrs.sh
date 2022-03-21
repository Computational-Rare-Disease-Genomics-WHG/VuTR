#!/bin/bash

# Runs the parser on all possible variants
# This splits each row of the annotated VEP table
# to a single row per variant consequence
# Also merges the chrom_separated table to a single table

# Flag : Hardcoded paths and version names
# to be replaced by worflow management


INPATH='../../data/pipeline/vep_data/output'
OUTPATH='../../data/pipeline/vep_data/output'
PIPLINE_OUTPATH='../../data/pipeline'
VEP_CACHE='/well/whiffin/shared/vep'
ASSEMBLY='GRCh38'
MANE_VERSION='1.0'
HEADER_LINES=47

for i in {1..22} X Y
do
    echo "Running parser on chr ${i}"
    python3 parse_vep_outputs.py  \
    --vep_file ${INPATH}/UTR_variants_vep_all_possible_${ASSEMBLY}_${MANE_VERSION}_chr${i}.txt \
    --output_file ${OUTPATH}/UTR_variants_all_possible_${ASSEMBLY}_${MANE_VERSION}_chr${i}_parsed_vep.txt \
    --header_lines $HEADER_LINES
done

# Merge all of the variants into a single file
head -n 1 ${OUTPATH}/UTR_variants_all_possible_${ASSEMBLY}_${MANE_VERSION}_chr1_parsed_vep.txt > \
 ${PIPLINE_OUTPATH}/UTR_variants_all_possible_${ASSEMBLY}_${MANE_VERSION}_all_chrs.tsv && \
 tail -n+2 -q ${OUTPATH}/UTR_variants_all_possible_${ASSEMBLY}_${MANE_VERSION}_chr*_parsed_vep.txt >> \
  ${PIPLINE_OUTPATH}/UTR_variants_all_possible_${ASSEMBLY}_${MANE_VERSION}_all_chrs.tsv
