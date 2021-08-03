#!/bin/sh

# To run an interactive terminal inside docker 
docker run -t -i \
 -v $(pwd)/UTRannotator/:/opt/vep/.vep/Plugins \
 -v $(pwd)/input_vep:/input_vep \
 -v $(pwd)/output_vep:/output_vep \
 ensemblorg/ensembl-vep \
 ./vep -i /input_vep/test_grch38.vcf \
 --assembly GRCh38\
 --force_overwrite\
 --species homo_sapiens\
 --tab \
 --database \
 -plugin UTRannotator,\
 ./UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt\
 -o /output_vep/testrun_vep.test_grch38.out