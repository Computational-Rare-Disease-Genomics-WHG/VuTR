#!/bin/sh

# To run an interactive terminal inside docker 
docker run -t -i \
 -v $(pwd)/vep_data:/opt/vep/.vep\
 -v /home/elston/Packages/vep_cache:/opt/vep/.vep/Cache\
 ensemblorg/ensembl-vep \
 ./vep \
 --assembly GRCh38\
 --force_overwrite\
 --species homo_sapiens\
 --cache \
 --offline\
 --tab \
 --dir_cache /opt/vep/.vep/Cache \
 --plugin UTRannotator, /opt/vep/.vep/Plugins/uORF_starts_ends_GRCh38_PUBLIC.txt\
 -i /opt/vep/.vep/input/test_grch38.vcf \
 -o /opt/vep/.vep/output/testrun_vep.test_grch38.out
