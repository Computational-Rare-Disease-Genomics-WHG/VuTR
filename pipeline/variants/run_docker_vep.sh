#!/bin/sh

# To run an interactive terminal inside docker
docker run -t -i \
 -v $(realpath ../../data/pipeline/vep_data):/opt/vep/.vep\
 ensemblorg/ensembl-vep \
 ./vep \
 --assembly GRCh38\
 --force_overwrite\
 --species homo_sapiens\
 --cache \
 --mane \
 --mane_select\
 --canonical\
 --tab \
 --dir_cache /opt/vep/.vep/Cache \
 -i /opt/vep/.vep/input/input_list.txt\
 -o /opt/vep/.vep/output/output_list.out
