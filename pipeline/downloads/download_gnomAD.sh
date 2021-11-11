#!/bin/sh


cd $(dirname "$0")

cd ../../

OUTPATH="/well/whiffin/projects/gnomAD/"

for i in {1..22} X Y
do
    echo "Downloading chr${i} from GCP Bucket"
    gsutil cp gs://gcp-public-data--gnomad/release/3.1.1/gnomad.genomes.v3.1.1.sites.chr${i}.vcf.bgz ${OUTPATH}
    gsutil cp gs://gcp-public-data--gnomad/release/3.1.1/gnomad.genomes.v3.1.1.sites.chr${i}.vcf.bgz.tbi ${OUTPATH}
done