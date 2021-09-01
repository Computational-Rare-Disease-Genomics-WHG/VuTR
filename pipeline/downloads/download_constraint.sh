#!/bin/bash

# downloads constraint metrics from gnomad's public bucket
# ensure that gsutil is installed

# Currently using exomes 2.1.1

# TODO : Update when the new constraint metrics will be out

cd $(dirname "$0")

cd ../../

gs_contraint_transcript="gs://gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz"
gs_contraint_gene="gs://gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

gsutil cp gs://gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz data/pipeline/GNOMAD
gsutil cp gs://gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz data/pipeline/GNOMAD
