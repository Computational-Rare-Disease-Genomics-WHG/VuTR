#!/bin/sh

cd $(dirname "$0")

cd ..

mkdir -p data/pipeline/vep_data


# Create the vep_data pipeline
mkdir -p data/pipeline/vep_data/input
mkdir -p data/pipeline/vep_data/output
mkdir -p data/pipeline/vep_data/Plugins

# Copy UTR Annotator source code to Plugins
