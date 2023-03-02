#!/bin/bash 

cd $(dirname "$0")

cd ..

mkdir -p data/pipeline/vep_data

# Create the vep_data pipeline
mkdir -p data/pipeline/vep_data/input
mkdir -p data/pipeline/vep_data/output
mkdir -p data/pipeline/vep_data/Plugins

# Copy UTR Annotator source code to Plugins
cp -R -u -p pipeline/UTR-Annotator/* data/pipeline/vep_data/Plugins

# Directories
mkdir -p data/pipeline/MANE

mkdir -p data/pipeline/CLINVAR
mkdir -p data/pipeline/CLINGEN
mkdir -p data/pipeline/GNOMAD
mkdir -p data/pipeline/REFSEQ
mkdir -p data/pipeline/SMORFS


mane_v=${mane_v:-1.0}

while [ $# -gt 0 ]; do

   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        # echo $1 $2 // Optional to see the parameter:value result
   fi

  shift
done

# Download MANE 
bash ./downloads/mane.sh --mane_v ${mane_v}

# Download Clingen
bash ./downloads/clingen.sh

# Download Clinvar 
bash ./downloads/clinvar.sh

# Download Constraint
bash ./downloads/constraint.sh

# Download Omim
bash ./downloads/omim.sh

# Download Reference Sequence
bash ./downloads/refseq.sh

# Download translational efficiency
bash ./downloads/te.sh


echo "Download completed. See ../../data/pipeline/* for the output"



# We aren't using sORFS.org anymore but mentioned here for archival purposes
# mkdir -p data/pipeline/SORFS
# # Download sORFS
# python3 downloads/download_sorf_query.py
