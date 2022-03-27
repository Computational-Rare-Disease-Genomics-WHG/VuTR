#!/bin/bash

# To run
# bash download_mane --mane_v 0.95

mane_v=${mane_v:-0.93}

while [ $# -gt 0 ]; do

   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        # echo $1 $2 // Optional to see the parameter:value result
   fi

  shift
done

# Base url
ftp_url="ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${mane_v}/*"

# CD to script directory
cd $(dirname "$0")

# Create MANE directory if it doesn't exist
mkdir -p ../../data/pipeline/MANE/${mane_v}

# Download to directory
wget --no-parent $ftp_url -P ../../data/pipeline/MANE/${mane_v}

echo "Completed downloading MANE ${mane_v}"
