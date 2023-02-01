#!/bin/bash 


mane_v=${mane_v:-0.93}

while [ $# -gt 0 ]; do

   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        # echo $1 $2 // Optional to see the parameter:value result
   fi

  shift
done

# Download MANE 
bash ./download_mane.sh --mane_v ${mane_v}

# Download Clingen
bash ./download_clingen.sh

# Download Clinvar 
bash ./download_clinvar.sh

# Download Constraint
bash ./download_constraint.sh

# Download sORFS
python3 download_sorf_query.py

echo "Download completed. See ../../data/pipeline/* for the output"