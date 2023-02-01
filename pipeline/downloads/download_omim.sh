#!/bin/bash 

# Downloads the omim map. 
cd $(dirname "$0")

cd ../../

wget "https://www.omim.org/static/omim/data/mim2gene.txt" 
mv mim2gene.txt data/pipeline
