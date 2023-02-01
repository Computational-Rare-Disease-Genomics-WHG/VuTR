#!/bin/bash 

# Downloads the omim map. 


wget "https://www.omim.org/static/omim/data/mim2gene.txt" 
mv mim2gene.txt ../../data/pipeline
