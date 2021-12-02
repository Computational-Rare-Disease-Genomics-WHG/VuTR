#!/bin/bash


# Download the translational efficiency from
# Noderer, William L., et al. "Quantitative analysis of mammalian translation initiation sites by FACS‐seq." Molecular systems biology 10.8 (2014): 748.


curl "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4299517/bin/msb0010-0748-SD3.txt" |\
 sed 1d |\
 tr [:upper:] [:lower:] |\
 sed "1d;s/u/t/g" > tmp.txt

echo $'context\tefficiency\tlower_bound\tupper_bound'|
 cat -  tmp.txt |\
 sed "s/^M/\r/g" >\
 ../../data/pipeline/translational_efficiency.txt

rm tmp.txt
