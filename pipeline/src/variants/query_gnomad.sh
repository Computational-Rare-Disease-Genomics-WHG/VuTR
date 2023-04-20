#!/bin/bash

# Check if the input file argument was provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

# Get the input file path from the command-line argument
input=$1

# Check if the input file exists
if [ ! -f $input ]; then
    echo "Input file not found: $input"
    exit 1
fi

# Get the output file path from the command-line argument
output=$2

# Use bcftools query to filter the VCF file for the specified fields
bcftools query -f '%CHROM\t%POS\t%ID\t%HGVSC\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' $input > $output

echo "Filtered VCF file saved to $output"
