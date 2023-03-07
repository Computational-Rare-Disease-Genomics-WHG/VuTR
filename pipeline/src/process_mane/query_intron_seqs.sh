#!/bin/bash

# To be added at a later stage to nextflow
# Define default values for the parameters
bed_file=""
reference_fasta=""
output_file=""

# Parse named arguments
while [ $# -gt 0 ]; do
    case "$1" in
        --bed-file=*)
            bed_file="${1#*=}"
            ;;
        --reference-fasta=*)
            reference_fasta="${1#*=}"
            ;;
        --output-file=*)
            output_file="${1#*=}"
            ;;
        *)
            echo "Invalid argument: $1"
            exit 1
            ;;
    esac
    shift
done

# Check that the required arguments are set
if [ -z "$bed_file" ] || [ -z "$reference_fasta" ]; then
    echo "Usage: $0 --bed-file=<bed_file> --reference-fasta=<reference_fasta> [--output-file=<output_file>]"
    exit 1
fi

# Check if the input files exist
if [ ! -f "$bed_file" ]; then
    echo "Bed file not found: $bed_file"
    exit 1
fi

if [ ! -f "$reference_fasta" ]; then
    echo "Reference fasta file not found: $reference_fasta"
    exit 1
fi

# Define the output file name if not provided
if [ -z "$output_file" ]; then
    echo "Output fasta file not found: $output_file"
    exit 1
fi

# Run bedtools getfasta with the input files
bedtools getfasta -fi "$reference_fasta" -bed "$bed_file" -fo "$output_file"