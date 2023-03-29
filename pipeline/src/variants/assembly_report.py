"""
assembly_report.py

This script takes the assembly report from NCBI and creates a VCF header file
and a tab delimited file with the contigs and their lengths.

python3 assembly_report.py 
    --input ../../../data/pipeline/GCA_000001405.15_GRCh38_assembly_report.txt
    --output assembly_report.txt 
"""

import argparse
import pandas


def main(args):
    """
    Main function
    """
    b = pandas.read_csv(
        args.input,
        sep="\t",
        skiprows=61)
    b = b[['# Sequence-Name', 'RefSeq-Accn',
           'UCSC-style-name', 'Assembly-Unit',
           'Sequence-Length']]
    b = b[b['Assembly-Unit'] == "Primary Assembly"]
    b = b[b['# Sequence-Name'].isin([str(i) for i in list(range(1, 23))+["X", "Y"]])]
    b.to_csv(args.output, index=False, sep="\t")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=""  # noqa: E501 # pylint: disable=C0301
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Input File",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output file",
    )
    parser.add_argument(
        "--header_output",
        type=str,
        help="Header for the VCF file",
    )
    main(args=parser.parse_args())
