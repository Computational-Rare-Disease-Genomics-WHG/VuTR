"""
assembly_report.py

This script takes the assembly report from NCBI and creates a VCF header file
and a tab delimited file with the contigs and their lengths.
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

    # Create a list of header lines for the contigs
    contig_lines = list(b.apply(
     lambda x: f'##contig=<ID={x["# Sequence-Name"]},length={x["Sequence-Length"]}>\n',
     axis=1
    ))
    # Create a VCF header file
    with open(args.header_output, "w", encoding='utf-8') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##fileDate=20190801\n")
        f.write("##source=NCBI\n")
        f.write("##reference=GRCh38\n")
        f.writelines(contig_lines)


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
