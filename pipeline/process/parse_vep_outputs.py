"""Parses the outputs from VEP running UTR annotator and saves it as a .tsv"""

import pandas as pd
import numpy as np
import argparse


def wide_to_long(df):
    """
    Make from wide to long. 
    @param df with 5'utr consequence concatentated by &.
    @returns long_df What we want is a 
            row per variant / transcript / five_prime_UTR_variant_annotation. 
    """
    long_df = pd.DataFrame(columns=df.columns.values.tolist())

    for index, row in df.iterrows():

        binding_df = pd.DataFrame()
        # If single just row bind it to long_df
        if "&" not in row["five_prime_UTR_variant_consequence"]:
            binding_df = pd.DataFrame(row.to_frame().T)
        else:
            # Replicate the rows per ORF consequences
            nrows = row["five_prime_UTR_variant_consequence"].count("&")
            consequences_split = row["five_prime_UTR_variant_consequence"].split("&")
            annotation_split = row["five_prime_UTR_variant_annotation"].split("&")
            binding_df = pd.concat([row.to_frame().T.drop(
                columns=["five_prime_UTR_variant_consequence",
                         "five_prime_UTR_variant_annotation"])]*(nrows+1),
                ignore_index=True)
            binding_df["five_prime_UTR_variant_consequence"] = consequences_split
            binding_df["five_prime_UTR_variant_annotation"] = annotation_split
        long_df = pd.concat([long_df, binding_df], axis=0, ignore_index=True)

    return long_df


def main(args):
    """
    Main entry point
    """
    vep_file_path = args.vep_file
    mane_summary_path = '../../data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.summary.txt.gz'
    write_path = args.output_file

    # Read clinvar file and the mane summary file
    vep_df = pd.read_csv(vep_file_path,
                         sep='\t',
                         skiprows=args.header_lines)

    mane_summary_df = pd.read_csv(mane_summary_path,
                                  sep='\t')

    # Filter to variants that impact uORFs (remove all other variants)
    no_utr_consequence = ['-', '', None, np.nan]
    vep_df = vep_df[~vep_df['five_prime_UTR_variant_annotation'].isin(
        no_utr_consequence)]

    # Remove version number
    mane_summary_df['transcript_id'] = mane_summary_df['Ensembl_nuc'].apply(
        lambda x: x[0:15])

    # Filter to consequence on the MANE transcript
    vep_df = vep_df[vep_df['Feature'].isin(
        mane_summary_df['transcript_id'])]

    # Convert wide to long data frame
    long_vep_df = wide_to_long(vep_df)

    # Write to file
    long_vep_df.to_csv(write_path, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Creates a TSV file that creates all possible UTR variants'
    )
    parser.add_argument(
        '--header_lines',
        required=True,
        type=int,
        help='Number of lines to skip when reading in the VEP file'
    )
    parser.add_argument(
        '--vep_file',
        required=True,
        type=str,
        help='Input file from VEP + UTR Annotator',
    )
    parser.add_argument(
        '--output_file',
        required=True,
        type=str,
        help='Output file name and location',
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose outputs',
    )
    main(args=parser.parse_args())
