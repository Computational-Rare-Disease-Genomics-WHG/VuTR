"""

Parses the outputs from VEP running UTR annotator and saves it as a .tsv
Potential performance improvement : Use Dask as the dataframe are really big
"""
import sys
import os
import argparse
import tqdm
import pandas as pd
import numpy as np


def wide_to_long(df):
    """
    Convert a DataFrame from wide to long format.

    Args:
        df (pandas.DataFrame): The input DataFrame, with 5'UTR consequence concatenated by "&".

    Returns:
        pandas.DataFrame: A new DataFrame with one row per variant/transcript/five_prime_UTR_variant_annotation.
    """
    # Split the columns into lists of values
    df = df.copy()
    df['five_prime_UTR_variant_consequence'] = df['five_prime_UTR_variant_consequence'].str.split('&')
    df['five_prime_UTR_variant_annotation'] = df['five_prime_UTR_variant_annotation'].str.split('&')

    # Define a boolean mask for the rows with multiple consequences
    has_multiple_consequences = df['five_prime_UTR_variant_consequence'].apply(lambda x: len(x) > 1)

    # Define the subsets of the DataFrame
    single_consequence_df = df.loc[~has_multiple_consequences]
    multiple_consequence_df = df.loc[has_multiple_consequences]

    # Melt the columns with multiple consequences
    melted_df = multiple_consequence_df.explode(['five_prime_UTR_variant_consequence', 'five_prime_UTR_variant_annotation'])

    # Concatenate the DataFrames
    single_consequence_df['five_prime_UTR_variant_consequence'] = single_consequence_df['five_prime_UTR_variant_consequence'].apply(lambda x: x[0])
    single_consequence_df['five_prime_UTR_variant_annotation'] = single_consequence_df['five_prime_UTR_variant_annotation'].apply(lambda x: x[0])

    long_df = pd.concat([single_consequence_df, melted_df], ignore_index=True)

    return long_df


def main(args):
    """
    Main entry point
    """
    vep_file_path = args.vep_file
    mane_summary_path = args.mane_file
    write_path = args.output_file

    print(f'Reading file {vep_file_path}')
    # Read clinvar file and the mane summary file
    vep_df = pd.read_csv(vep_file_path, sep='\t', skiprows=args.header_lines)

    print('Read completed')

    mane_summary_df = pd.read_csv(mane_summary_path, sep='\t')

    print('Filtering out 5 prime UTR variants with no consequence')
    # Filter to variants that impact uORFs (remove all other variants)
    no_utr_consequence = ['-', '', None, np.nan]
    vep_df = vep_df[
        ~vep_df['five_prime_UTR_variant_annotation'].isin(no_utr_consequence)
    ]

    # Remove version number
    mane_summary_df['transcript_id'] = mane_summary_df['Ensembl_nuc'].apply(
        lambda x: x[0:15]
    )

    print('Filtering to MANE consequences')
    # Filter to consequence on the MANE transcript
    vep_df = vep_df[vep_df['Feature'].isin(mane_summary_df['transcript_id'])]

    # Chunk the data frame
    df_groups = vep_df.groupby(vep_df.index // args.chunk_size)

    if os.path.isfile(write_path):
        # file already exists
        print('File already exists. Aborting.')
        sys.exit(1)

    print(f'Writing to file {write_path}')

    # Convert wide to long data frame and write to file
    header = True
    for _group_idx, df_group in tqdm.tqdm(df_groups):
        wide_to_long(df_group).to_csv(write_path, sep='\t', mode='a',
                                      header=header, index=False)
        header = False


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parses the vep output files')
    parser.add_argument(
        '--mane_file',
        required=True,
        type=str,
        help='Which mane version to use?',
    )
    parser.add_argument(
        '--chunk_size',
        required=True,
        type=int,
        default=1000,
        help='Chunk size to use when writing to file',
    )

    parser.add_argument(
        '--header_lines',
        required=True,
        type=int,
        help='Number of lines to skip when reading in the VEP file',
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
