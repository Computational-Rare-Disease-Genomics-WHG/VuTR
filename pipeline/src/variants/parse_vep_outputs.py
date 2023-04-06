"""

Parses the outputs from VEP running UTR annotator and saves it as a .tsv
Potential performance improvement : Use Dask as the dataframe are really big
"""
import argparse
import pandas as pd
import numpy as np
import tqdm


def wide_to_long(df):
    """
    Make from wide to long.
    @param df with 5'utr consequence concatentated by &.
    @returns long_df What we want is a
            row per variant / transcript / five_prime_UTR_variant_annotation.
    """
    # split the columns into lists of values
    df['five_prime_UTR_variant_consequence'] = df['five_prime_UTR_variant_consequence'].str.split('&')
    df['five_prime_UTR_variant_annotation'] = df['five_prime_UTR_variant_annotation'].str.split('&')

    # use pd.concat and pd.melt to transform the dataframe
    dfs = [df.loc[~df['five_prime_UTR_variant_consequence'].str.contains('&')], 
           pd.concat([df.loc[df['five_prime_UTR_variant_consequence'].str.contains('&')].reset_index(drop=True),
                      pd.melt(df.loc[df['five_prime_UTR_variant_consequence'].str.contains('&')], 
                              id_vars=df.columns[:-2], 
                              value_vars=['five_prime_UTR_variant_consequence', 'five_prime_UTR_variant_annotation'],
                              var_name='variable', 
                              value_name='value')], axis=1)]

    # concatenate the dataframes
    long_df = pd.concat(dfs, ignore_index=True)

    # set the index and drop unnecessary columns
    long_df = long_df.set_index([c for c in df.columns[:-2] if c != 'five_prime_UTR_variant_consequence']).drop(['variable', 'index'], axis=1)

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

    print(f'Writing to file {write_path}')
    # Convert wide to long data frame
    for _group_idx, df_group in tqdm(df_groups):
        wide_to_long(df_group).to_csv(write_path, sep='\t', mode='a',
                                      header=False, index=False)

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
