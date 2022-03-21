"""Parses the outputs from VEP running UTR annotator and saves it as a .tsv

Potential performance improvement : Use Dask as the dataframe are really big
"""
import argparse
import pandas as pd
import numpy as np


def wide_to_long(df):
    """
    Make from wide to long.
    @param df with 5'utr consequence concatentated by &.
    @returns long_df What we want is a
            row per variant / transcript / five_prime_UTR_variant_annotation.
    """
    long_df = pd.DataFrame(columns=df.columns.values.tolist())

    for _index, row in df.iterrows():

        binding_df = pd.DataFrame()
        # If single just row bind it to long_df
        if '&' not in row['five_prime_UTR_variant_consequence']:
            binding_df = pd.DataFrame(row.to_frame().T)
        else:
            # Replicate the rows per ORF consequences
            nrows = row['five_prime_UTR_variant_consequence'].count('&')
            consequences_split = row['five_prime_UTR_variant_consequence'].split('&')
            annotation_split = row['five_prime_UTR_variant_annotation'].split('&')
            binding_df = pd.concat(
                [
                    row.to_frame().T.drop(
                        columns=[
                            'five_prime_UTR_variant_consequence',
                            'five_prime_UTR_variant_annotation',
                        ]
                    )
                ]
                * (nrows + 1),
                ignore_index=True,
            )
            binding_df['five_prime_UTR_variant_consequence'] = consequences_split
            binding_df['five_prime_UTR_variant_annotation'] = annotation_split
        long_df = pd.concat([long_df, binding_df], axis=0, ignore_index=True)

    return long_df


def main(args):
    """
    Main entry point
    """
    vep_file_path = args.vep_file
    mane_version = args.mane_version
    mane_summary_path = f'../../data/pipeline/MANE/{mane_version}/MANE.GRCh38.v{mane_version}.summary.txt.gz'  # noqa: E501 # pylint: disable=C0301

    write_path = args.output_file

    print(f'Reading file {vep_file_path}')
    # Read clinvar file and the mane summary file
    vep_df = pd.read_csv(vep_file_path, sep='\t', skiprows=args.header_lines)

    print(f'Read completed')

    mane_summary_df = pd.read_csv(mane_summary_path, sep='\t')

    print(f'Filtering out 5 prime UTR variants with no consequence')
    # Filter to variants that impact uORFs (remove all other variants)
    no_utr_consequence = ['-', '', None, np.nan]
    vep_df = vep_df[
        ~vep_df['five_prime_UTR_variant_annotation'].isin(no_utr_consequence)
    ]

    # Remove version number
    mane_summary_df['transcript_id'] = mane_summary_df['Ensembl_nuc'].apply(
        lambda x: x[0:15]
    )

    print(f'Filtering to MANE consequences')
    # Filter to consequence on the MANE transcript
    vep_df = vep_df[vep_df['Feature'].isin(mane_summary_df['transcript_id'])]

    print(f'Transforming table to long')
    # Convert wide to long data frame
    long_vep_df = wide_to_long(vep_df)

    print(f'Writing to file')
    # Write to file
    long_vep_df.to_csv(write_path, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parses the vep output files')
    parser.add_argument(
        '--mane_version',
        default='1.0',
        help='Which mane version to use?',
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
