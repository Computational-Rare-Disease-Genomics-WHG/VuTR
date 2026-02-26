"""Process Clinvar variants file to add transcript_id column to be database ready"""

import os
import sys
import argparse
import logging
import pandas as pd  # noqa # pylint: disable=import-error


def main(args):
    """
    Main entry point for the script.
    """
    # Check if the input file exists
    if not os.path.exists(args.input):
        logging.error('Input file does not exist')
        sys.exit(1)

    # Check if the mane file exists
    if not os.path.exists(args.mane):
        logging.error('Mane file does not exist')
        sys.exit(1)

    # Check if the output file exists
    if os.path.exists(args.output):
        logging.error('Output file already exists')
        sys.exit(1)

    # Read input file
    df = pd.read_csv(args.input, sep='\t')
    mane = pd.read_csv(args.mane, sep='\t')

    # Filter mane to 5' UTR
    mane = mane.loc[mane['type'] == 'five_prime_UTR']

    df['transcript_id'] = df.apply(
        lambda x: mane.loc[
            (mane['start'] <= x['pos']) &
            (mane['end'] >= x['pos']) &
            (mane['seqid'] == x['chrom']), 'transcript_id'
        ].tolist(), axis=1)
    df = df.explode('transcript_id')

    # Change the transcript id by renomving the ensemble version
    df['transcript_id'] = df['transcript_id'].str.split('.').str[0]

    # Add variant_id column and chr column
    df['chr'] = df['chrom'].str.replace('chr', '')

    # Delete chrom column
    df = df.drop(columns=['chrom'])

    # Write output file
    df.to_csv(args.output, sep='\t', index=False)

    # Log the end of the script
    logging.info('End of script')


if __name__ == '__main__':
    # Read input and output files
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, type=str, help='Input Clinvar file')
    parser.add_argument('--mane', required=True, type=str, help='Mane file')
    parser.add_argument('--output', required=True, type=str, help='Output file')
    main(parser.parse_args())
