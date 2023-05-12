"""
Creates the visualisation intervals based on the variant data.

Usage: 
python3 create_intervals.py 
    -i <input_file> 
    -m <mane_file>
    -o <output_file>
"""

import argparse
import pandas
import tqdm
import os
import sys
import json

consequences = {
        'variant_id' : '',
        'annotation_id' : '',

        'transcript_start':  '',
        'transcript_end': '',

        '' : '',

        'intervals' : [],
    }


def ustop_lost():
    """
    """
    pass
    return consequences


def ustop_gained():
    """
    """
    pass


def uaug_lost():
    """
    
    """
    pass

def frameshift():
    """
    
    """
    pass

def uaug_gained():
    """
    
    """
    
    pass

def parse_five_prime_utr_variant_consequence(conseq_str):
    """
    Parses the consequence str into a keyed dictionary as per
    https://github.com/ImperialCardioGenetics/UTRannotator#the-detailed-annotation-for-each-consequence

    """
    return {
        annotation.split(':')[0]: annotation.split(':')[1]
        for annotation in conseq_str.split(',')
    }

def create_intervals(parsed_variants, transcript_sequences, mane_file, intron_size=20, buffer_size=40):
    """
    
    """
    
    # Filter sequences to 
    features = parsed_variants['Feature'].unique()
    transcript_sequences = transcript_sequences[transcript_sequences['Feature'].isin(features)]

    # Find the UTR exons
    exons = mane_file[mane_file['Feature_type'] == 'UTR_exon' && mane_file['Feature'].isin(features)]

    # Iterate over each variant
    consequence_functions = {
        'uSTOP_lost' : ustop_lost,
        'uSTOP_gained' : ustop_gained,
        'uAUG_lost' : uaug_lost,
        'uAUG_gained' : uaug_gained,
        'uFrameshift' : frameshift,
    }

    for variant in parsed_variants:

        # Get the strand 
        # Get the exon
        # Get the transcript sequence
        # Get the variant position

        # Parse the consequence string
        var_cons = parse_five_prime_utr_variant_consequence(
            variant['five_prime_UTR_variant_annotation']
        )

        visualised_cons = consequence_functions[variant['five_prime_UTR_variant_consequence']](

        )
        # Convert to json and serialise to string 
        variant['intervals'] = json.dumps(visualised_cons)
    return parsed_variants

def main(args):
    """Main entry point"""

    # Check if output file exists
    output_file = args.output_file
    if os.path.isfile(output_file):
        print('File already exists. Aborting.')
        sys.exit(1)

    intron_size = args.intron_size
    buffer_size = args.buffer_size

    # Check if files exist
    if not os.path.isfile(args.mane_file):
        raise FileNotFoundError(f"File {args.mane_file} does not exist")
    if not os.path.isfile(args.transcript_sequences):
        raise FileNotFoundError(f"File {args.transcript_sequences} does not exist")
    if not os.path.isfile(args.parsed_variants):
        raise FileNotFoundError(f"File {args.parsed_variants} does not exist")

    # Reading in transcript sequences  and keep transcript_id, start_site_pos, seq columns, 
    transcript_sequences = pandas.read_csv(args.transcript_sequences, sep="\t")
    transcript_sequences = transcript_sequences[['transcript_id', 'start_site_pos', 'seq']]

    # Reading in mane file
    mane_file = pandas.read_csv(args.mane_file, sep="\t")

    # Filter to exon features
    mane_file = mane_file[mane_file['type'] == 'exon']
    mane_file = mane_file[['transcript_id', 'exon_number', 'start', 'end', 'strand']]

    parsed_variants = pandas.read_csv(args.input_file, sep="\t")


    # Append new column to parsed_variants with intervals


    print(f'Writing to file {output_file}')

    df_groups = parsed_variants.groupby(
        'Feature'
    )
    header = True
    for _group_idx, df_group in tqdm.tqdm(df_groups):
        create_intervals(df_group).to_csv(output_file, sep='\t', mode='a',
                                      header=header, index=False)
        header = False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a VCF file that creates all possible UTR variants"  # noqa: E501 # pylint: disable=C0301
    )

    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=str,
        help="Input file from VEP + UTR Annotator",
    )
    parser.add_argument(
        "-m",
        "--mane_file",
        required=True,
        type=str,
        help="Which mane version to use?",
    )
    parser.add_argument(
        "-t",
        "--transcript_sequences",
        required=True,
        type=str,
        help="Transcript sequences file",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        type=str,
        help="Output file name and location",
    )
    parser.add_argument(
        "-is",
        "--intron_size",
        default=20,
        type=int,
        help="Size of intron to use when creating intervals",
    )
    parser.add_argument(
        "-bs",
        "--buffer_size",
        default=40,
        type=int,
        help="Size of buffer to use when creating intervals",
    )
    main(args=parser.parse_args())
