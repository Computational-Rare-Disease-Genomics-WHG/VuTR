"""
Creates the visualisation intervals based on the variant data.

Usage:
python create_intervals.py \
    -i <input_file> \
    -t <transcript_sequences> \
    --orf_file <orf_file> \
    -o <output_file> \

The output file will be a tab separated file. Namely it will have a column for the
ach variant-consequence which has the following structure
saved as a json string in the intervals column of the output file.
{
        'variant_id' : '',
        'annotation_id' : '',
        'peturbing_orf_id' : '',
        'context':{
            'ref' :'',
            'alt' : '',
        },
        'kozak_context': {
            'ref' : '',
            'alt' : '',
        }, 
        'intervals':  {
            'transcript' : {
                'start' : '',
                'end' : '',
            },
        },
    }
"""

import argparse
import pandas
import tqdm
import os
import sys
import json

def find_orf_impacted(variant_position, orfs):
    """
    """
    pass

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


def ustop_lost(**kwargs):
    """
    Calculates the consequences of a uSTOP_lost variant
    """

def parse_five_prime_utr_variant_consequence(conseq_str):
    """
    Parses the consequence str into a keyed dictionary as per
    https://github.com/ImperialCardioGenetics/UTRannotator#the-detailed-annotation-for-each-consequence

    """
    return {
        annotation.split(':')[0]: annotation.split(':')[1]
        for annotation in conseq_str.split(',')
    }


def create_intervals(
        parsed_variants,
        ensembl_transcript_id,
        transcript_sequences, 
        exon_structure,
        strand,
        start_site_pos,
        intron_size=20, 
        buffer_size=40):
    """
    Creates a dataframe with the intervals for each variant
    @param parsed_variants: Dataframe with parsed variants
    @param ensembl_transcript_id: Ensembl transcript id
    @param transcript_sequences: Dataframe with transcript sequences
    @param mane_file: Dataframe with mane file
    @param intron_size: Size of intron to use
    @param buffer_size: Size of buffer to use
    @return: Dataframe with intervals
    """

    # Define the consequence functions
    consequence_functions = {
        'uSTOP_lost' : ustop_lost,
        'uSTOP_gained' : ustop_gained,
        'uAUG_lost' : uaug_lost,
        'uAUG_gained' : uaug_gained,
        'uFrameshift' : frameshift,
    }
    # Iterate over each variant consequence
    for variant in parsed_variants:
        # Get the variant position
        # Get the ORFS impacted

        # Parse the consequence string
        var_cons = parse_five_prime_utr_variant_consequence(
            variant['five_prime_UTR_variant_annotation']
        )

        visualised_cons = consequence_functions[variant['five_prime_UTR_variant_consequence']](
            ...,
       )
        
        # Get plotting coordinates
        visualised_cons[
            'intervals'
        ] = get_intervals(
            visualised_cons['transcript']['start'],
            visualised_cons['transcript']['end'],
            exon_structure,
            intron_size,
            buffer_size,
            strand
        )

        # Convert to json and serialise to string 
        variant['intervals']['visualisation'] = json.dumps(visualised_cons)
    return parsed_variants


def process_exons(
        mane_file,
        transcript_sequences,
        intron_size,
        buffer_size):
    """
    Process the exons to add the visualisation coordinates
    @param exons: Dataframe with exons
    @return: Dataframe with exons with visualisation coordinates
    """

    # Find the start site position
    exons = mane_file[mane_file['type'] == 'exon']
    exons = exons[['transcript_id', 'exon_number', 'start', 'end', 'strand']]

    # Merge with transcript sequences to get the start site position
    exons = exons.merge(
        transcript_sequences[['transcript_id', 'start_site_pos']],
        left_on='transcript_id',
        right_on='transcript_id',
        how='left'
    )

    # Strip the version number from the transcript_id
    exons['transcript_id'] = exons['transcript_id'].apply(
        lambda x: x.split('.')[0])
    
    # |Calculate the width of the exon
    exons['width'] = exons['end'] - exons['start'] + 1

    # Order the mane file by transcript_id and exon_number
    exons = exons.sort_values(
        ['transcript_id', 'exon_number']
    )
    # Create a cumulative sum of the exon widths
    exons['tstart'] = exons.groupby('transcript_id',
                                            group_keys=False)['width'].apply(
        lambda x: x.shift(fill_value=1).cumsum())
    exons['tend'] = exons.groupby('transcript_id')['width'].cumsum()

    # Add the visualisation coordinates
    exons['vis_start'] = exons['tstart'] + (
        exons['exon_number'] - 1) * intron_size
    exons['vis_end'] = exons['tend'] + (
        exons['exon_number'] - 1) * intron_size

    # TODO: When strand is '-', the visualisation coordinates are reversed
    # CONSIDER NEGATIVE STRANDS
    # exons.loc[exons['strand'] == '-', ['vis_start', 'vis_end']] = 

    # Drop the genomic start and end columns
    exons = exons.drop(columns=['start', 'end'])
    return exons


def main(args):
    """Main entry point"""

    # Check if output file exists
    output_file = args.output_file
    if os.path.isfile(output_file):
        print('File already exists. Aborting.')
        sys.exit(1)

    # Check if files exist
    if not os.path.isfile(args.mane_file):
        raise FileNotFoundError(f"File {args.mane_file} does not exist")
    if not os.path.isfile(args.transcript_sequences):
        raise FileNotFoundError(f"File {args.transcript_sequences} does not exist")
    if not os.path.isfile(args.orf_file):
        raise FileNotFoundError(f"File {args.orf_file} does not exist")
    if not os.path.isfile(args.parsed_variants):
        raise FileNotFoundError(f"File {args.parsed_variants} does not exist")

    # Reading in transcript sequences  and keep only the
    # transcript_id, start_site_pos, seq columns
    # and remove the version number from the transcript_id
    transcript_sequences = pandas.read_csv(args.transcript_sequences, sep="\t")
    transcript_sequences = transcript_sequences[['transcript_id',
                                                 'start_site_pos', 'seq']]
    transcript_sequences['transcript_id'] = transcript_sequences['transcript_id'].apply(
        lambda x: x.split('.')[0]
    )

    # Reading in mane file and filter to exon features
    # and keep only the transcript_id, exon_number, start, end, strand columns
    # and remove the version number from the transcript_id
    # calculate the length of the exon
    mane_file = pandas.read_csv(args.mane_file, sep="\t")
    mane_file = process_exons(mane_file, transcript_sequences, intron_size, buffer_size)

    # Reading in ORF file and filter to exon features
    orfs = pandas.read_csv(args.orf_file, sep="\t")
    orfs['ensembl_transcript_id'] = orfs['ensembl_transcript_id'].apply(
        lambda x: x.split('.')[0]
    )

    # Reading in parsed variants from VEP + UTR Annotator
    parsed_variants = pandas.read_csv(args.input_file, sep="\t")

    # Append new column to parsed_variants with intervals
    print(f'Writing to file {output_file}')
    df_groups = parsed_variants.groupby(
        'Feature'
    )

    # Write header for output file
    header = True
    for transcript_id, variants_per_transcript in tqdm.tqdm(df_groups):

        # Filter transcript sequences and mane file to transcript
        sequence = transcript_sequences[
            transcript_sequences['transcript_id'] == transcript_id]
        exon_structure = mane_file[mane_file['transcript_id'] == transcript_id]

        # Get the strand of the transcript
        strand = exon_structure['strand'].unique()[0]
        start_site_pos = exon_structure['start_site_pos'].unique()[0]

        # Create intervals
        output = create_intervals(
            parsed_variants=variants_per_transcript,
            ensembl_transcript_id=transcript_id,
            transcript_sequences=sequence,
            exon_structure=exon_structure,
            strand=strand,
            start_site_pos=start_site_pos,
            intron_size=intron_size,
            buffer_size=buffer_size
        )

        # Write to file
        output.to_csv(
            output_file,
            sep='\t',
            mode='a',
            header=header,
            index=False
        )

        header = False


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates the visualisation intervals based on the variant data"  # noqa: E501 # pylint: disable=C0301
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
        "--orf_file",
        required=True,
        type=str,
        help="ORF file",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        type=str,
        help="Output file name and location",
    )
    main(args=parser.parse_args())
