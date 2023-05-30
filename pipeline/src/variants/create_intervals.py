"""
Creates the visualisation intervals based on the variant data.

Usage: 
python3 create_intervals.py 
    -i <input_file> 
    -m <mane_file>
    -o <output_file>

Each variant-consequence has the following structure
saved as a json string in the intervals column of the output file.  
{
        'variant_id' : '',
        'annotation_id' : '',
        'peturbing_orf_id' : '',
        'kozak_strength':{
            'ref' :'',
            'alt' : '',
        },

        'intervals':  {
            'transcript' : {
                'start' : '',
                'end' : '',
            },
            'visualisation': [
                {
                    'x' : '', // Visualisation coordinates 
                    'y' : '', 
                    ...
                }
            ] 
        },
    }
"""

import argparse
import pandas
import tqdm
import os
import sys
import json




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



def get_genomic_vis_intervals(
        transcript_start,
        transcript_end,
        exons,
        intron_size,
        buffer_size,
        strand
    ):
    """
    Converts the transcript coordinates a list of genomic intervals
    @param transcript_start: Start position of transcript
    @param transcript_end: End position of transcript
    @param exons: List of exons
    @return: List of dictionaries with genomic intervals keyed by x and y
    """
    if strand == '+':
        # Find the exon_number that the transcript_start and transcript_end fall between
        exon_number_start = exons[exons['start'] <= transcript_start]['exon_number'].max()
        exon_start = exons[exons['exon_number'] == exon_number_start]['start'].values[0]
        delta_start = transcript_start - exon_start

        # Get the exon start and end positions
        exon_number_end = exons[exons['start'] <= transcript_end]['exon_number'].max()
        exon_end = exons[exons['exon_number'] == exon_number_end]['end'].values[0]
        delta_end = transcript_end - exon_end
        
        if exon_number_start == exon_number_end:
            exon_coordinate = exons[exons['exon_number'] == exon_number_start]['vis_start'].values[0] 
            return [
                {
                    'x' : exon_coordinate + delta_start,
                    'y' : exon_coordinate + delta_end, 
                }
            ]
        else:
            # TODO: List the values between the start and end
            # Iterate over the exons from start to end.
            pass

    else:
        # TODO: Do the case for the reverse strand
        pass

    return []



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
        strand,
        mane_file, 
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
        # Get the strand
        # Get the exon
        # Get the transcript sequence
        # Get the variant position
        # Get the ORFS impacted

        # Parse the consequence string
        var_cons = parse_five_prime_utr_variant_consequence(
            variant['five_prime_UTR_variant_annotation']
        )

        visualised_cons = consequence_functions[variant['five_prime_UTR_variant_consequence']](
        )

        # Get plotting coordinates
        visualised_cons[
            'intervals'
        ] = get_genomic_vis_intervals(
            visualised_cons['transcript']['start'],
            visualised_cons['transcript']['end'],
            mane_file,
            intron_size,
            buffer_size,
            strand
        )

        # Convert to json and serialise to string 
        variant['intervals']['visualisation'] = json.dumps(visualised_cons)
    return parsed_variants


def main(args):
    """Main entry point"""

    # Check if output file exists
    output_file = args.output_file
    if os.path.isfile(output_file):
        print('File already exists. Aborting.')
        sys.exit(1)

    # Set intron and buffer size
    intron_size = args.intron_size
    buffer_size = args.buffer_size

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
    mane_file = mane_file[mane_file['type'] == 'exon']
    mane_file = mane_file[['transcript_id', 'exon_number', 'start', 'end', 'strand']]
    mane_file['transcript_id'] = mane_file['transcript_id'].apply(
        lambda x: x.split('.')[0])
    mane_file['width'] = mane_file['end'] - mane_file['start'] + 1

    # Order the mane file by transcript_id and exon_number
    mane_file = mane_file.sort_values(
        ['transcript_id', 'exon_number']
    )
    
    # TODO: Create a cumulative sum of the exon widths
    # TODO: Add intron size to the cumulative sum
    # TODO: Drop the genomic start and end columns


    # Reading in ORF file and filter to exon features
    orfs = pandas.read_csv(args.orf_file, sep="\t")
    orfs['ensembl_transcript_id'] = orfs['ensembl_transcript_id'].apply(
        lambda x: x.split('.')[0])

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

        # Create intervals
        output = create_intervals(
            parsed_variants=variants_per_transcript,
            ensembl_transcript_id=transcript_id,
            transcript_sequences=sequence,
            mane_file=exon_structure,
            strand=strand,
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
