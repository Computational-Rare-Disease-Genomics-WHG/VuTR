"""
Stores variants (and their consequences) in a key : value store database
"""


import sqlite3
import argparse
import os
import sys
import json
import tqdm

import pandas as pd


def parse_five_prime_utr_variant_consequence(conseq_str):
    """
    Parses the consequence str into a keyed dictionary as per
    https://github.com/ImperialCardioGenetics/UTRannotator#the-detailed-annotation-for-each-consequence

    """
    return {
        annotation.split(':')[0]: annotation.split(':')[1]
        for annotation in conseq_str.split(',')
    }


def convert_uploaded_variation_to_variant_id(uploaded_variation):
    """
    Replaces the uploaded variation in VEP to a gnomad-esq variant id
    """
    return uploaded_variation.replace('_', '-').replace('/', '-')


def main(args):
    """Entry point"""

    if os.path.isfile(args.db_name) and not args.overwrite:
        print(f'Database already exists at {args.db_name}')
        sys.exit(0)

    # Create connection and cursor with the database
    conn = sqlite3.connect(args.db_name)
    c = conn.cursor()

    print(f'Creating table variant table in database {args.db_name}')
    c.execute(
        """CREATE TABLE IF NOT EXISTS variant_annotations (
            ensembl_transcript_id varchar,
            variant_id varchar,
            cdna_pos int,
            five_prime_UTR_variant_consequence varchar,
            five_prime_UTR_variant_annotation data,
            annotations data)"""
    )
    conn.commit()
    print(f'Completed creating tables')

    variant_df = pd.read_csv(args.variant_file, sep='\t')
    variant_df = variant_df.drop_duplicates()
    for index, row in tqdm.tqdm(variant_df.iterrows()):
        variant_conseq = row.to_dict()
        c.execute(
            'insert into variant_annotations values (?, ?, ?,?, ?, ?)',
            [
                variant_conseq['Feature'],
                convert_uploaded_variation_to_variant_id(
                    variant_conseq['#Uploaded_variation']
                ),
                variant_conseq['cDNA_position'],
                variant_conseq['five_prime_UTR_variant_consequence'],
                json.dumps(
                    parse_five_prime_utr_variant_consequence(
                        variant_conseq['five_prime_UTR_variant_annotation']
                    )
                ),
                json.dumps(variant_conseq),
            ],
        )
        conn.commit()

    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Ingresses all variant data into a sqlite3 database'
    )
    parser.add_argument(
        '--db_name',
        required=True,
        type=str,
        help='Output sqlite3 file name and location',
    )
    parser.add_argument(
        '--variant_file',
        required=True,
        type=str,
        help='The variant tsv file to be injested',
    )
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Overwrite existing database if exists',
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose outputs',
    )
    main(args=parser.parse_args())
