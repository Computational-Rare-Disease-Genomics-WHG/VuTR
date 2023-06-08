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
    return {
        annotation.split(':')[0]: annotation.split(':')[1]
        for annotation in conseq_str.split(',')
    }

def convert_uploaded_variation_to_variant_id(uploaded_variation):
    return uploaded_variation.replace('_', '-').replace('/', '-')

def process_batch(batch_df, c):
    rows = []
    for _, row in batch_df.iterrows():
        variant_conseq = row.to_dict()
        rows.append((
            variant_conseq['Feature'],
            convert_uploaded_variation_to_variant_id(variant_conseq['#Uploaded_variation']),
            variant_conseq['cDNA_position'],
            variant_conseq['five_prime_UTR_variant_consequence'],
            json.dumps(parse_five_prime_utr_variant_consequence(variant_conseq['five_prime_UTR_variant_annotation'])),
            json.dumps(variant_conseq),
        ))
    c.executemany('INSERT INTO variant_annotations VALUES (?, ?, ?, ?, ?, ?)', rows)

def main(args):
    if os.path.isfile(args.db_name) and not args.overwrite:
        print(f'Database already exists at {args.db_name}')
        sys.exit(0)

    conn = sqlite3.connect(args.db_name)
    conn.execute('PRAGMA journal_mode = WAL')  # Enable WAL mode
    c = conn.cursor()

    print(f'Creating table variant table in database {args.db_name}')
    c.execute("""
        CREATE TABLE IF NOT EXISTS variant_annotations (
            ensembl_transcript_id varchar,
            variant_id varchar,
            cdna_pos int,
            five_prime_UTR_variant_consequence varchar,
            five_prime_UTR_variant_annotation data,
            annotations data
        )"""
    )
    print(f'Completed creating tables')

    batch_size = 10000  # Adjust the batch size based on available memory
    chunk_iter = pd.read_csv(args.variant_file, sep='\t', chunksize=batch_size)
    for chunk_df in tqdm.tqdm(chunk_iter):
        process_batch(chunk_df, c)
        conn.commit()

    conn.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Ingresses all variant data into a sqlite3 database')
    parser.add_argument('--db_name', required=True, type=str, help='Output sqlite3 file name and location')
    parser.add_argument('--variant_file', required=True, type=str, help='The variant tsv file to be ingested')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing database if exists')
    parser.add_argument('--verbose', action='store_true', help='Verbose outputs')
    main(args=parser.parse_args())
