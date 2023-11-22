"""
Stores variants (and their consequences) in a key : value store database
"""

import sqlite3
import argparse
import os
import sys
import json
import tqdm  # pylint: disable=import-error
import pandas as pd  # pylint: disable=import-error
import multiprocessing




def parse_five_prime_utr_variant_consequence(conseq_str):
    """
    Parses the five_prime_UTR_variant_consequence column into a dictionary
    """
    return {
        annotation.split(':')[0]: annotation.split(':')[1]
        for annotation in conseq_str.split(',')
    }

def convert_uploaded_variation_to_variant_id(uploaded_variation):
    """
    Converts the uploaded_variation column into a variant_id
    """
    return uploaded_variation.replace('_', '-').replace('/', '-')


def process_batch(conn, batch_df):
    """
    Processes a batch of variants and their consequences
    """
    c = conn.cursor()
    rows = []
    for _, row in batch_df.iterrows():
        variant_conseq = row.to_dict()
        rows.append((
            variant_conseq['Feature'],
            convert_uploaded_variation_to_variant_id(
                variant_conseq['#Uploaded_variation']
            ),
            variant_conseq['cDNA_position'],
            variant_conseq['five_prime_UTR_variant_consequence'],
            json.dumps(parse_five_prime_utr_variant_consequence(
                variant_conseq['five_prime_UTR_variant_annotation'])),
            json.dumps(variant_conseq),
        ))
    c.executemany('INSERT INTO variant_annotations VALUES (?, ?, ?, ?, ?, ?)', rows)
    conn.commit()


def process_gnomad_batch(batch_df, c):
    """
    Processes a batch of gnomAD variants
    """
    rows = []
    for _, row in batch_df.iterrows():
        variant = row.to_dict()
        chrom = variant['chr'].replace('chr', '')
        variant_id = f'{chrom}-{variant["pos"]}-{variant["ref"]}-{variant["alt"]}'
        rows.append((
            variant_id,
            variant['ensembl_transcript_id'],
            chrom,
            variant['pos'],
            variant['ref'],
            variant['alt'],
            variant['af'],
            variant['ac'],
            variant['an']
        ))
    c.executemany('INSERT INTO gnomad_variants VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', rows)


def process_clinvar_batch(batch_df, c):
    """
    Processes a batch of clinvar variants
    """
    rows = []
    for _, row in batch_df.iterrows():
        variant = row.to_dict()
        rows.append((
            variant['variant_id'],
            variant['ensembl_transcript_id'],
            variant['chrom'],
            variant['pos'],
            variant['ref'],
            variant['alt'],
            variant['clinvar_id'],
            variant['clnsig'],
            variant['clinvar_review'],
            variant['clnrevstat'],
            variant['clinvar_star'],
        ))
    c.executemany('INSERT INTO clinvar_variants VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', rows)


def main(args):
    """
    Main entry point
    """
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
        )""")
    print('Completed creating tables')

    # Set up process pool
    num_processes = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)

    batch_size = 80000  # Adjust the batch size based on available memory
    conn.execute('BEGIN TRANSACTION')
    chunk_iter = pd.read_csv(args.variant_file, sep='\t', chunksize=batch_size)
    results = []

    for chunk_df in tqdm.tqdm(chunk_iter):
        results.append(
            pool.apply_async(process_batch, (conn, chunk_df))
        )
    conn.commit()

    # Add gnomAD variants
    if args.gnomad:
        print('Adding gnomAD variants')
        c.execute("""
            CREATE TABLE IF NOT EXISTS gnomad_variants (
                variant_id varchar PRIMARY KEY,
                ensembl_transcript_id varchar,
                chr varchar,
                pos int,
                ref varchar,
                alt varchar,
                af float,
                ac int,
                an int
            )""")
        batch_size = 10000
        chunk_iter = pd.read_csv(args.gnomad, sep='\t', chunksize=batch_size)
        for chunk_df in tqdm.tqdm(chunk_iter):
            process_gnomad_batch(chunk_df, c)
            conn.commit()
    
    # Add clinvar variants
    if args.clinvar:
        print(f'Adding clinvar variants')
        c.execute("""
            CREATE TABLE IF NOT EXISTS clinvar_variants (
                variant_id varchar PRIMARY KEY,
                ensembl_transcript_id varchar,
                chr varchar,
                pos int,
                ref varchar,
                alt varchar,
                clinvar_id varchar,
                clinvar_clnsig varchar,
                clinvar_review varchar,
                clinvar_clnrevstat varchar,
                clinvar_star varchar);
        """)
        batch_size = 10000
        chunk_iter = pd.read_csv(args.clinvar, sep='\t', chunksize=batch_size)

        for chunk_df in tqdm.tqdm(chunk_iter):
            process_clinvar_batch(chunk_df, c)
            conn.commit()

    print(f'Completed adding variants to database {args.db_name}')

    # Create indexes
    print('Creating indexes')
    c.execute("""
CREATE INDEX idx_variant_annotations ON variant_annotations(ensembl_transcript_id);
CREATE INDEX idx_gnomad_variants ON gnomad_variants(ensembl_transcript_id);
CREATE INDEX idx_clinvar_variants ON clinvar_variants(ensembl_transcript_id);
    """)
    print('Completed creating indexes')

    conn.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Ingresses all variant data into a sqlite3 database')
    parser.add_argument(
        '--db_name', required=True, type=str,
        help='Output sqlite3 file name and location'
        )
    parser.add_argument(
        '--variant_file', required=True, type=str, 
        help='The variant tsv file to be ingested')
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Overwrite existing database if exists')
    parser.add_argument('--verbose', action='store_true', help='Verbose outputs')
    parser.add_argument('--gnomad', type=str, help='gnomAD variants tsv file')
    parser.add_argument('--clinvar', type=str, help='Clinvar variants tsv file')
    main(args=parser.parse_args())
