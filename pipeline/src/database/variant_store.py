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
from concurrent.futures import ThreadPoolExecutor


CLINVAR_TABLE_SQL_QUERY = """
    CREATE TABLE IF NOT EXISTS clinvar_variants (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        variant_id VARCHAR,
        chr VARCHAR,
        pos INTEGER,
        ref VARCHAR,
        alt VARCHAR,
        clinvar_id VARCHAR,
        alleleid VARCHAR,
        clndn VARCHAR,
        clndnincl VARCHAR,
        clndisdb VARCHAR,
        clndisdbincl VARCHAR,
        clnhgvs VARCHAR,
        clnrevstat VARCHAR,
        clnsig VARCHAR,
        clnsigconf VARCHAR,
        clnsigincl VARCHAR,
        clnvc VARCHAR,
        clnvco VARCHAR,
        clnvi VARCHAR,
        dbvarid VARCHAR,
        geneinfo VARCHAR,
        mc VARCHAR,
        origin VARCHAR,
        rs VARCHAR,
        ensembl_transcript_id VARCHAR
);
"""

GNOMAD_TABLE_SQL_QUERY = """
    CREATE TABLE IF NOT EXISTS gnomad_variants (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        variant_id VARCHAR,
        chr VARCHAR,
        pos INTEGER,
        ref VARCHAR,
        alt VARCHAR,
        rsid VARCHAR,
        ac INTEGER,
        an INTEGER,
        af REAL,
        variant_type VARCHAR,
        ensembl_transcript_id VARCHAR
    );
"""

VARIANT_ANNOTATION_SQL_QUERY = """
    CREATE TABLE IF NOT EXISTS variant_annotations (
        ensembl_transcript_id varchar,
        variant_id varchar,
        cdna_pos int,
        dataset varchar,
        five_prime_UTR_variant_consequence varchar,
        five_prime_UTR_variant_annotation data,
        annotations data
    );
"""

def parse_five_prime_utr_variant_consequence(conseq_str):
    """
    Parses the five_prime_UTR_variant_consequence column into a dictionary
    """
    return {annotation.split(':')[0]: annotation.split(':')[1] for annotation in conseq_str.split(',')}

def convert_uploaded_variation_to_variant_id(uploaded_variation):
    """
    Converts the uploaded_variation column into a variant_id
    """
    return uploaded_variation.replace('_', '-').replace('/', '-')


def variant_annotations_db_writer(db_conn, rows):
    """
    Writes the rows to the database for simulated
    """
    query = """
        INSERT INTO variant_annotations (
            ensembl_transcript_id, variant_id, cdna_pos,
            five_prime_UTR_variant_consequence,
            five_prime_UTR_variant_annotation, annotations, dataset
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
    """
    c = db_conn.cursor()
    c.executemany(query, rows)

def process_batch(batch_df, dataset):
    """
    Processes a batch of variants and their consequences
    """
    rows = []
    convert_variant_id = convert_uploaded_variation_to_variant_id
    parse_consequence = parse_five_prime_utr_variant_consequence

    for _, row in batch_df.iterrows():
        variant_conseq = row.to_dict()
        feature = variant_conseq['Feature']
        uploaded_variation = convert_variant_id(variant_conseq['#Uploaded_variation'])
        cDNA_position = variant_conseq['cDNA_position']
        five_prime_UTR_variant_consequence = variant_conseq['five_prime_UTR_variant_consequence']
        utr_variant_annotation = parse_consequence(variant_conseq['five_prime_UTR_variant_annotation'])
        variant_conseq_json = json.dumps(variant_conseq)
        rows.append((feature, uploaded_variation, cDNA_position, five_prime_UTR_variant_consequence, json.dumps(utr_variant_annotation), variant_conseq_json, dataset))
    return rows
    

def process_gnomad_batch(batch_df, c):
    """
    Processes a batch of gnomAD variants
    """
    rows = [tuple(row) for row in batch_df[[
        'variant_id', 'chr','pos', 'ref', 'alt', 'rsid', 'ac', 'an', 'af', 'variant_type',
        'transcript_id'
    ]].to_records(index=False)]

    # Using parameterized SQL to avoid SQL injection
    c.executemany(
        '''
        INSERT INTO gnomad_variants (
            variant_id, chr, pos,
            ref, alt, rsid, ac, an, af,
            variant_type, ensembl_transcript_id
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''',
        rows
    )


def process_clinvar_batch(batch_df, c):
    """
    Processes a batch of clinvar variants
    """
    rows = [tuple(row) for row in batch_df[[
        'variant_id', 'chr', 'pos', 'ref', 'alt', 'clinvar_id', 'alleleid', 'clndn',
        'clndnincl', 'clndisdb', 'clndisdbincl', 'clnhgvs', 'clnrevstat', 'clnsig',
        'clnsigconf', 'clnsigincl', 'clnvc', 'clnvco', 'clnvi', 'dbvarid',
        'geneinfo', 'mc', 'origin', 'rs', 'transcript_id'
    ]].to_records(index=False)]

    # Using parameterized SQL to avoid SQL injection
    c.executemany(
        '''
        INSERT INTO clinvar_variants (
            variant_id, chr, pos, ref,
            alt, clinvar_id, alleleid, clndn,
            clndnincl, clndisdb, clndisdbincl, clnhgvs,
            clnrevstat, clnsig, clnsigconf, clnsigincl,
            clnvc, clnvco, clnvi, dbvarid,
            geneinfo, mc, origin, rs, ensembl_transcript_id
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''',
        rows
    )

def main(args):
    """
    Main entry point
    """
    if os.path.isfile(args.db_name) and not args.overwrite:
        print(f'Database already exists at {args.db_name}')
        sys.exit(0)

    with sqlite3.connect(args.db_name) as conn:
        conn.execute('PRAGMA journal_mode = WAL')  # Enable WAL mode
        conn.execute('PRAGMA synchronous = OFF')  # Disable synchronous mode
        conn.execute('PRAGMA locking_mode=EXCLUSIVE')
        conn.execute('PRAGMA mmap_size=30000000000')
        c = conn.cursor()

        # Add clinvar variants
        if args.clinvar:
            print(f'Adding clinvar variants schema to database {args.db_name}')
            c.execute(CLINVAR_TABLE_SQL_QUERY)
            process_clinvar_batch(pd.read_csv(args.clinvar, sep='\t'), c)


        # Add gnomAD variants
        if args.gnomad:
            print('Adding gnomAD variants')
            c.execute(GNOMAD_TABLE_SQL_QUERY)
            batch_size = 10000
            chunk_iter = pd.read_csv(args.gnomad, sep='\t', chunksize=batch_size)
            for chunk_df in tqdm.tqdm(chunk_iter):
                process_gnomad_batch(chunk_df, c)


        print(f'Creating table variant table in database {args.db_name}')
        c.execute(VARIANT_ANNOTATION_SQL_QUERY)

        # Add high impact variants
        print('Adding high impact variants')
        clinvar_rows = process_batch(pd.read_csv(args.clinvar_himpact, sep='\t'), 'clinvar')
        gnomad_rows = process_batch(pd.read_csv(args.gnomad_himpact, sep='\t'), 'gnomad')
        variant_annotations_db_writer(conn, clinvar_rows+gnomad_rows)

        # Add simulated variants
        batch_size = 5000  # Adjust the batch size based on available memory
        chunk_iter = pd.read_csv(args.variant_file, sep='\t', chunksize=batch_size)
        
        for chunk in tqdm.tqdm(chunk_iter):
            rows = process_batch(chunk, 'simulated')
            variant_annotations_db_writer(conn, rows)
        
        # Commit the transaction after all processes have completed
        print(f'Completed adding variants to database {args.db_name}')

        # Create indexes
        print('Creating indexes')
        c.execute('CREATE INDEX idx_variant_annotations ON variant_annotations(ensembl_transcript_id);')
        c.execute('CREATE INDEX idx_gnomad_variants ON gnomad_variants(ensembl_transcript_id);')
        c.execute('CREATE INDEX idx_clinvar_variants ON clinvar_variants(ensembl_transcript_id);')

        print('Completed creating indexes')
        conn.commit()

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
    parser.add_argument('--gnomad', type=str, help='gnomAD UTR variants file')
    parser.add_argument('--clinvar', type=str, help='Clinvar UTR variants tsv file')
    parser.add_argument('--gnomad-himpact', type=str, help='gnomAD high impact variants file')
    parser.add_argument('--clinvar-himpact', type=str, help='Clinvar high impact variants file')

    main(args=parser.parse_args())
