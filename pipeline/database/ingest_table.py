"""Takes a table and ingests it into the sqlite3 database"""


import sqlite3
import argparse
import sys
import os
import pandas as pd
from .model import tbl_models


def main(args):
    """Entry point"""
    # check db name
    if os.path.isfile(args.db_name) and not args.overwrite:
        print(f'Database already exists at {args.db_name}')
        sys.exit(0)

    # Create connection and cursor with the database
    conn = sqlite3.connect(args.db_name)
    # c = conn.cursor()

    # List of tables to create
    tbl_list = []
    if args.single:
        tbl_list = [args.tbl_name]
    elif args.all:
        tbl_list = tbl_models.keys()

    # read the database
    for tbl in tbl_list:

        # Read the sql table
        df = pd.read_csv('')

        # Write to SQLtables
        df.to_sql(
            tbl,
            conn,
            if_exists='fail',
            index=False,
            chunksize=1000,
            dtype=tbl_models[tbl],
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Ingresses all variant data into a sqlite3 database'
    )
    parser.add_argument(
        '--db_name',
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
