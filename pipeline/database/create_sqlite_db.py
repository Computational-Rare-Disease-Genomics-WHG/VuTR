"""
Creates table names in the sqlite3 database
"""

# E. D'Souza
# Final part of the pipeline
# Takes the data generated by VEP from ClinVar and gnomAD
# Creates a sqlite3 database with the data

# from sqlalchemy import create_engine
# from sqlalchemy.schema import CreateTable


import sqlite3
import argparse
import os
import sys
from model import tbl_models  # pylint: disable=E0401


def main(args):
    """Entry point"""

    if os.path.isfile(args.db_name) and not args.overwrite:
        print(f'Database already exists at {args.db_name}')
        sys.exit(0)

    # Create connection and cursor with the database
    conn = sqlite3.connect(args.db_name)
    c = conn.cursor()

    if args.table_name == 'all':
        options = tbl_models.keys()
    else:
        options = [args.table_name]

    for opt in options:
        print(f'Creating table {opt} in database {args.db_name}')
        c.execute(tbl_models[opt])
        conn.commit()
    print(f'Completed creating tables')
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Creates a table that stores utr features'
    )
    parser.add_argument(
        '--table_name',
        required=True,
        type=str,
        default='all',
        help='Name of the table that we wish to create (Default: all)',
    )
    parser.add_argument(
        '--db_name',
        required=True,
        type=str,
        help='Output sqlite3 file name and location',
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
