"""
Takes a table and ingests it into the sqlite3 database

Usage : python3 ingest_table.py --db_name db.sqlite --all --verbose --overwrite

TODO : Dtype specification issue doesn't work
    either through sqlcol or by manual specification in model.py
"""
from pathlib import Path
import sqlite3
import argparse
import sys
import os
import sqlalchemy

from model import tbl_models  # pylint: disable=E0401
import pandas as pd

SCRIPT_PATH = Path(__file__).parent


def sqlcol(dfparam):
    """
    Create a dictionary for native SQL alchemy types
    from a dataframe column for use in dtype paramater
    in pandas DataFrame.to_sql

    Note : Doesn't work for string for some reason
    hence not used in main()
    """

    dtypedict = {}
    for i, j in zip(dfparam.columns, dfparam.dtypes):
        if 'object' in str(j):
            dtypedict.update({i: sqlalchemy.types.NVARCHAR()})

        if 'datetime' in str(j):
            dtypedict.update({i: sqlalchemy.types.DateTime()})

        if 'float' in str(j):
            dtypedict.update({i: sqlalchemy.types.Float(precision=3, asdecimal=True)})

        if 'int' in str(j):
            dtypedict.update({i: sqlalchemy.types.INT()})

    return dtypedict


def strip_version_identifiers(df, id_names):
    """
    Use grep to remove trailing ENSTXXX[.3] - >  ENSTXXX
    @param df
    """

    df.loc[:, id_names] = df.loc[:, id_names].replace(
        to_replace='\.\d+$', value='', regex=True  # pylint: disable=W1401 # noqa: W605
    )
    return df


def main(args):
    """Entry point"""
    # check db name
    if os.path.isfile(args.db_name) and not args.overwrite:
        print(f'Database already exists at {args.db_name}')
        sys.exit(0)

    # Create connection and cursor with the database
    conn = sqlite3.connect(args.db_name)

    # List of tables to create based on query
    tbl_list = []
    if args.single:
        tbl_list = [args.tbl_name]
    elif args.all:
        tbl_list = tbl_models.keys()

    for tbl in tbl_list:

        # Read the sql table
        print(f'Reading table {tbl} into pandas')
        df = pd.read_csv(
            SCRIPT_PATH / ('../../data/pipeline' + tbl_models[tbl]['location']),
            sep=tbl_models[tbl]['separator'],
        )

        # Select columns specified
        df = df.loc[:, tbl_models[tbl]['col_mappings'].keys()]
        # Replace column names with the col_mappings
        df = df.rename(columns=tbl_models[tbl]['col_mappings'])

        # Strip ensembl_identifiers for the columns
        if tbl_models[tbl]['remove_ensembl_id_version_numbers']:
            df = strip_version_identifiers(df, tbl_models[tbl]['ensembl_ids'])

        # Write to SQLtables
        print(f'Writing to SQLite DB {args.db_name}')
        df.to_sql(tbl, conn, if_exists='fail', index=False, chunksize=1000)
        print(f'Ingestion for table {tbl} complete')


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
        '--single',
        action='store_true',
        help='Single table',
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='All tables',
    )
    parser.add_argument(
        '--tbl_name',
        required='single' in sys.argv,
        type=str,
        help='The name of the table to be ingested',
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
