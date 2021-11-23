"""
Takes a table and ingests it into the sqlite3 database

# Inspiration 
https://towardsdatascience.com/dramatically-improve-your-database-inserts-with-a-simple-upgrade-6dfa672f1424

"""


from model import tbl_queries, fils_locs
import dask.dataframe as dd
import sqlite3
import sqlalchemy


# TODO
