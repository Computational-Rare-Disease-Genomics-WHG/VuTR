"""
Functions to access variant store database
"""


import sqlite3
from flask import current_app, g  # pylint: disable=E0401


def init_app(app):
    """
    Init app
    """
    app.teardown_appcontext(close_db)


def get_db():
    """
    Get variant database
    """
    if 'fdb' not in g:
        g.fdb = sqlite3.connect(
            current_app.config['VARIANT_DATABASE'], detect_types=sqlite3.PARSE_DECLTYPES
        )
        g.fdb.row_factory = sqlite3.Row
    return g.fdb


def close_db(e=None):  # pylint: disable=W0613
    """
    Close the db
    """
    fdb = g.pop('fdb', None)
    if fdb is not None:
        fdb.close()
