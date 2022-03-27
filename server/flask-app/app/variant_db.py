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
    if 'vdb' not in g:
        g.vdb = sqlite3.connect(
            current_app.config['VARIANT_DATABASE'], detect_types=sqlite3.PARSE_DECLTYPES
        )
        g.vdb.row_factory = sqlite3.Row
    return g.vdb


def close_db(e=None):  # pylint: disable=W0613
    """
    Close the db
    """
    vdb = g.pop('vdb', None)
    if vdb is not None:
        vdb.close()
