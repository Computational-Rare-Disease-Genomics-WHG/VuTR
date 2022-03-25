"""
Functions to access features store database
"""

import sqlite3
from flask import current_app, g  # pylint: disable=E0401


def init_app(app):
    """
    Init app
    """
    app.teardown_appcontext(close_db)


def dict_factory(cursor, row):
    """
    converts rows into dictionaries
    """
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d


def get_db():
    """
    Get variant database
    """
    if 'db' not in g:
        g.db = sqlite3.connect(
            current_app.config['FEATURES_DATABASE'],
            detect_types=sqlite3.PARSE_DECLTYPES,
        )
        g.db.row_factory = dict_factory
    return g.db


def close_db(e=None):  # pylint: disable=W0613
    """
    Close the db
    """
    db = g.pop('db', None)
    if db is not None:
        db.close()
