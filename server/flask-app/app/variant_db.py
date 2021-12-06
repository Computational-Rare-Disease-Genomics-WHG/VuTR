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
    if 'db' not in g:
        g.db = sqlite3.connect(
            current_app.config['VARIANT_DATABASE'], detect_types=sqlite3.PARSE_DECLTYPES
        )
        g.db.row_factory = sqlite3.Row
    return g.db


def close_db(e=None):  # pylint: disable=W0613
    """
    Close the db
    """
    db = g.pop('db', None)
    if db is not None:
        db.close()
