import sqlite3

from flask import current_app, g


def init_app(app):
    app.teardown_appcontext(close_db)


def get_db():
    if 'db' not in g:
        g.db = sqlite3.connect(
            current_app.config['VARIANT_DATABASE'],
            detect_types=sqlite3.PARSE_DECLTYPES
        )
        g.db.row_factory = sqlite3.Row

    return g.db


def close_db(e=None):
    db = g.pop('db', None)

    if db is not None:
        db.close()
