# ___init__.py

from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from .config import config_by_name

import sqlite3
# db = SQLAlchemy() # To add database when we have to put in the tables and the cache


def create_app(runtime_environment):
    app = Flask(__name__)
    # create the app through the app configuration
    app.config.from_object(config_by_name["development"])

    # TODO No database for now
    # db.init(app)

    # Variant DB as a key value store of variant consequences
    from . import variant_db
    variant_db.init_app(app)

    # Register blueprints
    from .viewer import viewer as viewer_blueprint
    app.register_blueprint(viewer_blueprint)

    from .main import main as main_blueprint
    app.register_blueprint(main_blueprint)

    return app
