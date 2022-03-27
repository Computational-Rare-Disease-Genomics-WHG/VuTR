""" ___init__.py"""

from os import environ
from flask import Flask  # pylint: disable=E0401

from .config import config_by_name

# Variant DB as a key value store of variant consequences
# Add the features database as well
from . import variant_db
from . import features_db

# Register blueprints
from .viewer import viewer as viewer_blueprint
from .main import main as main_blueprint


def create_app():
    """
    Create the flask app
    """
    app = Flask(__name__)
    # create the app through the app configuration
    app.config.from_object(config_by_name[environ.get('FLASK_ENV')])

    features_db.init_app(app)
    variant_db.init_app(app)
    app.register_blueprint(viewer_blueprint)
    app.register_blueprint(main_blueprint)

    return app
