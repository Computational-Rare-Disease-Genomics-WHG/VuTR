# ___init__.py 

from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from .config import config_by_name 

# db = SQLAlchemy() # To add database when we have to put in the tables and the cache

def create_app(runtime_environment): 
    app = Flask(__name__)
    # create the app through the app configuration 
    print(runtime_environment)
    app.config.from_object(config_by_name["development"])
    
    # TODO No database for now
    #db.init(app)

    # Register blueprints 
    from .viewer import viewer as viewer_blueprint
    app.register_blueprint(viewer_blueprint)

    from .main import main as main_blueprint 
    app.register_blueprint(main_blueprint)

    return app
