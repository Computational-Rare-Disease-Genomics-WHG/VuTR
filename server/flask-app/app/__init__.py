# ___init__.py 

from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from .confg import config_by_name 


db = SQLAlchemy()
runtime_app_config = config_by_name["development"]

def create_app(runtime_environment): 
    app = Flask(__name__)

    # create the app through the app configuration 
    runtime_app_config = config_by_name["development" if runtime_environment == "production" else "production"]
    app.config.from_object(runtime_app_config)
    db.init(app)


    return app
