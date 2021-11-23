# config.py
# E D'Souza

import os

# To be implemented from pipeline
database = "sqlite:///../../db/db.sqlite"


class Config:
    SECRET_KEY = "jtg90458hgy258h02tg85u"
    FLASK_APP = "app"
    DEBUG = True
    SQLALCHEMY_TRACK_MODIFICATIONS = False


class DevelopmentConfig (Config):
    DEBUG = True
    SQLALCHEMY_DATABASE_URI = database
    URL = "http://127.0.0.1:5000"
    VARIANT_DATABASE = 'sqlite:////../../../pipeline/database/test.db'
    PORT = 5000


class ProductionConfig(Config):
    DEBUG = True
    SQLALCHEMY_DATABASE_URI = "sqlite:////data.db"  # Hosted on docker container
    DEBUG = False
    VARIANT_DATABASE = 'sqlite:////../../../pipeline/database/test.db'


config_by_name = dict(production=ProductionConfig, development=DevelopmentConfig)
