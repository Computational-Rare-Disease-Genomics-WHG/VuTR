"""Flask configuration"""


class Config:
    """Main config"""

    SECRET_KEY = 'jtg90458hgy258h02tg85u'
    FLASK_APP = 'app'
    DEBUG = True
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    VARIANT_DATABASE = 'sqlite:////../../../pipeline/database/test.db'
    FEATURES_DATABASE = 'sqlite:////../../../pipeline/database/features.db'


class DevelopmentConfig(Config):
    """Dev local config"""

    DEBUG = True
    URL = 'http://127.0.0.1:5000'
    PORT = 5000


class ProductionConfig(Config):
    """Prod config"""

    DEBUG = False


config_by_name = dict(production=ProductionConfig, development=DevelopmentConfig)
