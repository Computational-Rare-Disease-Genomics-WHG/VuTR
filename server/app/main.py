# main.py 
# E. D'Souza 



import os 

base_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

class Config : 
    SECRET_KEY=""
    FLASK_APP="app"
    DEBUG=True
    SQLALCHEMY_TRACK_MODIFICATIONS=False



class DevelopmentConfig (Config): 
    DEBUG=True 
    SQLALCHEMY_DATABASE_URI="sqlite:///"+os.path.join(base_dir, 'data.db')
    URL="http://127.0.0.1:5000"
    PORT=5000
    LOCAL_IP="127.0.0.1"


class ProductionConfig(Config):
    DEBUG=True 
    SQLALCHEMY_DATABASE_URI="sqlite:////data.db" # Hosted on docker container
    URL="http://127.0.0.1:5000"
    PORT=5000
    LOCAL_IP="127.0.0.1"


config_by_name = dict(production = ProductionConfig, development=DevelopmentConfig)