# main.py 
# E. D'Souza 

from flask import (
    Blueprint, 
    render_template, 
    url_for, 
    redirect, 
    session 
)

import itertools
import datetime 

main = Blueprint('main', __name__)

# Home Page 
@main.route("/")
def index(): 
    return render_template("index.html")

@main.route("/index")
def index_alt(): 
    return render_template("index.html")

@main.route("/about")
def about(): 
    return render_template("about.html")

@main.route("/help")
def help(): 
    return render_template("help.html")

@main.route("/contact")
def contact(): 
    return render_template("contact.html")

@main.route("/change_log")
def change(): 
    return render_template("change_log.html")


@main.route("/download")
def download(): 
    return render_template("")