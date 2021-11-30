# main.py
# E. D'Souza

from flask import (
    Blueprint,
    render_template,
    url_for,
    request,
    jsonify, render_template,
    redirect)

from utr_utils.tools.utils import convert_betweeen_identifiers

main = Blueprint('main', __name__)


@main.route("/search", methods=["POST"])
def search():
    # TODO
    pass


@main.route("/gene_search", methods=["POST"])
def gene_search():
    gene_name = request.form["gene_q"]
    ensembl_transcript_id = convert_betweeen_identifiers(
        gene_name,
        "hgnc_symbol",
        "ensembl_transcript")
    print(ensembl_transcript_id)
    ensembl_transcript_id = ensembl_transcript_id[0:15]
    return redirect(url_for('viewer.viewer_page', ensembl_transcript_id=ensembl_transcript_id))


@main.route("/")
def index():
    return render_template("index.html")
    # return redirect(url_for("viewer.viewer_page", ensembl_transcript_id="ENST00000504921"))


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
