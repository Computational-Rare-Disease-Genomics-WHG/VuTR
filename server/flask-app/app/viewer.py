# viewer.py 
# E. D'Souza 

from flask import (
    Blueprint, 
    render_template, 
    url_for, 
    redirect, 
    jsonify,
    session 
)
import json
from .tools.gnomAD import gnomad_search_by_transcript_id, gnomad_search_by_gene_id

viewer = Blueprint('viewer', __name__)

@viewer.route("/viewer/<ensembl_transcript_id>")
def viewer_page(ensembl_transcript_id):
    # Check if transcript id is in MANE 

    # Filter to 5' UTR

    gnomad_data = gnomad_search_by_transcript_id(ensembl_transcript_id)["data"]
    return render_template("viewer.html", 
                            ensembl_transcript_id=ensembl_transcript_id, 
                            gnomad_data=gnomad_data)

#@viewer.route("/viewer/<ensembl_transcript_id>/<variant>")
#def viewer_variant(ensembl_transcript_id, variant):
#    return render_template()
