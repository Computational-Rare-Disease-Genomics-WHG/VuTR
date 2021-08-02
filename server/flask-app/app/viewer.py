# viewer.py 
# E. D'Souza 

from flask import (
    Blueprint, 
    render_template, 
    url_for, 
    redirect, 
    session 
)
import json
from .tools.gnomAD import gnomad_search_by_transcript_id, gnomad_search_by_gene_id

viewer = Blueprint('viewer', __name__)

@viewer.route("/viewer/<ensembl_transcript_id>")
def viewer_page(ensembl_transcript_id):
    print(gnomad_search_by_transcript_id("ENST00000357985")["data"]["transcript"])
    return render_template("viewer.html", ensembl_transcript_id=ensembl_transcript_id)

#@viewer.route("/viewer/<ensembl_transcript_id>/<variant>")
#def viewer_variant(ensembl_transcript_id, variant):
#    return render_template()
