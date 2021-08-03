# viewer.py 
# E. D'Souza 

from flask import (
    Blueprint, 
    render_template, 
    url_for, 
    redirect, 
    flash,
    jsonify,
    session 
)
import json
from .tools.gnomAD import gnomad_search_by_transcript_id, gnomad_search_by_gene_id

viewer = Blueprint('viewer', __name__)

@viewer.route("/viewer/<ensembl_transcript_id>")
def viewer_page(ensembl_transcript_id):
    # Check if transcript id is in MANE 
    gnomad_data = gnomad_search_by_transcript_id(ensembl_transcript_id)["data"]


    # This needs in a specific function to search if the gene / transcript of interest actually exists
    if gnomad_data["transcript"] == None: 
        flash("Transcript not found")
        return redirect(url_for("viewer.viewer_page", ensembl_transcript_id="ENST00000380152"))

    # Filter to 5' UTR (Make this into a specific function) for both population and clinvar variants
    gnomad_data["transcript"]["clinvar_variants"] = [clinvar for clinvar in gnomad_data["transcript"]["clinvar_variants"] if clinvar["major_consequence"] == "5_prime_UTR_variant"]
    gnomad_data["transcript"]["variants"] = [var for var in gnomad_data["transcript"]["variants"] if var["transcript_consequence"]["major_consequence"] == "5_prime_UTR_variant"]


    return render_template("viewer.html", 
                            ensembl_transcript_id=ensembl_transcript_id, 
                            gnomad_data=gnomad_data)

#@viewer.route("/viewer/<ensembl_transcript_id>/<variant>")
#def viewer_variant(ensembl_transcript_id, variant):
#    return render_template()
