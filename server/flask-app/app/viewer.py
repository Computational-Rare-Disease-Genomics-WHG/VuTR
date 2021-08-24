# viewer.py 
# E. D'Souza 
from flask import (
    Blueprint, 
    render_template, 
    url_for, 
    redirect, 
    flash)

import json
import pandas as pd

from .tools.gnomAD import gnomad_search_by_transcript_id
from .tools.mane import genomic_features_by_ensg, get_utr_stats
from .tools.sorfs import find_sorfs_by_ensg

viewer = Blueprint('viewer', __name__)

@viewer.route("/viewer/<ensembl_transcript_id>")
def viewer_page(ensembl_transcript_id):
    # Check if transcript id is in MANE 
    gnomad_data = gnomad_search_by_transcript_id(ensembl_transcript_id)
    # This needs in a specific function to search if the gene / transcript of interest actually exists
    if gnomad_data["transcript"] == None: 
        flash("Transcript not found")
        return redirect(url_for("viewer.viewer_page", ensembl_transcript_id="ENST00000504921"))
    
    
    # Filter to 5' UTR (Make this into a specific function) for both population and clinvar variants
    gnomad_data["transcript"]["clinvar_variants"] = [clinvar for clinvar in gnomad_data["transcript"]["clinvar_variants"] if clinvar["major_consequence"] == "5_prime_UTR_variant"]
    gnomad_data["transcript"]["variants"] = [var for var in gnomad_data["transcript"]["variants"] if var["transcript_consequence"]["major_consequence"] == "5_prime_UTR_variant"]

    # Find ENSG by ENST 

    ensg="ENSG00000081189"
    gene_features = genomic_features_by_ensg(ensg)
    five_prime_utr_stats = get_utr_stats(ensg)
    sorfs = find_sorfs_by_ensg(ensg)

    return render_template("viewer.html", 
                            ensembl_transcript_id=ensembl_transcript_id, 
                            gnomad_data=gnomad_data,
                            gene_features=gene_features,
                            five_prime_utr_stats=five_prime_utr_stats)