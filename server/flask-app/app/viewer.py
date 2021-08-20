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


viewer = Blueprint('viewer', __name__)


def get_mane_genomic_features_by_ensg(ensembl_gene_id):
    """
    Get the MANE genomic features for a given gene. 

        Parameters: 
            ensembl_gene_id (str) : A stable ensembl gene identifier (ensure that this is in MANE) e.g. ENSG00000081189 
    
        Returns: 
            genomic_features (dict) : genomic features such as start, end, cds, and features which is a set of genomic records. 
    """

    # Load up MANE 
    mane = pd.read_csv("../db/MANE/0.93/MANE.genomic.csv")
    mane['ensembl_stable_gene_id']=mane['gene_id'].apply(lambda x: str(x)[0:15])
    gene_data = mane [mane["ensembl_stable_gene_id"] == ensembl_gene_id]
    genomic_features = {}
    genomic_features["gene_start"] = gene_data[gene_data["type"]=="gene"]["start"].item()
    genomic_features["gene_end"] = gene_data[gene_data["type"]=="gene"]["end"].item()
    genomic_features["start_codon"] = gene_data[gene_data["type"]=="start_codon"]["start"].item()
    genomic_features["strand"] = gene_data[gene_data["type"]=="gene"]["strand"].item()
    genomic_features["features"] = gene_data[gene_data["type"]!="gene"].to_dict('records')

    # Add sequence as well

    return genomic_features


@viewer.route("/viewer/<ensembl_transcript_id>")
def viewer_page(ensembl_transcript_id):
    # Check if transcript id is in MANE 
    gnomad_data = gnomad_search_by_transcript_id(ensembl_transcript_id)["data"]
    # This needs in a specific function to search if the gene / transcript of interest actually exists
    if gnomad_data["transcript"] == None: 
        flash("Transcript not found")
        return redirect(url_for("viewer.viewer_page", ensembl_transcript_id="ENST00000504921"))
    # Filter to 5' UTR (Make this into a specific function) for both population and clinvar variants
    gnomad_data["transcript"]["clinvar_variants"] = [clinvar for clinvar in gnomad_data["transcript"]["clinvar_variants"] if clinvar["major_consequence"] == "5_prime_UTR_variant"]
    gnomad_data["transcript"]["variants"] = [var for var in gnomad_data["transcript"]["variants"] if var["transcript_consequence"]["major_consequence"] == "5_prime_UTR_variant"]

    gene_features = get_mane_genomic_features_by_ensg("ENSG00000186575")


    return render_template("viewer.html", 
                            ensembl_transcript_id=ensembl_transcript_id, 
                            gnomad_data=gnomad_data,
                            gene_features=gene_features)