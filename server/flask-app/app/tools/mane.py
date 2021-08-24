# mane.py 
# E. D'Souza 

# Provides utility function to get access to the feature track data from MANE.

from Bio import SeqIO
from gtfparse import read_gtf
import gffutils
import pandas as pd 

def read_mane (ensembl_gene_id) : 
    # TODO : This needs to be updated to Pipeline
    mane = pd.read_csv("../db/MANE/0.93/MANE.genomic.csv")
    mane['ensembl_stable_gene_id']=mane['gene_id'].apply(lambda x: str(x)[0:15])
    gene_data = mane [mane["ensembl_stable_gene_id"] == ensembl_gene_id]
    return gene_data

def get_utr_stats (ensembl_gene_id):
    """
    Gets the MANE UTR Statistics for a given ENGS
        Paramters: 
            ensembl_gene_id (str) : A stable ensembl gene identifier (ensure that this is in MANE) e.g. ENSG00000081189 
        Returns: 
            utr_stats (dict) : Statistics of the five prime utr
    
    """
    gene_data = read_mane(ensembl_gene_id)
    five_prim_utrs = gene_data[gene_data["type"]=="five_prime_UTR"]
    five_prim_utrs["width"] =five_prim_utrs["end"]-five_prim_utrs["start"]+1
    utr_stats = {}
    utr_stats["count"] = five_prim_utrs.shape[0]
    utr_stats["5_prime_utr_length"] = sum(five_prim_utrs["width"])
    return utr_stats

def genomic_features_by_ensg(ensembl_gene_id):
    """
    Get the MANE genomic features for a given gene. 

        Parameters: 
            ensembl_gene_id (str) : A stable ensembl gene identifier (ensure that this is in MANE) e.g. ENSG00000081189 
    
        Returns: 
            genomic_features (dict) : genomic features such as start, end, cds, and features which is a set of genomic records. 
    """

    # Load up MANE 
    gene_data = read_mane(ensembl_gene_id)
    gene_data["width"] = gene_data["end"]-gene_data["start"]+1
    genomic_features = {}
    genomic_features["gene_start"] = gene_data[gene_data["type"]=="gene"]["start"].item()
    genomic_features["gene_end"] = gene_data[gene_data["type"]=="gene"]["end"].item()
    genomic_features["start_codon"] = gene_data[gene_data["type"]=="start_codon"]["start"].item()
    genomic_features["strand"] = gene_data[gene_data["type"]=="gene"]["strand"].item()
    genomic_features["features"] = gene_data[gene_data["type"]!="gene"].to_dict('records')

    # TODO : Add sequence as well
    return genomic_features