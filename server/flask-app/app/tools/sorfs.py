# sorfs.py

# Access to sorfs.org data
import pandas as pd 
from .mane import get_gene_features

def read_sorfs():
    """
    Reading the sorfs downloaded from the SOAP API in the pipeline 
        
        Parameters :
            None 
        Returns : 
            data (Dataframe) : A dataframe with the list of sorfs. 
    """

    # TODO : Change this to a SQL Alchemy model
    data = pd.read_csv("../db/SORFS/sorfs.csv", sep=",")
    return data

def find_sorfs_by_ensg(ensembl_gene_id):
    """
    Filters the sorfs datasets to the sorfs that appear in the 5' utrs of the given ENSG. 

        Parameters : 
            ensembl_gene_id (str) : A stable ensembl gene identifier (ensure that this is in MANE)
             e.g. ENSG00000081189 

        Returns : 
            filtered_sorfs (dict) : A list of dictionaries each with a sorf within the genomic coordinates of the ensembl_gene_id
    
    """
    sorfs_db = read_sorfs()
    gene_features = get_gene_features(ensembl_gene_id)

    # TODO: Mapping stranges from the two different formats used by the datasets.
    strand_mapping = {
        "+" : 1,
        "-" : -1  
    }

    # TODO : Change this to a SQL Alchemy model filtering query

    ## Filter to genomic location 
    filtered_sorfs = sorfs_db[(sorfs_db["Chromosome"] == str(gene_features["seqid"])[3:]) & 
                    (sorfs_db["Sorf end"] < gene_features["end"]) & 
                    (sorfs_db["Sorf start"] > gene_features["start"]) &
                    (sorfs_db["Strand"]== strand_mapping[gene_features["strand"]])]


    return filtered_sorfs.to_dict("records")