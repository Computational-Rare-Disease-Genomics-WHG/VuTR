# 


import pandas as pd # Replace with database

def get_utrs_with_evidence(hgnc):
    data = pd.read_csv("../../../db/uORF_5UTR_GRCh38_PUBLIC.txt" )  # Replace with database
    return data[data["gene"] ==hgnc]