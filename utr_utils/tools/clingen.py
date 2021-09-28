"""
Provides access to ClinGen data
"""
import pandas as pd


def get_clingen_curation(hgnc):
    """
    Searches through clingen for the haploinsufficency score

    @param hgnc (str) : HGNC Gene Symbol
    @returns clingen_curation_record (dict) : The dictionary\
         associated with the clingen record
    """
    clingen_curation = pd.read_csv(
        '../../data/pipeline/CLINGEN/ClinGen_gene_curation_list_GRCh38.tsv',
        sep='\t',
        skiprows=5,
    )
    if hgnc in clingen_curation['#Gene Symbol']:
        clingen_curation_record = clingen_curation[
            clingen_curation['#Gene Symbol'] == hgnc
        ].to_dict('records')[0]
        return clingen_curation_record
    return None
