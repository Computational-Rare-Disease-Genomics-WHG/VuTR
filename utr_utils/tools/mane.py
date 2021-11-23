"""
Provides utility function to get access to the feature track data from MANE
"""

from .utils import (
    read_mane_transcript,
    read_mane_genomic_features,
    read_oorf_features
)


def get_transcript_features(ensembl_transcript_id):

    # read through the transcript sequence
    transcript_feats = {}

    transcript_entry = read_mane_transcript(ensembl_transcript_id=ensembl_transcript_id)
    transcript_feats["full_seq"] = transcript_entry["seq"].values[0]

    transcript_feats['orfs'] = read_oorf_features(
        ensembl_transcript_id).to_dict('records')

    return transcript_feats


def get_utr_stats(ensembl_gene_id):
    """
    Gets the MANE UTR Statistics for a given ENGS

    @params ensembl_gene_id (str) : A stable ensembl gene
                identifier (ensure that this is in MANE)
                e.g. ENSG00000081189
    @returns utr_stats (dict) : Statistics of the five prime utr

    """
    gene_data = read_mane_genomic_features(ensembl_gene_id)
    five_prim_utrs = gene_data[gene_data['type'] == 'five_prime_UTR']
    five_prim_utrs['width'] = five_prim_utrs.end - five_prim_utrs.start + 1
    utr_stats = {}
    utr_stats['count'] = five_prim_utrs.shape[0]
    utr_stats['5_prime_utr_length'] = sum(five_prim_utrs['width'])
    utr_stats['utr_region'] = five_prim_utrs.to_dict('records')

    return utr_stats


def get_gene_features(ensembl_gene_id):
    """
    Gets the features for a given gene by ensembl_gene_id.

    @params ensembl_gene_id (str) : A stable ensembl gene
                identifier (ensure that this is in MANE)
                 e.g. ENSG00000081189

    @returns: gene_features (dict) : A dictionary for the
     genomic features of the specified ensembl_gene_id.
    """
    gene = read_mane_genomic_features(ensembl_gene_id)
    gene_features = gene[gene['type'] == 'gene'].to_dict('records')[0]
    return gene_features


def genomic_features_by_ensg(ensembl_gene_id):
    """
    Get all MANE genomic features for a given gene.

    @params ensembl_gene_id (str) : A stable ensembl
            gene identifier (ensure that this is in MANE)
             e.g. ENSG00000081189

    @returns genomic_features (dict) : genomic features such as
            start, end, cds, and features which is a set of genomic records.
    """

    # Load up MANE
    gene_data = read_mane_genomic_features(ensembl_gene_id)
    gene_data['width'] = gene_data['end'] - gene_data['start'] + 1
    genomic_features = {}
    genomic_features['gene_start'] = gene_data[gene_data['type'] == 'gene'][
        'start'
    ].item()
    genomic_features['gene_end'] = gene_data[gene_data['type'] == 'gene']['end'].item()
    genomic_features['start_codon'] = gene_data[gene_data['type'] == 'start_codon'][
        'start'
    ].item()
    genomic_features['strand'] = gene_data[gene_data['type'] == 'gene']['strand'].item()
    genomic_features['features'] = gene_data[gene_data['type'] != 'gene'].to_dict(
        'records'
    )

    return genomic_features
