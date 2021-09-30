"""
Provides utility function to get access to the feature track data from MANE
"""

from .utils import (
    find_uorfs_in_transcript,
    convert_betweeen_identifiers,
    read_mane_transcript,
    read_mane_genomic_features
)


def get_transcript_features(ensembl_transcript_id):

    # read through the transcript sequence
    transcript_feats = {}

    transcript_entry = read_mane_transcript(ensembl_transcript_id=ensembl_transcript_id)

    # get the start and  end points of the transcript.
    gene_id = convert_betweeen_identifiers(
        ensembl_transcript_id, "ensembl_transcript", "ensembl_gene")
    utr_stats = get_utr_stats(gene_id)

    transcript_feats["full_seq"] = transcript_entry["seq"].values[0]

    # find the start sites
    transcript_feats["start_site"] = utr_stats["5_prime_utr_length"]

    transcript_feats["utr_stats"] = utr_stats

    # find the stop sites

    # find all uORFs
    transcript_feats["uORF"] = find_uorfs_in_transcript(
        seq=transcript_feats["full_seq"],
        start_site=transcript_feats["start_site"],
        ensembl_transcript_id=ensembl_transcript_id)

    # find oORFS
    # transcript_feats["oORFs"] = find_oorf_in_transcript()
    # to be implemented

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
    five_prim_utrs['width'] = five_prim_utrs['end'] - five_prim_utrs['start'] + 1
    utr_stats = {}
    utr_stats['count'] = five_prim_utrs.shape[0]
    utr_stats['5_prime_utr_length'] = sum(five_prim_utrs['width'])
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

    # TODO : Add sequence as well

    ensembl_transcript_id = convert_betweeen_identifiers(
        ensembl_gene_id, "ensembl_gene", "ensembl_transcript")
    transcript_features = get_transcript_features(ensembl_transcript_id)

    # Add uORFS
    # Add oORFS

    return genomic_features
