"""
A set of sqlite3 helper functions
"""
import json

from utr_utils.tools.utils import (
    get_lookup_df,
    add_tloc_to_dict,
)

# import the datasets
from . import variant_db
from . import features_db


def convert_between_ids(from_id, from_entity, to_entity):
    """
    Converts between different entity
    """
    list_possible_cols = [
        'ensembl_transcript_id',
        'ensembl_gene_id',
        'ensembl_protein_id',
        'ncbi_gene_id',
        'refseq_transcript_id',
        'refseq_protein_id',
        'mane_status',
        'name',
        'hgnc_symbol',
        'hgnc_id',
    ]
    if from_entity in list_possible_cols and to_entity in list_possible_cols:
        rows = features_db.get_db()
        query = f"SELECT {to_entity} FROM mane_summary WHERE {from_entity}='{from_id}'"
        results = rows.execute(query).fetchone()
        features_db.close_db()
        return results[to_entity]
    return None


def get_genomic_features(ensg):
    """
    Gets the genomic features tab
    """
    cursor = features_db.get_db()

    # Query and search for results
    query = cursor.execute(
        'SELECT * FROM mane_genomic_features WHERE ensembl_gene_id=? ', [ensg]
    )
    result = query.fetchall()
    features_db.close_db()
    return result


def find_all_high_impact_utr_variants(ensembl_transcript_id):
    """
    Finds all possible UTR variants for a
     given transcript id from the database
    """
    db = variant_db.get_db()
    cursor = db.execute(
        'SELECT variant_id FROM variant_annotations WHERE ensembl_transcript_id=?',
        [ensembl_transcript_id],
    )
    rows = cursor.fetchall()
    return [i[0] for i in rows]


def get_possible_variants(ensembl_transcript_id):
    """
    Searches the database for variants
    """
    var_db = variant_db.get_db()
    cursor = var_db.execute(
        'SELECT annotations FROM variant_annotations WHERE ensembl_transcript_id =?',
        [ensembl_transcript_id],
    )
    rows = cursor.fetchall()
    variants = [json.loads(row[0]) for row in rows]
    return variants


def process_gnomad_data(gnomad_data, ensembl_transcript_id):
    """
    Get the gnomAD data and find their transcript coordinates
    and filter to 5' UTR variants
    """
    # Check if transcript id is in MANE
    glookup_table = get_lookup_df(ensembl_transcript_id=ensembl_transcript_id)
    # Filtering to SNVs for now
    gnomad_data['clinvar_variants'] = [
        add_tloc_to_dict(clinvar, glookup_table, ensembl_transcript_id)
        for clinvar in gnomad_data['clinvar_variants']
        if clinvar['major_consequence'] == '5_prime_UTR_variant'
        and len(clinvar['ref']) == 1
        and len(clinvar['alt']) == 1
    ]

    gnomad_data['variants'] = [
        add_tloc_to_dict(var, glookup_table, ensembl_transcript_id)
        for var in gnomad_data['variants']
        if var['transcript_consequence']['major_consequence'] == '5_prime_UTR_variant'
        and len(var['ref']) == 1
        and len(var['alt']) == 1
    ]
    gnomad_variants_list = [var['variant_id'] for var in gnomad_data['variants']]

    clinvar_variants_list = [
        var['variant_id'] for var in gnomad_data['clinvar_variants']
    ]

    return gnomad_data, gnomad_variants_list, clinvar_variants_list


def get_transcript_features(ensembl_transcript_id):
    """
    Get transcript features
    """
    db = features_db.get_db()
    cursor = db.execute(
        'SELECT * FROM mane_transcript_features WHERE ensembl_transcript_id=?',
        [ensembl_transcript_id],
    )
    rows = cursor.fetchone()
    features_db.close_db()
    return rows


def get_genome_to_transcript_intervals(ensembl_transcript_id, tpos):
    """
    Retrieves all of the features of the uorfs / uorfs
    for the native architechure of the gene
    """
    db = features_db.get_db()
    cursor = db.execute(
        'SELECT genomic_pos FROM genome_to_transcript_coordinates WHERE ensembl_transcript_id=? AND transcript_pos=?',  # noqa: E501 # pylint: disable=C0301
        [ensembl_transcript_id, tpos],
    )
    result = cursor.fetchone()
    features_db.close_db()
    return result['genomic_pos']


def get_all_orfs_features(ensembl_transcript_id):
    """
    Retrieves all of the features of the uorfs / uorfs
    for the native architechure of the gene
    """
    db = features_db.get_db()
    cursor = db.execute(
        'SELECT * FROM orf_features WHERE ensembl_transcript_id=?',
        [ensembl_transcript_id],
    )
    rows = cursor.fetchall()
    features_db.close_db()

    # Add genomic coordinates : TODO Optimize this in the
    # by adding this in the pipeline
    for u in rows:
        u['uorf_start_genome'] = get_genome_to_transcript_intervals(
            u['ensembl_transcript_id'], u['orf_start_codon']
        )
        u['uorf_stop_genome'] = get_genome_to_transcript_intervals(
            u['ensembl_transcript_id'], u['orf_stop_codon']
        )

    return rows


def get_constraint_score(ensembl_gene_id):
    """
    Get constraint score from the features db
    @param ensembl_gene_id
    @returns constraint score (double)
    """
    db = features_db.get_db()
    cursor = db.execute(
        'SELECT loeuf FROM loeuf_constraint WHERE ensembl_gene_id=?',
        [ensembl_gene_id],
    )
    result = cursor.fetchone()
    features_db.close_db()
    return result['loeuf']
