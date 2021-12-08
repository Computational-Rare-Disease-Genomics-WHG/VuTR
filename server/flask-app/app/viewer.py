"""
Flask blueprint to define the core viewer page.
"""

import json
from flask import (  # pylint: disable=E0401
    Blueprint,
    render_template,
    request,
)  # pylint disable=E0401

# import from the packages
from utr_utils.tools.gnomad import (
    get_constraint_by_ensg,
    get_gnomad_variants_in_utr_regions,
)
from utr_utils.tools.mane import (
    genomic_features_by_ensg,
    get_transcript_features,
    get_utr_stats,
)
from utr_utils.tools.sorfs import find_sorfs_by_ensg
from utr_utils.tools.clingen import get_clingen_curation
from utr_utils.tools.utils import (
    convert_betweeen_identifiers,
    get_lookup_df,
    add_tloc_to_dict,
)
from utr_utils.tools.utr_annotation import get_utr_annotation_for_list_variants

from . import variant_db

viewer = Blueprint('viewer', __name__)


def find_all_high_impact_utr_variants(ensembl_transcript_id):
    """
    Finds all possible UTR variants for a given transcript id from the database
    """
    db = variant_db.get_db()
    cursor = db.execute(
        'SELECT variant_id FROM variant_annotations WHERE ensembl_transcript_id=?',
        [ensembl_transcript_id],
    )
    rows = cursor.fetchall()
    return [i[0] for i in rows]


@viewer.route('/viewer/utr_impact', methods=['GET', 'POST'])
def get_utr_impacts():
    """
    A JSON API resource to get the 5' UTR annotation for a supplied variant
    @param variant_id e.g. 5-150904976-T-A
    @param ensembl_transcript_id e.g. ENST00000274599
    """
    variant_id = request.args['variant_id']
    ensembl_transcript_id = request.args['ensembl_transcript_id']
    try:
        db = variant_db.get_db()
        cursor = db.execute(
            'SELECT annotations FROM variant_annotations'
            'WHERE ensembl_transcript_id =? AND variant_id=?',
            [ensembl_transcript_id, variant_id],
        )
        rows = cursor.fetchall()
        variants = [json.loads(row[0]) for row in rows]
        response_object = {'status': 'Success', 'message': 'Ok', 'data': variants}
        return response_object, 200
    except Exception as e:  # pylint: disable=W0703
        response_object = {
            'status': 'Failure',
            'message': f'Unable to fetch utr consequence error {str(e)}',
        }
        return response_object, 400


def get_possible_variants(ensembl_transcript_id):
    """
    Searches the database for variants
    """
    db = variant_db.get_db()
    cursor = db.execute(
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


@viewer.route('/viewer/<ensembl_transcript_id>')
def viewer_page(ensembl_transcript_id):
    """
    Collects data for a given ENST
    @param ensembl_transcript_id
    """

    # Find ENSG by ENST
    ensembl_gene_id = convert_betweeen_identifiers(
        ensembl_transcript_id, 'ensembl_transcript', 'ensembl_gene'
    )
    hgnc = convert_betweeen_identifiers(
        ensembl_transcript_id, 'ensembl_transcript', 'hgnc_symbol'
    )
    name = convert_betweeen_identifiers(
        ensembl_transcript_id, 'ensembl_transcript', 'name'
    )
    refseq_match = convert_betweeen_identifiers(
        ensembl_transcript_id, 'ensembl_transcript', 'refseq_mrna'
    )

    # Get features
    gene_features = genomic_features_by_ensg(ensembl_gene_id)
    five_prime_utr_stats = get_utr_stats(ensembl_gene_id)
    transcript_features = get_transcript_features(ensembl_transcript_id)

    sorfs = find_sorfs_by_ensg(ensembl_gene_id)
    constraint = get_constraint_by_ensg(ensembl_gene_id)
    clingen_curation_record = get_clingen_curation(hgnc)
    buffer = 140
    start_site = five_prime_utr_stats['5_prime_utr_length'] + 1

    possible_variants = get_possible_variants(
        ensembl_transcript_id=ensembl_transcript_id
    )

    # This needs in a specific functon
    # to search if the gene / transcript
    # of interest actually exists
    gnomad_data, gnomad_variants_list, clinvar_variants_list = process_gnomad_data(
        get_gnomad_variants_in_utr_regions(five_prime_utr_stats['utr_region']),
        ensembl_transcript_id,
    )

    # Filter to 5' UTR (Make this into
    # a specific function) for both population and clinvar variants

    gnomad_utr_impact = get_utr_annotation_for_list_variants(
        gnomad_variants_list, possible_variants, start_site, buffer
    )
    clinvar_utr_impact = get_utr_annotation_for_list_variants(
        clinvar_variants_list, possible_variants, start_site, buffer
    )

    #
    all_possible_variants = find_all_high_impact_utr_variants(
        ensembl_transcript_id=ensembl_transcript_id
    )
    print(all_possible_variants)
    # Render template
    return render_template(
        'viewer.html',
        ensembl_transcript_id=ensembl_transcript_id,
        ensembl_gene_id=ensembl_gene_id,
        hgnc=hgnc,
        name=name,
        refseq_match=refseq_match,
        gnomad_data=gnomad_data,
        constraint=constraint,
        clingen_curation_record=clingen_curation_record,
        sorfs=sorfs,
        gene_features=gene_features,
        five_prime_utr_stats=five_prime_utr_stats,
        transcript_features=transcript_features,
        gnomad_utr_impact=gnomad_utr_impact,
        clinvar_utr_impact=clinvar_utr_impact,
        all_possible_variants=all_possible_variants,
    )
