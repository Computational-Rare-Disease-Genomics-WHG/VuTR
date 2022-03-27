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
    get_gnomad_variants_in_utr_regions,
)

# from utr_utils.tools.sorfs import find_sorfs_by_ensg
from utr_utils.tools.clingen import get_clingen_curation
from utr_utils.tools.utr_annotation import (
    get_utr_annotation_for_list_variants,
    find_intervals_for_utr_consequence,
)
from .helpers import (
    get_possible_variants,
    process_gnomad_data,
    find_all_high_impact_utr_variants,
    get_genomic_features,
    get_transcript_features,
    convert_between_ids,
    get_constraint_score,
    get_all_orfs_features,
)

from . import variant_db

viewer = Blueprint('viewer', __name__)


@viewer.route('/viewer/utr_impact', methods=['GET', 'POST'])
def get_utr_impacts():
    """
    A JSON API resource to get the 5' UTR annotation for a supplied variant
    @param variant_id e.g. 5-150904976-T-A
    @param ensembl_transcript_id e.g. ENST00000274599
    """
    variant_id = request.args['variant_id']
    ensembl_transcript_id = request.args['ensembl_transcript_id']
    start_site = request.args['start_site']
    buffer = request.args['buffer']
    try:
        db = variant_db.get_db()
        cursor = db.execute(
            'SELECT annotations, five_prime_UTR_variant_annotation FROM variant_annotations WHERE ensembl_transcript_id =? AND variant_id=?',  # noqa: E501 # pylint: disable=C0301
            [ensembl_transcript_id, variant_id],
        )
        rows = cursor.fetchall()
        variant = [json.loads(row[0]) for row in rows][0]
        annotation = [json.loads(row[1]) for row in rows][0]
        response_object = {
            'status': 'Success',
            'message': 'Ok',
            'data': {
                'variant': {**{'variant_id': variant['variant_id']}, **annotation},
                'intervals': find_intervals_for_utr_consequence(
                    var_id=variant['variant_id'],
                    conseq_type=variant['five_prime_UTR_variant_consequence'],
                    conseq_dict=variant['five_prime_UTR_variant_annotation'],
                    cdna_pos=variant['cDNA_position'],
                    start_site=start_site,
                    buffer_length=buffer,
                ),
            },
        }
        return response_object, 200
    except Exception as e:  # pylint: disable=W0703
        response_object = {
            'status': 'Failure',
            'message': f'Unable to fetch utr consequence error {str(e)}',
        }
        return response_object, 400


@viewer.route('/viewer/<ensembl_transcript_id>')
def viewer_page(ensembl_transcript_id):
    """
    Collects data for a given ENST
    @param ensembl_transcript_id
    """
    buffer = 40

    # Find ENSG by ENST
    ensembl_gene_id = convert_between_ids(
        ensembl_transcript_id, 'ensembl_transcript_id', 'ensembl_gene_id'
    )
    hgnc = convert_between_ids(
        ensembl_transcript_id, 'ensembl_transcript_id', 'hgnc_symbol'
    )
    name = convert_between_ids(ensembl_transcript_id, 'ensembl_transcript_id', 'name')
    refseq_match = convert_between_ids(
        ensembl_transcript_id, 'ensembl_transcript_id', 'refseq_transcript_id'
    )

    # Get features
    gene_features = get_genomic_features(ensembl_gene_id)
    five_prime_utr_stats = get_transcript_features(ensembl_transcript_id)
    transcript_features = get_all_orfs_features(ensembl_transcript_id)

    # sorfs = find_sorfs_by_ensg(ensembl_gene_id)
    constraint = get_constraint_score(ensembl_gene_id)
    clingen_curation_record = get_clingen_curation(hgnc)
    start_site = five_prime_utr_stats['start_site_pos']

    possible_variants = get_possible_variants(
        ensembl_transcript_id=ensembl_transcript_id
    )

    # Find UTR regions
    utr_regions = [i for i in gene_features if i['type'] == 'five_prime_UTR']

    # This needs in a specific functon
    # to search if the gene / transcript
    # of interest actually exists
    gnomad_data, gnomad_variants_list, clinvar_variants_list = process_gnomad_data(
        get_gnomad_variants_in_utr_regions(utr_regions=utr_regions),
        ensembl_transcript_id,
    )

    # Get the annotations for these values.
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
        gene_features=gene_features,
        five_prime_utr_stats=five_prime_utr_stats,
        transcript_features=transcript_features,
        gnomad_utr_impact=gnomad_utr_impact,
        clinvar_utr_impact=clinvar_utr_impact,
        all_possible_variants=all_possible_variants,
    )
