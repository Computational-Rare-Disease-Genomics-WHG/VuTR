"""
Flask blueprint to define the core viewer page.
"""

from flask import (  # pylint: disable=E0401
    Blueprint,
    render_template,
    url_for,
    redirect,
    flash,
)  # pylint disable=E0401

# import from tools
from .tools.gnomad import get_constraint_by_ensg, gnomad_search_by_transcript_id
from .tools.mane import genomic_features_by_ensg, get_utr_stats
from .tools.sorfs import find_sorfs_by_ensg
from .tools.clingen import get_clingen_curation

viewer = Blueprint('viewer', __name__)


@viewer.route('/viewer/<ensembl_transcript_id>')
def viewer_page(ensembl_transcript_id):
    """
    Collects data for a given ENST
    @param ensembl_transcript_id
    """
    # Check if transcript id is in MANE
    gnomad_data = gnomad_search_by_transcript_id(ensembl_transcript_id)
    # This needs in a specific functon
    # to search if the gene / transcript
    # of interest actually exists
    if gnomad_data['transcript'] is None:
        flash('Transcript not found')
        return redirect(
            url_for('viewer.viewer_page', ensembl_transcript_id='ENST00000504921')
        )

    # Filter to 5' UTR (Make this into
    # a specific function) for both population and clinvar variants
    gnomad_data['transcript']['clinvar_variants'] = [
        clinvar
        for clinvar in gnomad_data['transcript']['clinvar_variants']
        if clinvar['major_consequence'] == '5_prime_UTR_variant'
    ]
    gnomad_data['transcript']['variants'] = [
        var
        for var in gnomad_data['transcript']['variants']
        if var['transcript_consequence']['major_consequence'] == '5_prime_UTR_variant'
    ]

    # Find ENSG by ENST

    ensg = 'ENSG00000081189'
    hgnc = 'MEF2C'
    gene_features = genomic_features_by_ensg(ensg)
    five_prime_utr_stats = get_utr_stats(ensg)
    sorfs = find_sorfs_by_ensg(ensg)
    constraint = get_constraint_by_ensg(ensg)
    clingen_curation_record = get_clingen_curation(hgnc)

    return render_template(
        'viewer.html',
        ensembl_transcript_id=ensembl_transcript_id,
        gnomad_data=gnomad_data,
        constraint=constraint,
        clingen_curation_record=clingen_curation_record,
        sorfs=sorfs,
        gene_features=gene_features,
        five_prime_utr_stats=five_prime_utr_stats,
    )
