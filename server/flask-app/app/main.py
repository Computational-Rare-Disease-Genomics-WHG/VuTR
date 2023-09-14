"""Main Blueprint and routes"""

import re
from flask import (  # pylint: disable=E0401
    Blueprint,
    render_template,
    url_for,
    request,
    redirect,
    current_app
)

from .helpers import (
    convert_between_ids,
    search_enst_by_transcript_id,
    get_transcript_features,
    get_genomic_features,
    find_transcript_ids_by_gene_id
)

main = Blueprint('main', __name__)


@main.route('/gene_search', methods=['POST'])
def gene_search():
    """Search by gene"""
    query = request.form['gene_q'].upper().strip()

    # detect query type using grep
    if re.search(r'\d-\d+-\w-\w', query):
        return redirect(url_for('main.search_variant', variant=query))

    elif query[0:4] == 'ENSG':
        ensembl_transcript_id = convert_between_ids(
        query, 'ensembl_gene_id', 'ensembl_transcript_id'
    )   
    else: 
        ensembl_transcript_id = convert_between_ids(
            query, 'hgnc_symbol', 'ensembl_transcript_id'
        )
        
    
    if ensembl_transcript_id is not None:
        return redirect(
            url_for('viewer.viewer_page', ensembl_transcript_id=ensembl_transcript_id)
        )
    return redirect(url_for('main.not_found'))


@main.route('/variant_search/<variant>')
def search_variant(variant):
    """
    Route to search by variant to see all of the transcripts
    """
    # Search the list of transcript ids that this variant
    # falls under
    # Extract Genomic position
    variant_list = search_enst_by_transcript_id(variant)
    variant_dat_list = []
    for enst in variant_list :
        print(enst)
        variant_dat = {}
        variant_dat['ensembl_transcript_id'] = enst['ensembl_transcript_id']
        variant_dat['ensembl_gene_id'] = convert_between_ids(enst['ensembl_transcript_id'],
        'ensembl_transcript_id', 'ensembl_gene_id')
        variant_dat['gene_features'] = get_genomic_features(
            variant_dat['ensembl_gene_id'])
        variant_dat['five_prime_utr_stats'] = get_transcript_features(enst['ensembl_transcript_id'])
        variant_dat_list.append(variant_dat)

    impact_url = current_app.config['IMPACT_URL']
    return render_template(
        'variant.html',
        variant_list=variant_list,
        variant=variant,
        variant_dat_list=variant_dat_list,
        impact_url=impact_url
    )


@main.route('/gene_not_found')
def not_found():
    """Not found page"""
    return render_template('404.html')


@main.route('/')
def index():
    """Index home page"""
    return render_template('index.html')


@main.route('/about')
def about():
    """About page"""
    return render_template('about.html')


@main.route('/help')
def help_page():
    """Help Page"""
    return render_template('help.html')


@main.route('/changelog')
def changelog():
    """Change Log"""
    return render_template('changelog.html')
