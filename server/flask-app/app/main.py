"""Main Blueprint and routes"""

from flask import (  # pylint: disable=E0401
    Blueprint,
    render_template,
    url_for,
    request,
    redirect,
)

from .helpers import convert_between_ids

main = Blueprint('main', __name__)


@main.route('/gene_search', methods=['POST'])
def gene_search():
    """Search by gene"""
    gene_name = request.form['gene_q'].upper().strip()

    ensembl_transcript_id = convert_between_ids(
        gene_name, 'hgnc_symbol', 'ensembl_transcript_id'
    )
    if ensembl_transcript_id is not None:
        return redirect(
            url_for('viewer.viewer_page', ensembl_transcript_id=ensembl_transcript_id)
        )
    return redirect(url_for('main.not_found'))


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
