"""
Flask blueprint to define the core viewer page.
"""

import json
from sqlite3 import Error as SQLiteError  # pylint: disable=E0401
from flask import (  # pylint: disable=E0401
    Blueprint,
    render_template,
    request,
    jsonify,
    current_app,
)  # pylint disable=E0401

from .helpers import (
    get_possible_variants,
    process_gnomad_data,
    find_all_high_impact_utr_variants,
    get_genomic_features,
    get_transcript_features,
    convert_between_ids,
    get_smorfs,
    get_constraint_score,
    get_all_orfs_features,
    get_utr_annotation_for_list_variants,
    find_intervals_for_utr_consequence,
    get_gnomad_variants_in_utr_regions,
    get_omim_id,
    get_clingen_entry,
)
from . import variant_db


viewer = Blueprint("viewer", __name__)


def get_first_variants(ensembl_transcript_id):
    """
    Get the first 5 variants for a given ENST
    @param ensembl_transcript_id
    @returns a list of the first 5 variants
    """
    db = variant_db.get_db()
    cursor = db.execute(
        """
            SELECT variant_id
            FROM variant_annotations
            WHERE ensembl_transcript_id = ?
            LIMIT 5;
        """,
        [ensembl_transcript_id],
    )
    rows = [u['variant_id'] for u in cursor.fetchall()]
    return rows



@viewer.route("/viewer/possible_variants", methods=["GET"])
def get_possible_variants_api():
    """
    A JSON API resource to get the possible variants for a supplied search query
    @param search_term e.g. 5-150904976-T-A
    @param ensembl_transcript_id e.g. ENST00000274599
    @returns a list of matches
    """
    try:
        search_term = request.args.get("search_term")
        ensembl_transcript_id = request.args.get("ensembl_transcript_id")

        if not search_term or not ensembl_transcript_id:
            return (
                jsonify(
                    {
                        "message": "Both search_term and ensembl_transcript_id are required.",
                        "data": [],
                    }
                ),
                400,
            )

        db = variant_db.get_db()  # pylint: disable=C0103
        cursor = db.execute(
            """
                SELECT variant_id
                FROM variant_annotations
                WHERE variant_id LIKE ? COLLATE NOCASE
                AND ensembl_transcript_id = ?
                LIMIT 10;
            """,
            ["%" + search_term + "%", ensembl_transcript_id],
        )
        rows = cursor.fetchall()
        # Convert row objects to dictionaries
        rows_as_dict = [{'text': row['variant_id'], 'id' : row['variant_id']} for row in rows]
        response_object = {
            "message": "Ok",
            "data": rows_as_dict,
        }
        return jsonify(response_object), 200
    except SQLiteError as error:
        return (
            jsonify(
                {
                    "message": "Database error occurred: {}".format(str(error)),
                    "data": [],
                }
            ),
            500,
        )
    except Exception as error:  # pylint: disable=W0703
        return (
            jsonify(
                {
                    "message": "An unexpected error occurred: {}".format(str(error)),
                    "data": [],
                }
            ),
            500,
        )


@viewer.route("/viewer/utr_impact", methods=["GET"])
def get_utr_impacts():
    """
    A JSON API resource to get the 5' UTR annotation for a supplied variant
    @param variant_id e.g. 5-150904976-T-A
    @param ensembl_transcript_id e.g. ENST00000274599
    """
    try:
        variant_id = request.args.get("variant_id")
        ensembl_transcript_id = request.args.get("ensembl_transcript_id")
        start_site = request.args.get("start_site")
        buffer = request.args.get("buffer")

        if not all([variant_id, ensembl_transcript_id, start_site, buffer]):
            return (
                jsonify(
                    {
                        "message": "All parameters variant_id, ensembl_transcript_id, start_site and buffer are required.",
                        "data": {},
                    }
                ),
                400,
            )

        db = variant_db.get_db()
        cursor = db.execute(
            """
                SELECT annotations, five_prime_UTR_variant_annotation 
                FROM variant_annotations 
                WHERE ensembl_transcript_id = ? AND variant_id = ?
            """,
            [ensembl_transcript_id, variant_id],
        )
        rows = cursor.fetchall()

        if not rows:
            return (
                jsonify(
                    {
                        "message": "No data found for the given variant_id and ensembl_transcript_id",
                        "data": {},
                    }
                ),
                404,
            )

        variant = [json.loads(row[0]) for row in rows][0]
        annotation = [json.loads(row[1]) for row in rows][0]
        intervals = find_intervals_for_utr_consequence(
            var_id=variant_id,
            conseq_type=variant["five_prime_UTR_variant_consequence"],
            conseq_dict=variant["five_prime_UTR_variant_annotation"],
            cdna_pos=variant["cDNA_position"],
            start_site=start_site,
            buffer_length=buffer,
            annotation_id=variant_id,
        )

        response_object = {
            "message": "Ok",
            "data": {
                "variant": {**{"variant_id": variant_id}, **annotation},
                "intervals": intervals,
            },
        }

        return jsonify(response_object), 200

    except SQLiteError as error:
        return (
            jsonify(
                {"message": "Database error occurred: {}".format(str(error)), "data": {}}
            ),
            500,
        )

    except Exception as error:  # pylint: disable=W0703
        return (
            jsonify(
                {
                    "message": "An unexpected error occurred: {}".format(str(error)),
                    "data": {},
                }
            ),
            500,
        )


@viewer.route("/viewer/<ensembl_transcript_id>")
def viewer_page(ensembl_transcript_id):
    """
    Collects data for a given ENST
    @param ensembl_transcript_id
    """
    buffer = 40

    # Convert IDs and get basic data
    ensembl_gene_id = convert_between_ids(
        ensembl_transcript_id, "ensembl_transcript_id", "ensembl_gene_id"
    )
    hgnc = convert_between_ids(
        ensembl_transcript_id, "ensembl_transcript_id", "hgnc_symbol"
    )
    name = convert_between_ids(ensembl_transcript_id, "ensembl_transcript_id", "name")
    refseq_match = convert_between_ids(
        ensembl_transcript_id, "ensembl_transcript_id", "refseq_transcript_id"
    )

    # Get features
    gene_features = get_genomic_features(ensembl_gene_id)
    five_prime_utr_stats = get_transcript_features(ensembl_transcript_id)
    transcript_features = get_all_orfs_features(ensembl_transcript_id)

    # Score and position data
    constraint = get_constraint_score(ensembl_gene_id)
    start_site = five_prime_utr_stats["start_site_pos"]

    # Clinical data
    clingen_entry = get_clingen_entry(hgnc)
    omim_id = get_omim_id(ensembl_gene_id)

    # smORF data
    smorfs = get_smorfs(ensembl_transcript_id)

    # Variants
    possible_variants = get_possible_variants(ensembl_transcript_id)
    all_possible_variants = find_all_high_impact_utr_variants(ensembl_transcript_id)
    few_possible_variants = get_first_variants(ensembl_transcript_id)

    # UTR regions
    utr_regions = [i for i in gene_features if i["type"] == "five_prime_UTR"]

    # Gnomad data
    gnomad_data, gnomad_variants_list, clinvar_variants_list = process_gnomad_data(
        get_gnomad_variants_in_utr_regions(utr_regions), ensembl_transcript_id
    )

    gnomad_utr_impact = get_utr_annotation_for_list_variants(
        gnomad_variants_list, possible_variants, start_site, buffer
    )
    clinvar_utr_impact = get_utr_annotation_for_list_variants(
        clinvar_variants_list, possible_variants, start_site, buffer
    )

    # URLs for external services
    impact_url = current_app.config["IMPACT_URL"]
    search_url = current_app.config["SEARCH_URL"]

    return render_template(
        "viewer.html",
        ensembl_transcript_id=ensembl_transcript_id,
        ensembl_gene_id=ensembl_gene_id,
        hgnc=hgnc,
        name=name,
        refseq_match=refseq_match,
        smorfs=smorfs,
        clingen_entry=clingen_entry,
        omim_id=omim_id,
        impact_url=impact_url,
        few_possible_variants=few_possible_variants,
        search_url=search_url,
        gnomad_data=gnomad_data,
        constraint=constraint,
        gene_features=gene_features,
        five_prime_utr_stats=five_prime_utr_stats,
        transcript_features=transcript_features,
        gnomad_utr_impact=gnomad_utr_impact,
        clinvar_utr_impact=clinvar_utr_impact,
        all_possible_variants=all_possible_variants,
    )

