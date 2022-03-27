"""Deals with the possible variants"""
# pylint: skip-file
# flake8: noqa
# TODO: Clean utils
import pandas as pd
import json
from pathlib import Path


def get_all_possible_variants(ensembl_transcript_id, cursor):
    cursor.execute(
        'SELECT annotations FROM variant_annotations WHERE ensembl_transcript_id =?',
        [ensembl_transcript_id],
    )
    rows = cursor.fetchall()
    variants = [json.loads(row[0]) for row in rows]
    return variants
