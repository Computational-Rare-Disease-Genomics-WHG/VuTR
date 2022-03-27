"""Deals with the possible variants"""
# pylint: skip-file
# flake8: noqa
# TODO: Clean utils
import pandas as pd
import json
from pathlib import Path


def parse_values(val, start_site, buffer_length):
    """
    Converts the val into an int
    """
    if val == "NA":
        return start_site + buffer_length
    else:
        return int(val)


def get_utr_annotation_for_list_variants(
    list_variants, possible_variants_dict, start_site, buffer_length
):

    high_impact_utr_variants = list(
        set(
            [
                v['variant_id']
                for v in possible_variants_dict
                if v['variant_id'] in list_variants
            ]
        )
    )

    if len(high_impact_utr_variants) > 0:
        return [
            find_intervals_for_utr_consequence(
                var_id=v['variant_id'],
                conseq_type=v['five_prime_UTR_variant_consequence'],
                conseq_dict=v['five_prime_UTR_variant_annotation'],
                cdna_pos=v['cDNA_position'],
                start_site=start_site,
                buffer_length=buffer_length,
            )
            for v in possible_variants_dict
            if v['variant_id'] in high_impact_utr_variants
        ]
    return []


def parse_five_prime_UTR_variant_consequence(conseq_str):
    """
    Parses the consequence str into a keyed dictionary as per
    https://github.com/ImperialCardioGenetics/UTRannotator#the-detailed-annotation-for-each-consequence

    """
    return {
        annotation.split(':')[0]: annotation.split(':')[1]
        for annotation in conseq_str.split(',')
    }


def find_intervals_for_utr_consequence(
    var_id, conseq_type, conseq_dict, cdna_pos, start_site, buffer_length
):
    """
    Parses the output of UTR annotator to a dictionary of intervals [start, end] for the visualization
    """
    intervals = {}
    intervals['variant_id'] = var_id
    conseq_dict = parse_five_prime_UTR_variant_consequence(conseq_dict)
    print(conseq_dict)
    if conseq_type == 'uAUG_gained':
        # Done
        intervals['start'] = cdna_pos
        intervals['end'] = cdna_pos + parse_values(
            conseq_dict['uAUG_gained_DistanceToStop'], start_site, buffer_length
        )
        intervals['viz_type'] = 'New Feature'
        intervals['viz_color'] = 'main'
        intervals['type'] = 'uAUG_gained'
        intervals['kozak_strength'] = conseq_dict['uAUG_gained_KozakStrength']

    elif conseq_type == 'uAUG_lost':
        # Done
        intervals['start'] = int(conseq_dict['uAUG_lost_CapDistanceToStart'])
        intervals['end'] = start_site - parse_values(
            conseq_dict['uAUG_lost_DistanceToCDS'], start_site, buffer_length
        )
        intervals['viz_type'] = 'New Feature'
        intervals['viz_color'] = 'null'
        intervals['type'] = 'uAUG_lost'
        intervals['kozak_strength'] = conseq_dict['uAUG_lost_KozakStrength']

    elif conseq_type == 'uSTOP_lost':
        intervals['start'] = cdna_pos
        intervals['end'] = cdna_pos
        intervals['viz_type'] = 'New Feature'
        intervals['viz_color'] = 'main'
        intervals['type'] = 'uSTOP_lost'
        intervals['kozak_strength'] = conseq_dict['uSTOP_lost_KozakStrength']

    elif conseq_type == 'uSTOP_gained':
        # Get the cdna position of the start site
        intervals['start'] = start_site - int(
            conseq_dict['uSTOP_gained_ref_StartDistanceToCDS']
        )
        # cDNA position of the new stop gained
        intervals['end'] = start_site - parse_values(
            conseq_dict['uSTOP_gained_newSTOPDistanceToCDS'], start_site, buffer_length
        )
        intervals['viz_type'] = 'New Feature'
        intervals['viz_color'] = 'main'
        intervals['type'] = 'uSTOP_gained'
        intervals['kozak_strength'] = conseq_dict['uSTOP_gained_KozakStrength']

    # Once we have indels as well
    elif conseq_type == 'uFrameshift':
        pass
    intervals.update(conseq_dict)

    return intervals


def get_all_possible_variants(ensembl_transcript_id, cursor):
    cursor.execute(
        'SELECT annotations FROM variant_annotations WHERE ensembl_transcript_id =?',
        [ensembl_transcript_id],
    )
    rows = cursor.fetchall()
    variants = [json.loads(row[0]) for row in rows]
    return variants
