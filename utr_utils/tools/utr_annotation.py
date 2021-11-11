"""Deals with the possible variants"""

import pandas as pd
import numpy as np
import re
from pathlib import Path

# Define script path
script_path = Path(__file__).parent


def load_possible_variants_df(ensembl_transcript_id):
    # Currentl
    possible_df = pd.read_csv(script_path /
                              "../../data/pipeline/vep_data/output/UTR_variants_vep_all_possible_GRCh38_0.93_chr5_parsed.txt", sep="\t")
    possible_df = possible_df.set_index('variant_id')
    return possible_df[ensembl_transcript_id == possible_df.Feature]


def parse_five_prime_UTR_variant_consequence(conseq_str):
    """
    Parses the consequence str into a keyed dictionary as per 
    https://github.com/ImperialCardioGenetics/UTRannotator#the-detailed-annotation-for-each-consequence

    """
    return {annotation.split(':')[0]: annotation.split(':')[1]
            for annotation in conseq_str.split(',')}


def parse_values(val, start_site, buffer_length):
    """
    Converts the val into an int
    """
    if val == "NA":
        return start_site + buffer_length
    else:
        return int(val)


def get_utr_annotation_for_list_variants(list_intervals, possible_variants_df, start_site, buffer_length):

    # filter list_intervals to those in possible_variants_df

    high_impact_utr_variants = list(set([
        var for var in list_intervals if var in possible_variants_df.index.to_list()]))
    if len(high_impact_utr_variants) > 0:
        return list(possible_variants_df.loc[high_impact_utr_variants]
                    .drop_duplicates()
                    .apply(
            lambda var: find_intervals_for_utr_consequence(
                var_id=var.name,
                conseq_type=var.five_prime_UTR_variant_consequence,
                conseq_str=var.five_prime_UTR_variant_annotation,
                cdna_pos=var.cDNA_position,
                start_site=start_site,
                buffer_length=buffer_length
            ),
            axis=1))
    return []


def find_intervals_for_utr_consequence(
        var_id,
        conseq_type,
        conseq_str,
        cdna_pos,
        start_site,
        buffer_length):
    """
    Parses the output of UTR annotator to a dictionary of intervals [start, end] for the visualization
    """
    conseq_dict = parse_five_prime_UTR_variant_consequence(conseq_str)

    intervals = {}
    intervals['variant_id'] = var_id
    if conseq_type == 'uAUG_gained':
        # Done
        intervals['start'] = cdna_pos
        intervals['end'] = cdna_pos+parse_values(
            conseq_dict['uAUG_gained_DistanceToStop'],
            start_site,
            buffer_length)
        intervals['viz_type'] = 'New Feature'
        intervals['viz_color'] = 'main'
        intervals['type'] = 'uAUG_gained'
        intervals['kozak_strength'] = conseq_dict['uAUG_gained_KozakStrength']

    elif conseq_type == 'uAUG_lost':
        # Done
        intervals['start'] = int(conseq_dict['uAUG_lost_CapDistanceToStart'])
        intervals['end'] = start_site-parse_values(
            conseq_dict['uAUG_lost_DistanceToCDS'],
            start_site,
            buffer_length)
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
        intervals['start'] = start_site - \
            int(conseq_dict['uSTOP_gained_ref_StartDistanceToCDS'])
        # cDNA position of the new stop gained
        intervals['end'] = start_site - parse_values(
            conseq_dict['uSTOP_gained_newSTOPDistanceToCDS'],
            start_site,
            buffer_length)
        intervals['viz_type'] = 'New Feature'
        intervals['viz_color'] = 'main'
        intervals['type'] = 'uSTOP_gained'
        intervals['kozak_strength'] = conseq_dict['uSTOP_gained_KozakStrength']

    # Once we have indels as well
    elif conseq_type == 'uFrameshift':
        pass

    return intervals
