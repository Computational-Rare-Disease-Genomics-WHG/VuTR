"""
Module to get access to gnomAD data

Can test this in in gnomad
https://gnomad.broadinstitute.org/api
"""

# pylint: skip-file
# flake8: noqa
# TODO: Clean utils

import json
import requests  # pylint: disable=E0401


import pandas as pd

from pathlib import Path
from functools import reduce
from utr_utils.tools.utils import convert_uploaded_variation_to_variant_id

script_path = Path(__file__).parent


def gnomad_search_by_transcript_id(transcript_id):
    """
    For prototyping purposes
    Access gnomad through the GraphQL API
    @param transcript_id (str) : Ensembl transcript identifier
    (ensure this is the canonical MANE id)
    @returns response (dict) : The GraphQL object on response as a dictionary.
    """

    transcript_query = """
  query get_data($transcript_id: String!){
    transcript(transcript_id: $transcript_id, reference_genome: GRCh38) {
      transcript_id
      transcript_version
      strand
      start
      stop
      gene {
        gene_id
        name
        start
        stop
        chrom
        symbol
        hgnc_id
        omim_id
        mane_select_transcript{
          refseq_id
          ensembl_id

        }

      }
      clinvar_variants {
        transcript_id
        ref
        pos
        alt
        in_gnomad
        clinvar_variation_id
        gold_stars
        variant_id
        review_status
        hgvsc
        clinical_significance
        major_consequence
      }
      variants(dataset: gnomad_r3) {
        ref
        pos
        variant_id
        alt
        hgvsc
        genome {
          af
          an
          ac
        }
        transcript_consequence {
          is_mane_select
          major_consequence
          sift_prediction
          polyphen_prediction
          is_mane_select_version
        }
        variantId
      }
    }
  }
  """

    response = requests.post(
        'https://gnomad.broadinstitute.org/api',
        data=json.dumps(
            {
                'query': transcript_query,
                'variables': {'transcript_id': transcript_id},
            }
        ),
        headers={
            'Content-Type': 'application/json',
        },
    ).json()
    return response['data']


def gnomad_search_by_gene_id(hgnc):
    """
    For prototyping purposes

    Access gnomad through the GraphQL API

      Params :
        gene_id (str) : HGNC gene name

      Returns :
        response (dict) : The GraphQL object on response as a dictionary.

    """
    gene_query = """
  query get_data ($hgnc:String!)
    {
      gene(gene_symbol: $hgnc, reference_genome:GRCh38){
          gene_id
          name
          omim_id
        transcripts{
          transcript_version
          transcript_id

        }

        clinvar_variants {
          transcript_id
          ref
          pos
          alt
          in_gnomad
          clinvar_variation_id
          variant_id
          gold_stars
          review_status
          hgvsc
          clinical_significance
        }
        variants(dataset: gnomad_r3) {
          ref
          pos
          alt
          genome {
            af
            an
            ac
          }
          transcript_consequence {
            is_mane_select
            sift_prediction
            polyphen_prediction
            is_mane_select_version
          }
          variantId
        }

      }
    }
  """
    response = requests.post(
        'https://gnomad.broadinstitute.org/api',
        data=json.dumps(
            {
                'query': gene_query,
                'variables': {'hgnc': hgnc},
            }
        ),
        headers={
            'Content-Type': 'application/json',
        },
    ).json()

    return response['data']


def get_constraint_by_ensg(ensg):
    """
    For prototyping purposes

    Looks through the gnomad constraint
    from 2.1.1 to get constraint by gene id.

    @param ensg (str) : Ensembl gene id (stable)
    @returns gene_constraint_records (dict) : Constraint records for the gene (if found)
    """
    constraint = pd.read_csv(
        script_path
        / '../../data/pipeline/GNOMAD/gnomad.v2.1.1.lof_metrics.by_gene.txt',
        sep='\t',
    )
    # remove version number
    ensg = ensg[0:15]
    gene_constraint_records = constraint[constraint['gene_id'] == ensg].to_dict(
        'records'
    )[0]
    return gene_constraint_records


def get_gnomad_variants_in_utr_regions(utr_regions):
    searches = [
        gnomad_api_search_by_region(
            chrom=ur['chr'][3:], start=ur['start'], stop=ur['end']
        )['region']
        for ur in utr_regions
    ]
    data = {}
    # Append
    data['variants'] = sum([s['variants'] for s in searches], [])
    data['clinvar_variants'] = sum([s['clinvar_variants'] for s in searches], [])
    return data


def gnomad_api_search_by_region(chrom, start, stop):
    """
    For prototyping purposes
    """
    region_variant_query = """


query get_data ($chrom : String!,
								$start : Int!,
								$stop : Int!){
  region(chrom: $chrom, start:$start, stop:$stop, reference_genome:GRCh38){
  clinvar_variants {
        transcript_id
        ref
        pos
        alt
        in_gnomad
        clinvar_variation_id
        gold_stars
        variant_id
        review_status
        hgvsc
        clinical_significance
        major_consequence
      }
      variants(dataset: gnomad_r3) {
        ref
        pos
        alt
        hgvsc
        variant_id
        genome {
          af
          an
          ac
        }
        transcript_consequence {
          is_mane_select
          major_consequence
          sift_prediction
          polyphen_prediction
          is_mane_select_version
        }
      }
  }

  }
  """
    response = requests.post(
        'https://gnomad.broadinstitute.org/api',
        data=json.dumps(
            {
                'query': region_variant_query,
                'variables': {'start': start, 'stop': stop, 'chrom': chrom},
            }
        ),
        headers={
            'Content-Type': 'application/json',
        },
    ).json()

    return response['data']
