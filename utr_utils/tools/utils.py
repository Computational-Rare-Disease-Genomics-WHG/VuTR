# pylint: skip-file
# flake8: noqa
# TODO: Clean utils
"""Utility functions for dealing with cdna sequences"""

import pandas as pd
import numpy as np
import re
import json
from pathlib import Path

script_path = Path(__file__).parent


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


def read_summary_mane():
    summary = pd.read_csv(
        script_path / '../../data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.summary.txt.gz',
        sep='\t',
    )
    return summary


def convert_uploaded_variation_to_variant_id(uploaded_variation):
    """
    Replaces the uploaded variation in VEP to a gnomad-esq variant id
    """
    return uploaded_variation.replace('_', '-').replace('/', '-')


def find_transcript_location(glookup_table, variant_pos, ensembl_transcript_id):
    """
    Filters the glookup_table to find the genome position's transcript cdna location
    """
    return int(
        glookup_table[
            (variant_pos == glookup_table['gpos'])
            & (ensembl_transcript_id == glookup_table['ensembl_transcript_id'])
        ]['tpos'].values[0]
    )


def add_tloc_to_dict(variant_item, glookup_table, ensembl_transcript_id):
    """
    Wrapper function that adds a tpos field to a given dictionary element
    """
    variant_item['tpos'] = find_transcript_location(
        glookup_table, variant_item['pos'], ensembl_transcript_id
    )
    return variant_item


def get_lookup_df(ensembl_transcript_id):

    glookup_table = pd.read_csv(
        script_path / "../../data/pipeline/UTR_Genome_Transcript_Coordinates.tsv",
        sep="\t",
    )
    return glookup_table[ensembl_transcript_id == glookup_table.ensembl_transcript_id]


def read_mane_genomic_features(ensembl_gene_id):
    """
    Reads mane for a specific gene_id from the genomic feature file

    @param ensembl_gene_id (str) : A stable ensembl
     gene identifier (ensure that this is in MANE) e.g. ENSG00000081189
    @returns gene_data (DataFrame) : MANE filtered to that region
    """
    # TODO : This needs to be updated to Pipeline
    mane = pd.read_csv(
        script_path
        / '../../data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.select_ensembl_genomic.tsv',
        sep='\t',
    )
    mane['ensembl_stable_gene_id'] = mane['gene_id'].apply(lambda x: str(x)[0:14])
    gene_data = mane[mane['gene_id'] == ensembl_gene_id]
    return gene_data


def read_mane_transcript(ensembl_transcript_id):
    """
    Reads through the MANE transcript set to get the transcript sequences

    @param ensembl_transcript_id (str) : Full transcript, id (including version number)
    @returns transcript_df
    """

    transcript_df = pd.read_csv(
        script_path / "../../data/pipeline/MANE/0.93/MANE_transcripts_v0.93.tsv",
        sep="\t",
    )
    transcript_df = transcript_df[
        transcript_df["ensembl_transcript_id"].str.contains(ensembl_transcript_id)
    ]

    return transcript_df


def read_oorf_features(ensembl_transcript_id):

    orf_df = pd.read_csv(
        script_path / "../../data/pipeline/ORFS_Features_0.93.tsv", sep='\t'
    )

    return orf_df[orf_df.ensembl_transcript_id.str.contains(ensembl_transcript_id)]


def convert_betweeen_identifiers(id, from_type, to_type):
    """
    Convert between different identifiers

    @param id (str) : The id to transform
    @from_type (str) : One of the following 'ensembl_transcript', 'ensembl_protein' etc..
    @to_type (str) : The desired identifier 'ensembl_transcript', 'refseq_protein' etc..
    @returns transformed_id (str) : The identifier in the to_type
    """
    colmappings = {
        'ensembl_transcript': 'Ensembl_nuc',
        'ensembl_protein': 'Ensembl_prot',
        'ensembl_gene': 'Ensembl_Gene',
        'refseq_protein': 'RefSeq_prot',
        'refseq_mrna': 'RefSeq_nuc',
        'hgnc_symbol': 'symbol',
        'hgnc_id': 'HGNC_ID',
        'ncbi_gene_id': '#NCBI_GeneID',
        'name': 'name',
        'mane_status': 'MANE_status',
    }
    # reading mane
    summary = pd.read_csv(
        script_path / '../../data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.summary.txt.gz',
        sep='\t',
    )
    try:
        transformed_id = summary[summary[colmappings[from_type]].str.contains(id)][
            colmappings[to_type]
        ].values[0]
        return transformed_id

    except Exception as e:
        print(e)


def get_relative_frame(start_site, comparison_site):
    """
    Relative frame of i compared to start site
    @param start_site (int) : the "canonical" start_site we want to compare to
    @param comparison_site (int) : the position of the atg we want to compare it against.
    @returns (str) : Frame of i compared to start_site
    """

    if (start_site - (comparison_site + 1)) % 3 == 0:
        return "In-frame"
    else:
        return "Out-of-frame"


def find_genomic_interval(ensembl_transcript_id, transcript_start, transcript_end):
    """
    Finds the genomic interval
    """
    ensembl_gene_id = convert_betweeen_identifiers(
        ensembl_transcript_id, "ensembl_transcript", "ensembl_gene"
    )

    mane_features = read_mane_genomic_features(ensembl_gene_id)

    # get strand
    strand = mane_features[mane_features['type'] == 'gene']['strand'].item()

    # Swap if start is less than end, find genomic location of the interval
    if strand == '+':
        start_genome = convert_transcript_coordinates_to_genomic(
            mane_features, transcript_start
        )
        end_genome = convert_transcript_coordinates_to_genomic(
            mane_features, transcript_end
        )
    else:
        start_genome = convert_transcript_coordinates_to_genomic(
            mane_features, transcript_end
        )
        end_genome = convert_transcript_coordinates_to_genomic(
            mane_features, transcript_start
        )

    # Now check if there are introns between the two position
    exons = mane_features[mane_features['type'] == 'exon']
    overlapping_exons = exons[
        ((exons['start'] <= start_genome) & (start_genome <= exons['end']))
        | ((exons['start'] <= end_genome) & (end_genome <= exons['end']))
    ]

    # If there is an intron between the interval
    if overlapping_exons.shape[0] > 1:

        # Sort exons by exon number
        overlapping_exons = overlapping_exons.sort_values(
            by=['exon_number'], ascending=True
        )

        # Find the first start-end segment (end of the first exon)
        first_feature = [
            {"start": start_genome, "end": overlapping_exons['end'].values[0]}
        ]

        # If there a more than one intron, add them in the middle
        if overlapping_exons.shape[0] > 3:
            middle_features = [
                overlapping_exons.loc[1:-1, ["start", "end"]].to_dict("records")
            ]
            first_feature = first_feature + middle_features

        last_feature = [
            {
                "start": overlapping_exons['start'].values[
                    overlapping_exons.shape[0] - 1
                ],
                "end": end_genome,
            }
        ]
        # Add anything in the middle

        return first_feature + last_feature

    return [{"start": start_genome, "end": end_genome}]


def convert_transcript_coordinates_to_genomic(mane_features, pos):
    """
    Converts the transcript position to genomic coordinates based
    on the mane-features of the given pos (corresponding to the gene of interest)

    @param ensembl_transcript_id (str)
    @param pos (int)
    @returns dictionary
    """

    # get strand
    strand = mane_features[mane_features['type'] == 'gene']['strand'].item()

    # Filter to exons
    exons = mane_features[mane_features["type"] == "exon"]

    # Sort based on exon_number
    exons = exons.sort_values(by=['exon_number'], ascending=True)

    exons = exons[["start", "end", "exon_number"]]

    # Calculate the exon to map to transcript coordinates
    exons["width"] = exons["end"] - exons["start"] + 1
    exons["t_end"] = exons.width.cumsum()

    exon_starts = np.roll(exons["t_end"] + 1, 1)
    exon_starts[0] = 1
    exons["t_start"] = exon_starts

    # Found genomic exon
    found = exons.loc[(exons["t_start"] <= pos) & (pos <= exons["t_end"])]
    ref_exon_position = found["start"].values[0]

    delta = (
        (pos - found["t_start"].values[0])
        if strand == "+"
        else (found["t_end"].values[0] - pos)
    )
    genomic_pos = ref_exon_position + delta
    return genomic_pos


def find_uorfs_in_transcript(seq, start_site, ensembl_transcript_id):
    """
    @param seq : Transcript sequence (as cdna)
    @param start_site : The position of thst start site
    @return List[dicts] of uorfs and their relative cdna atg pos/frame/context/strength as keys
    """
    if start_site - 1 == 0:
        return None

    five_prime_utr_seq = seq[0:start_site]

    # find all downstream atgs
    atg_pos = [
        i
        for i in range(len(five_prime_utr_seq))
        if five_prime_utr_seq.startswith("atg", i)
    ]

    # filter to those that have in-frame stop
    # so we have a start_end tuple
    # e.g. [(start_of_uORF in transcript, end_of_uORF in transcript)...]
    start_end = [
        (i, i + find_stop(five_prime_utr_seq, i) + 3)
        for i in atg_pos
        if not find_stop(five_prime_utr_seq, i) is None
    ]

    # Create a dictionary of the element
    uorfs = [
        {
            "atg_pos": i[0],
            'stop_codon': five_prime_utr_seq[i[1] - 3 : i[1]],
            "stop_codon_pos": i[1] - 3,
            'start': i[0],
            'end': i[1],
            'genome_features': find_genomic_interval(ensembl_transcript_id, i[0], i[1]),
            "utr_seq": five_prime_utr_seq[i[0] : i[1]],
            'length': len(five_prime_utr_seq[i[0] : i[1]]),
            "frame": get_relative_frame(start_site, i[0]),
            "kozak_context": get_kozak_context(seq, i[0]),
            "kozak_strength": get_kozak_strength(get_kozak_context(seq, i[0])),
        }
        for i in start_end
    ]

    return uorfs


def find_stop(seq, i):
    """
    Finds the position of an in-frame stop codon
    @param seq : Transcript sequence search space
    @param i : position to start looking for
    @returns stop_pos
    """
    stop_codons = ['tga', 'taa', 'tag']
    stop_grep_str = "|".join(stop_codons)

    # Get downstream of seq only
    seq = seq[i:]

    stop_codon_start = [
        p.start()
        for p in re.finditer(stop_grep_str, seq)
        if get_relative_frame(i, p.start()) == "In-frame"
    ]
    if len(stop_codon_start) > 0:
        return stop_codon_start[0]
    else:
        return None


def find_oorf_in_transcript(seq, start_site):
    """"""
    pass


def get_kozak_context(seq, pos):
    """
    Gets the 3 base context around pos inside
    @param seq (str): A cdna sequence
    @param pos (int) : position of the sequence
    """
    if pos - 3 > 0 and pos + 4 < len(seq):
        return seq[pos - 3 : pos + 4]
    else:
        return None


def get_kozak_strength(context):
    """
    Returns the strength of the kozak
    @param context (str) : the +/-3 base context
    @returns strength (str) : The strength of the context (returns 'Not valid context' if not)
    """

    # Check lengths
    if not len(context) == 7:
        return 'Not valid context'

    fb = context[0]
    lb = context[6]

    if (fb == 'a' or fb == 'g') and lb == 'g':
        return 'Strong'
    elif (fb == 'a' or fb == 'g') and lb == 'g':
        return 'Moderate'
    else:
        return 'Weak'
