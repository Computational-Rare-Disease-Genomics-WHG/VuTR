"""
Utility functions for dealing with cdna sequences
"""

import pandas as pd
import re


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
        'ncbi_gene_id': '#NCBI_GeneID'
    }
    # reading mane
    summary = pd.read_csv(
        '../../data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.summary.txt.gz', sep='\t')
    try:
        transformed_id = summary[summary[colmappings[from_type]
                                         ].str.contains(id)][colmappings[to_type]].values[0]
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

    if (start_site-comparison_site) % 3 == 0:
        return "In-frame"
    else:
        return "Out-of-frame"


def find_uorfs_in_transcript(
        seq,
        start_site):
    """
    @param seq : Transcript sequence (as cdna)
    @param start_site : The position of thst start site
    @return List[dicts] of uorfs and their relative cdna atg pos/frame/context/strength as keys
    """
    if start_site-1 == 0:
        return None

    five_prime_utr_seq = seq[0:start_site]  # Check points

    # find all downstream atgs
    atg_pos = [i for i in range(len(five_prime_utr_seq))
               if five_prime_utr_seq.startswith("atg", i)]

    # filter to those that have in-frame stop
    start_end = [(i, i+find_stop(five_prime_utr_seq, i)+3)
                 for i in atg_pos if not find_stop(five_prime_utr_seq, i) is None]

    uorfs = [{"atg_pos": i[0],
              'stop_codon': five_prime_utr_seq[i[1]-3:i[1]],
              "stop_codon_pos": i[1]-3,
              'start': i[0],
              'end': i[1],
              "seq": five_prime_utr_seq[i[0]:i[1]],
              'length': len(five_prime_utr_seq[i[0]:i[1]]),
              "frame": get_relative_frame(start_site, i[0]),
              "kozak_context": get_kozak_context(seq, i[0]),
              "kozak_strength": get_kozak_strength(get_kozak_context(seq, i[0]))}
             for i in start_end]

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

    stop_codon_start = [p.start() for p in re.finditer(
        stop_grep_str, seq) if get_relative_frame(i, p.start()) == "In-frame"]
    if len(stop_codon_start) > 0:
        return stop_codon_start[0]
    else:
        return None


def find_oorf_in_transcript(seq, start_site):
    """
    """
    pass


def get_kozak_context(seq, pos):
    """
    Gets the 3 base context around pos inside 
    @param seq (str): A cdna sequence 
    @param pos (int) : position of the sequence
    """
    if pos-3 > 0 and pos+4 < len(seq):
        return seq[pos-3:pos+4]
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
    pass
