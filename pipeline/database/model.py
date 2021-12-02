"""The declarative schema for the databases"""
from sqlalchemy.types import BigInteger, Integer, String, Float  # pylint: disable=E0401


MANE_VERSION = 0.93
ENSEMBL_VERSION = 103

tbl_models = {
    'mane_summary': {
        'location': '',
        'remove_version_numbers': True,
        'dtype': {
            'ensembl_transcript_id': String,
            'ensembl_gene_id': String,
            'mane_status': String,
            'refseq_match': String,
            'hgnc_name': String,
            'hgnc_id': String,
            'name': String,
        },
    },
    'mane_genomic_features': {
        'location': '',
        'remove_version_numbers': True,
        'dtype': {
            'start': BigInteger,
            'stop': BigInteger,
            'chr': String,
            'strand': String,
            'type': String,
        },
    },
    'mane_transcript_features': {
        'location': '',
        'remove_version_numbers': True,
        'dtype': {
            'ensembl_transcript_id': String,
            'five_prime_utr_length': Integer,
            'three_prime_utr_length': Integer,
            'num_five_prime_utr_exons': Integer,
            'start_site_pos': Integer,
            'cds_start': Integer,
            'cds_end': Integer,
            'cds_length': Integer,
        },
    },
    'orf_features': {
        'location': '',
        'remove_version_numbers': True,
        'dtype': {
            'ensembl_transcript_id': String,
            'orf_start_codon': Integer,
            'orf_stop_codon': Integer,
            'orf_seq': String,
            'orf_type': String,
            'frame': String,
            'kozak_context': String,
            'kozak_consensus_strength': String,
            'orf_id': String,
        },
    },
    'translational_efficiencies': {
        'location': '',
        'remove_version_numbers': False,
        'dtype': {
            'sequence': String,
            'efficiency': Integer,
            'lower_bound': Integer,
            'upper_bound': Integer,
        },
    },
    'constraint': {
        'location': '',
        'remove_version_numbers': False,
        'dtype': {
            'ensembl_gene_name': String,
            'ensembl_transcript_id': String,
            'loeuf': Float,
        },
    },
}
