"""
The declarative schema for the databases

tbl_models = {
    'table_name' : {
        'location' : # location
        'separator' : # separator between cols
        'col_mappings' :  #renaming columns, ensure all cols
                        # that wish to be in the db are mentioned here
                        # even if the name isn't changed
            {
                'old_name' : 'new_name', # new name
                'old_name' : 'old_name' # retain name

            }
        'remove_ensembl_id_version_numbers': True, #Whether there are
        # ensembl identifier with trailing version numbers

        'ensembl_ids': [ # A list of ensembl_ids in this dataframe
                         # (USING THE NEW COL NAMES)
            'ensembl_transcript_id',
            'ensembl_gene_id',
            'ensembl_protein_id',
        ],

    }
}

"""
from sqlalchemy.types import (
    VARCHAR,
    Integer,
    Float,
)  # pylint: disable=E0401

# Replace with config/config.yml values
MANE_VERSION = 0.93
ENSEMBL_VERSION = 104
ASSEMBLY = 'GRCh38'
GNOMAD_VERSION = '2.1.1'

tbl_models = {
    'mane_summary': {
        'location': f'/MANE/{MANE_VERSION}/MANE.{ASSEMBLY}.v{MANE_VERSION}.summary.txt.gz',  # pylint: disable=C0301 # noqa: E501
        'separator': '\t',
        'col_mappings': {
            'Ensembl_nuc': 'ensembl_transcript_id',
            'Ensembl_Gene': 'ensembl_gene_id',
            'Ensembl_prot': 'ensembl_protein_id',
            '#NCBI_GeneID': 'ncbi_gene_id',
            'RefSeq_nuc': 'refseq_transcript_id',
            'RefSeq_prot': 'refseq_protein_id',
            'MANE_status': 'mane_status',
            'name': 'name',
            'symbol': 'hgnc_symbol',
            'HGNC_ID': 'hgnc_id',
        },
        'remove_ensembl_id_version_numbers': True,
        'ensembl_ids': [
            'ensembl_transcript_id',
            'ensembl_gene_id',
            'ensembl_protein_id',
        ],
        'dtype': {
            'ensembl_transcript_id': VARCHAR(length=30),
            'ensembl_gene_id': VARCHAR(length=30),
            'ensembl_protein_id': VARCHAR(length=30),
            'ncbi_gene_id': VARCHAR(length=30),
            'refseq_transcript_id': VARCHAR(length=30),
            'refseq_protein_id': VARCHAR(length=30),
            'mane_status': VARCHAR(length=30),
            'name': VARCHAR(length=30),
            'hgnc_symbol': VARCHAR(length=30),
            'hgnc_id': VARCHAR(length=30),
        },
    },
    'mane_transcript_features': {
        'location': f'/MANE_transcript_features_v{MANE_VERSION}.tsv',  # pylint: disable=C0301  # noqa: E501
        'separator': '\t',
        'col_mappings': {
            'five_prime_utr_length': 'five_prime_utr_length',
            'three_prime_utr_length': 'three_prime_utr_length',
            'num_five_prime_utr_exons': 'num_five_prime_utr_exons',
            'start_site_pos': 'start_site_pos',
            'cds_start': 'cds_start',
            'cds_end': 'cds_end',
            'cds_length': 'cds_length',
            'ensembl_transcript_id': 'ensembl_transcript_id',
        },
        'remove_ensembl_id_version_numbers': True,
        'ensembl_ids': ['ensembl_transcript_id'],
        'dtype': {
            'ensembl_transcript_id': VARCHAR(length=30),
            'five_prime_utr_length': Integer(),
            'three_prime_utr_length': Integer(),
            'num_five_prime_utr_exons': Integer(),
            'start_site_pos': Integer(),
            'cds_start': Integer(),
            'cds_end': Integer(),
            'cds_length': Integer(),
        },
    },
    'orf_features': {
        'location': f'/ORFS_Features_{MANE_VERSION}.tsv',  # pylint: disable=C0301  # noqa: E501
        'separator': '\t',
        'col_mappings': {
            'ensembl_transcript_id': 'ensembl_transcript_id',
            'orf_start_codon': 'orf_start_codon',
            'orf_stop_codon': 'orf_stop_codon',
            'orf_seq': 'orf_seq',
            'orf_type': 'orf_type',
            'frame': 'frame',
            'kozak_context': 'kozak_context',
            'kozak_consensus_strength': 'kozak_consensus_strength',
            'orf_id': 'orf_id',
        },
        'remove_ensembl_id_version_numbers': True,
        'ensembl_ids': ['ensembl_transcript_id', 'orf_id'],
        'dtype': {
            'ensembl_transcript_id': VARCHAR(length=30),
            'orf_start_codon': Integer(),
            'orf_stop_codon': Integer(),
            'orf_seq': VARCHAR(length=10000),
            'orf_type': VARCHAR(length=30),
            'frame': VARCHAR(length=30),
            'kozak_context': VARCHAR(length=30),
            'kozak_consensus_strength': VARCHAR(length=30),
            'orf_id': VARCHAR(length=30),
        },
    },
    'translational_efficiencies': {
        'location': f'/translational_efficiency.txt',
        'separator': '\t',
        'col_mappings': {
            'context': 'context',
            'efficiency': 'efficiency',
            'lower_bound': 'lower_bound',
            'upper_bound': 'upper_bound',
        },
        'remove_ensembl_id_version_numbers': False,
        'ensembl_ids': None,
        'dtype': {
            'context': VARCHAR(length=30),
            'efficiency': Integer(),
            'lower_bound': Integer(),
            'upper_bound': Integer(),
        },
    },
    'constraint': {
        'location': f'/GNOMAD/gnomad.v{GNOMAD_VERSION}.lof_metrics.by_transcript.txt',  # pylint: disable=C0301  # noqa: E501
        'separator': '\t',
        'col_mappings': {
            'gene': 'hgnc_name',
            'transcript': 'ensembl_transcript_id',
            'gene_id': 'ensembl_gene_id',
            'oe_lof_upper': 'loeuf',
        },
        'remove_ensembl_id_version_numbers': False,
        'ensembl_ids': None,
        'dtype': {
            'hgnc_name': VARCHAR(length=30),
            'ensembl_gene_id': VARCHAR(length=30),
            'ensembl_transcript_id': VARCHAR(length=30),
            'loeuf': Float(),
        },
    },
}

# MANE Genomic features needs a little bit of thought on how to replace this
#     'mane_genomic_features': {
#        'location':  f'/MANE/{MANE_VERSION}/MANE.{ASSEMBLY}.v{MANE_VERSION}.select_ensembl_genomic.tsv',  # pylint: disable=C0301  # noqa: E501
#        'separator': '\t',
#        'col_mappings': {
#
#        },
#        'remove_ensembl_id_version_numbers': True,
#        'ensembl_ids': ['ensembl_transcript_id', 'ensembl_gene_id', 'ensembl_protein_id'],  # pylint: disable=C0301  # noqa: E501
#
#        'dtype': {
#            'chr': VARCHAR(length=30),
#            'start': BigInteger(),
#            'stop': BigInteger(),
#            'strand': VARCHAR(length=30),
#            'type': VARCHAR(length=30),
#        },
#    },
