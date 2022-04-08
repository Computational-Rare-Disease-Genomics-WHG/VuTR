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
MANE_VERSION = 1.0
ENSEMBL_VERSION = 104
ASSEMBLY = 'GRCh38'
GNOMAD_VERSION = '2.1.1'

tbl_models = {
    'mane_summary': {
        'location': f'MANE/{MANE_VERSION}/MANE.{ASSEMBLY}.v{MANE_VERSION}.summary.txt.gz',  # pylint: disable=C0301 # noqa: E501
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
        'location': f'MANE/{MANE_VERSION}/MANE_transcripts_v{MANE_VERSION}.tsv',  # pylint: disable=C0301 # noqa: E501
        'separator': '\t',
        'col_mappings': {
            'five_prime_utr_length': 'five_prime_utr_length',
            'three_prime_utr_length': 'three_prime_utr_length',
            'num_five_prime_utr_exons': 'num_five_prime_utr_exons',
            'start_site_pos': 'start_site_pos',
            'cds_start': 'cds_start',
            'cds_end': 'cds_end',
            'strand': 'strand',
            'gene_symbol': 'hgnc_symbol',
            'seq': 'seq',
            'cds_length': 'cds_length',
            'ensembl_transcript_id': 'ensembl_transcript_id',
            # Need to add exon features somehow
        },
        'remove_ensembl_id_version_numbers': True,
        'ensembl_ids': ['ensembl_transcript_id'],
        'dtype': {
            'ensembl_transcript_id': VARCHAR(length=30),
            'five_prime_utr_length': Integer(),
            'three_prime_utr_length': Integer(),
            'num_five_prime_utr_exons': Integer(),
            'strand': VARCHAR(length=30),
            'seq': VARCHAR(length=100000),
            'start_site_pos': Integer(),
            'cds_start': Integer(),
            'cds_end': Integer(),
            'cds_length': Integer(),
        },
    },
    'orf_features': {
        'location': f'ORFS_Features_{MANE_VERSION}.tsv',  # pylint: disable=C0301  # noqa: E501
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
        'location': f'translational_efficiency.txt',
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
    'loeuf_constraint': {
        'location': f'GNOMAD/gnomad.v{GNOMAD_VERSION}.lof_metrics.by_transcript.txt',  # pylint: disable=C0301  # noqa: E501
        'separator': '\t',
        'col_mappings': {
            'gene': 'hgnc_symbol',
            'transcript': 'ensembl_transcript_id',
            'gene_id': 'ensembl_gene_id',
            'oe_lof_upper': 'loeuf',
        },
        'remove_ensembl_id_version_numbers': False,
        'ensembl_ids': None,
        'dtype': {
            'hgnc_symbol': VARCHAR(length=30),
            'ensembl_gene_id': VARCHAR(length=30),
            'ensembl_transcript_id': VARCHAR(length=30),
            'loeuf': Float(),
        },
    },
    'mane_genomic_features': {
        'location': f'MANE/{MANE_VERSION}/MANE.{ASSEMBLY}.v{MANE_VERSION}.select_ensembl_genomic.tsv',  # pylint: disable=C0301  # noqa: E501
        'separator': '\t',
        'col_mappings': {
            'seqid': 'chr',
            'source': 'source',
            'type': 'type',
            'start': 'start',
            'end': 'end',
            'score': 'score',
            'strand': 'strand',
            'phase': 'phase',
            'ID': 'ID',
            'gene_id': 'ensembl_gene_id',
            'gene_type': 'gene_type',
            'gene_name': 'hgnc_symbol',
            'Parent': 'parent',
            'transcript_id': 'ensembl_transcript_id',
            'transcript_type': 'transcript_type',
            'transcript_name': 'transcript_name',
            'tag': 'tag',
            'protein_id': 'ensembl_protein_id',
            'Dbxref': 'Dbxref',
            'exon_number': 'exon_number',
            'exon_id': 'exon_id',
        },
        'remove_ensembl_id_version_numbers': True,
        'ensembl_ids': [
            'ensembl_transcript_id',
            'ensembl_gene_id',
            'ensembl_protein_id',
            'exon_id',
        ],  # pylint: disable=C0301  # noqa: E501
        'dtype': None,  # To be finalized once we have things sorted out
    },
    'genome_to_transcript_coordinates': {
        'location': f'UTR_Genome_Transcript_Coordinates.tsv',  # pylint: disable=C0301  # noqa: E501
        'separator': '\t',
        'col_mappings': {
            'seqid': 'chr',
            'ensembl_transcript_id': 'ensembl_transcript_id',
            'strand': 'strand',
            'exon_number': 'exon_number',
            'gpos': 'genomic_pos',
            'tpos': 'transcript_pos',
        },
        'remove_ensembl_id_version_numbers': False,
        'ensembl_ids': None,
        'dtype': None,  # To be finalized once we have things sorted out
    },
    'orf_locations': {
        'location': f'UORFS_Genomic_Positions.tsv',
        'separator': '\t',
        'col_mappings': {
            'transcript_id': 'ensembl_transcript_id',
            'orf_id': 'orf_id',
            'strand': 'strand',
            'start': 'start',
            'end': 'end',
        },
        'remove_ensembl_id_version_numbers': True,
        'ensembl_ids': ['ensembl_transcript_id', 'orf_id'],
        'dtype': None,
    },
    'clingen': {
        'location': f'CLINGEN/ClinGen_gene_curation_list_{ASSEMBLY}.tsv',
        'separator': '\t',
        'col_mappings': {
            '#Gene Symbol': 'hgnc_symbol',
            'Gene ID': 'gene_id',
            'cytoBand': 'cytoband',
            'Genomic Location': 'genomic_location',
            'Haploinsufficiency Score': 'haplo_score',
            'Haploinsufficiency Description': 'haplo_description',
            'Haploinsufficiency PMID1': 'haplo_pmid1',
            'Haploinsufficiency PMID2': 'haplo_pmid2',
            'Haploinsufficiency PMID3': 'haplo_pmid3',
            'Haploinsufficiency PMID4': 'haplo_pmid4',
            'Haploinsufficiency PMID5': 'haplo_pmid5',
            'Haploinsufficiency PMID6': 'haplo_pmid6',
            'Triplosensitivity Score': 'triplo_score',
            'Triplosensitivity Description': 'triplo_description',
            'Triplosensitivity PMID1': 'triplo_pmid1',
            'Triplosensitivity PMID2': 'triplo_pmid2',
            'Triplosensitivity PMID3': 'triplo_pmid3',
            'Triplosensitivity PMID4': 'triplo_pmid4',
            'Triplosensitivity PMID5': 'triplo_pmid5',
            'Triplosensitivity PMID6': 'triplo_pmid6',
            'Date Last Evaluated': 'date_last_evaluated',
            'Loss phenotype OMIM ID': 'loss_phenotype_omim_id',
            'Triplosensitive phenotype OMIM ID': 'triplo_phenotype_omim_id',
        },
        'remove_ensembl_id_version_numbers': False,
        'ensembl_ids': None,
        'dtype': None,  # To be finalized once we have things sorted out
    },
}
