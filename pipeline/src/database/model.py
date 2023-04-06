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
        'dtype': None
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
        'dtype': None
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
            'efficiency': 'efficiency',
            'lower_bound': 'lower_bound',
            'upper_bound': 'upper_bound',
            'orf_start_codon_genome': 'orf_start_codon_genome',
            'orf_stop_codon_genome': 'orf_stop_codon_genome',
            'context': 'context',
        },
        'remove_ensembl_id_version_numbers': True,
        'ensembl_ids': ['ensembl_transcript_id', 'orf_id'],
        'dtype': None
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
        'dtype': None
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
        'dtype': None
    },
    'mane_genomic_features': {
        'location': f'MANE/{MANE_VERSION}/MANE.{ASSEMBLY}.v{MANE_VERSION}.ensembl_genomic.tsv',  # pylint: disable=C0301  # noqa: E501
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
        'dtype': None,  
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
        'dtype': None,  
    },
    'smorf_locations': {
        'location': f'SMORFS/smorfs_locations.tsv',  # pylint: disable=C0301  # noqa: E501
        'separator': '\t',
        'col_mappings': {
            'source': 'source',
            'orf_type': 'smorf_orf_type',
            'iORF_id': 'smorf_iorf_id',
            'genome_start': 'genome_start',
            'genome_end': 'genome_end',
            'transcript_id': 'ensembl_transcript_id',
            'transcript_start': 'transcript_start',
            'transcript_end': 'transcript_end',
        },
        'remove_ensembl_id_version_numbers': True,
        'ensembl_ids': ['ensembl_transcript_id'],
        'dtype': None, 
    },

    'smorf_features': {
        'location': f'SMORFS/final_orfs_filt_v5.txt',  # pylint: disable=C0301  # noqa: E501
        'separator': '\t',
        'col_mappings': {
            'ORF_id': 'smorf_orf_id',
            'iORF_id': 'smorf_iorf_id',
            'iORF_type': 'smorf_iorf_type',
            'gene_id': 'ensembl_gene_id',
            'iORF_pept': 'smorf_iORF_pept',
            'gene_name': 'gene_name',
            'gene_biotype': 'gene_biotype',
            'source': 'source',

            'reads_used': 'reads_used',
            'pct_reads_inframe': 'pct_reads_inframe',
            'pct_codons_inframe': 'pct_codons_inframe',
            'len': 'len',
            'strand': 'strand',
            'dropoff_score': 'dropoff_score',
            'starts': 'starts',
            'No._of_tools': 'num_tools',
            'Peptide_seq': 'peptide_seq',
        },
        'remove_ensembl_id_version_numbers': False,
        'ensembl_ids': None,
        'dtype': None
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

    'omim' : {
        'location' : 'processed_omim.txt',
        'separator' : '\t',
        'col_mappings' : {
            'omim_entry' : 'omim_entry',
            'type' : 'type',
            'entrez' : 'entrez',
            'hgnc' : 'hgnc_id',
            'ensembl_gene_id' : 'ensembl_gene_id'
        },
        'remove_ensembl_id_version_numbers' : False,
        'ensembl_ids' : None,
        'dtype' : None,
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
            'Date Last Evaluated': 'date_last_evaluated'},
        'remove_ensembl_id_version_numbers': False,
        'ensembl_ids': None,
        'dtype': None, 
    }
}