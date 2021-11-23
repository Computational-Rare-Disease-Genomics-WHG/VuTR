# model.py

# Elston D'Souza
# The declarative Model / Schema for the databases


gnomad_variants_query = '''
CREATE TABLE IF NOT EXISTS gnomad_variants 
    (variant_id varchar, 
    pop_ac int,
    pop_af float,
    ref varchar,
    alt varchar
    major_consequence varchar)'''

clinvar_variants_query = '''
CREATE TABLE IF NOT EXISTS clinvar_variants 
    (variant_id varchar, 
    ref varchar, 
    alt varchar,
    allele_id int,
    review_status varchar,
    clinsig varchar)'''

contraint_query = '''
CREATE TABLE IF NOT EXISTS constraint 
    (ensembl_gene_name varchar,
    ensembl_transcript_id varchar,
    loeuf float)'''

possible_utr_variants_query = '''
CREATE TABLE IF NOT EXISTS possible_utr_variants 
    (variant_id varchar,
    five_prime_utr_consequence varchar,
    five_prime_utr_annotation data,
    intervals data)'''


mane_summary_query = '''
CREATE TABLE IF NOT EXISTS mane_summary (
    ensembl_transcript_id varchar, 
    ensembl_gene_id varchar, 
    mane_status varchar, 
    refseq_match varchar,
    hgnc_name varchar
)
'''

mane_features_query = '''
CREATE TABLE IF NOT EXISTS mane_genomic_features (
    start int, 
    stop int, 
    chr int, 
    strand varchar
    type varchar,
)
'''

mane_transcript_seqs_query = '''
CREATE TABLE IF NOT EXISTS mane_transcript_seqs (
    ensembl_transcript_id  varchar, 
    annotations varchar,
    seq varchar, 
    type varchar, 
    gene varchar,
    gene_biotype varchar, 
    transcript_biotype varchar, 
    gene_symbol varchar, 
    description varchar, 
    build varchar, 
    chr varchar, 
    start varchar, 
    end varchar, 
    strand varchar,
    utr_stats data, 
    uorfs data, 
    oorfs data,
)
'''


tbl_queries = {
    "": {
        'gnomad_variants': gnomad_variants_query,
        'clinvar_variants': clinvar_variants_query,
        'constraint': contraint_query,
        'possible_utr_variants': possible_utr_variants_query,
        'mane_summary': mane_summary_query,
        'mane_genomic_features': mane_features_query,
        'mane_transcript_features': mane_transcript_seqs_query,
    }
}
