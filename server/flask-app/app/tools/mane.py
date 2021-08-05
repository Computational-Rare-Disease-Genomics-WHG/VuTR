
from Bio import SeqIO
from gtfparse import read_gtf
import gffutils
 
mane_gff_db = None


def create_mane_gff_database(db_name):
    fn = gffutils.example_filename('../../../db/MANE/MANE.GRCh38.v0.93.select_ensembl_genomic.gff.gz')
    mane_gff_db = gffutils.create_db(fn, dbfn='mane_gff_db.db', force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
    return mane_gff_db

def connect_to_gff_database(db_name):
    pass



def load_genomic_sequence (transcript_id):  
    pass    

def get_sequence_features(transcript_id):
    df = read_gtf("../../../db/MANE/MANE.GRCh38.v0.93.select_ensembl_genomic.gtf.gz")

    pass

def load_exons(transcript_id):
    pass