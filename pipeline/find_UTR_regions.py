# find_UTR_regions.py 
# E.D'Souza

from gtfparse import read_gtf

# returns GTF with essential columns such as "feature", "seqname", "start", "end"
# alongside the names of any optional keys which appeared in the attribute column
df = read_gtf("../db/MANE/MANE.GRCh38.v0.93.select_ensembl_genomic.gtf.gz")


# Looks through the MANE file used and then 

# PATH to change using snakemake config