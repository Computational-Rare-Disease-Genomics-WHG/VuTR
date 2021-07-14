# gene_classes.py
# J. Martin

class Gene: #Note, add functionality to get 5'UTR Sequence

    def __init__(self, tempvars):
        self.transcript_id = tempvars["Temp_Transcript_ID"]
        self.chrom = tempvars["Temp_Chr"]
        self.gene_start_coord = tempvars["Temp_Gene_Start_Coord"]
        self.gene_end_coord = tempvars["Temp_Gene_End_Coord"]
        self.gene_strand = tempvars["Temp_Gene_Strand"]
        self.gene_id = tempvars["Temp_Gene_ID"]
        self.gene_symbol = tempvars["Temp_Gene_Symbol"]
        self.hgnc_symbol = tempvars["Temp_HGNC_Symbol"]
        self.ncbi_symbol = tempvars["Temp_NCBI_Symbol"]
        self.gene_description = tempvars["Temp_Gene_Description"]
        self.sequence = tempvars["Temp_cDNA_Sequence"]