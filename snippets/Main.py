# Main.py 
# J. Martin

import gene_classes as gc
import gene_functions as gf
import operator
import sys

search = "MEF2C" 

file = "Data/data_files/MANE.GRCh38.v0.91.select_ensembl_rna.fna"
x = gf.read_MANE(file)

#Note here that the read_MANE function is used to generate all the gene summaries, but returns a list with all of these, rather than 
#the specific gene required - helpful to build a program that does this differently, so that the string that is searched for is returned

for gene in x:
    if gene.gene_symbol == search:
        print(gene.gene_id)
        print(gene.transcript_id)
        print(gene.chrom)
        print(gene.gene_start_coord)
        print(gene.gene_end_coord)
        print(gene.gene_strand)
        print(gene.gene_id)
        print(gene.gene_symbol)
        print(gene.hgnc_symbol)
        print(gene.ncbi_symbol)
        print(gene.gene_description)
        print(gene.sequence)
        seq = gene.sequence
        print(len(gene.sequence))


name = "Data/data_files/MANE.GRCh38.v0.92.select_ensembl_genomic.gff"

structure = gf.compile_gene_structure(name, search) #Note to self - in visualisation, show splice sites for segments using faint vertical line etc
print(structure)
# To find the co-ordinate of a nucleotide of given base
# for each dict in the gene_structure list
# if segment's index start <= nucleotide_index <=index end
# nucleotide coordinate = (nucleotide_index - index start (+/-1?)) + Beginning_co-ord

reading_frames = gf.generate_frames(seq)


# uORFs_unsorted = [] #Note, this will have to go a layer above when functionality expanded to incorporate all reading frames

unsorted_ORFs = []
frame0_ORFs = gf.generate_ORF_list(reading_frames, 0, seq, structure)
frame1_ORFs = gf.generate_ORF_list(reading_frames, 1, seq, structure)
frame2_ORFs = gf.generate_ORF_list(reading_frames, 2, seq, structure)
unsorted_ORFs.extend(frame0_ORFs)
unsorted_ORFs.extend(frame1_ORFs)
unsorted_ORFs.extend(frame2_ORFs)
ORFs = sorted(unsorted_ORFs, key=lambda k: k[-1])
uORFs = []
for ORF in ORFs:
    region_code = (ORF[-2])[0] + (ORF[-2])[1] + (ORF[-2])[2] + (ORF[-2])[3] + (ORF[-2])[4]
    if (region_code == "5'UTR"):
        uORFs.append(ORF)
# print(ORFs, sep = "\n")
print(*uORFs, sep = "\n")
# print(*uORFs, sep = "\n")
# print(*ORFs, sep = "\n")

# print("FRAME 0")
# print(*frame0_ORFs, sep = "\n") 
# print("FRAME 1")
# print(*frame1_ORFs, sep = "\n") 
# print("FRAME 2")
# print(*frame2_ORFs, sep = "\n") 

# def sort_uORFs():
#     # Must add functionality to decide how to create the overall picture - given that whether a start is translated or not will 
#     # alter the landscape of the uORFs 
#     uORFs = sorted()
#     return(uORFs)

# generate_uORF_list(1)
# generate_uORF_list(1)
# generate_uORF_list(2)

# print(seq)