# filter_to_UTR.R
# E. D'Souza

# Takes the MANE genomic features file and then
# creates a regions file for tabix to run with the -R
# flag

library(data.table)
library(magrittr)
library(rtracklayer)
library(stringr)
library(stringi)

setwd("../..")

mane_features <- readGFF("data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.select_ensembl_genomic.gff.gz") %>% as.data.table

mane_features  %<>% . [type == "five_prime_UTR"]
mane_features %<>% .[, .(seqid, start, end)]
names(mane_features) <- c("CHROM", "POS", "POS_TO")
mane_features[, CHROM:= as.character(CHROM)]
mane_features[, CHROM := substr(CHROM, 4, nchar(CHROM))]

fwrite(mane_features, "data/pipeline/UTR_regions.tsv", sep="\t", col.names=F)
