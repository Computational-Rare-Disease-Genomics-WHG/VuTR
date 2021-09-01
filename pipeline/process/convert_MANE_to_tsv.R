# convert_MANE_to_csv.R
# E. D'Souza

# Converts the gff features to CSV
# for easy reading

library(data.table)
library(rtracklayer)
library(magrittr)

setwd("../..")
mane <- readGFF("data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.select_ensembl_genomic.gff.gz") %>% as.data.table
fwrite(mane, "data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.select_ensembl_genomic.csv", sep="\t")
