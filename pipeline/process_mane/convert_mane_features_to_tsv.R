# convert_mane_features_to_tsv.R
# E. D'Souza

# Converts the gff features to CSV
# for easy reading

library("data.table")
library("rtracklayer")
library("magrittr")

setwd("../..")
mane <- "data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.select_ensembl_genomic.gff.gz" %>% # nolint
    readGFF() %>%
    as.data.table()
"data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.select_ensembl_genomic.tsv" %>%
    fwrite(mane, ., sep = "\t")