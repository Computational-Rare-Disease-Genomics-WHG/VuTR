# convert_mane_features_to_tsv.R
# E. D'Souza

# Converts the gff features to CSV
# for easy reading

library("data.table")
library("rtracklayer")
library("magrittr")
library("optparse")

option_list <- list(
    make_option(c("-m", "--mane_version"),
        type = "character", default = "1.0",
        help = "dataset file name", metavar = "character"
    ),
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
mane_version <- opt$mane_version

setwd("../..")
mane <- "data/pipeline/MANE/%s/MANE.GRCh38.v%s.select_ensembl_genomic.gff.gz" %>% # nolint
    sprintf(., mane_version) %>%
    readGFF() %>%
    as.data.table()
"data/pipeline/MANE/%s/MANE.GRCh38.v%s.select_ensembl_genomic.tsv" %>%
    sprintf(., mane_version) %>%
    fwrite(mane, ., sep = "\t")
