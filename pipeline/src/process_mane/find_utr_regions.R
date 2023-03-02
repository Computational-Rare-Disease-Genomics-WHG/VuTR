# find_utr_regions.R
# E. D'Souza

# Takes the MANE genomic features file and then
# creates a regions file for tabix to run with the -R
# flag



library("data.table")
library("magrittr")
library("stringr")
library("stringi")
library("optparse")
setwd("../../../")
option_list <- list(
    make_option(c("-m", "--mane_version"),
        type = "character", default = "1.0",
        help = "dataset file name", metavar = "character"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
mane_version <- opt$mane_version

mane_features <- "data/pipeline/MANE/%s/MANE.GRCh38.v%s.ensembl_genomic.tsv" %>%
    sprintf(., mane_version, mane_version) %>%
    fread(., sep = "\t")

mane_features %<>% .[type == "five_prime_UTR"]
mane_features %<>% .[, .(seqid, start, end)]
names(mane_features) <- c("CHROM", "POS", "POS_TO")
mane_features[, CHROM := as.character(CHROM)]


# For use with gnomAD
fwrite(mane_features,
    "data/pipeline/UTR_regions_with_chr_prefix.tsv",
    sep = "\t",
    col.names = F
)


mane_features[, CHROM := substr(CHROM, 4, nchar(CHROM))]

# For use with clinvar
fwrite(mane_features,
    "data/pipeline/UTR_regions.tsv",
    sep = "\t",
    col.names = F
)
