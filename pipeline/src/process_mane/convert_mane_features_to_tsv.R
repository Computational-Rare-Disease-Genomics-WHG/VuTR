# convert_mane_features_to_tsv.R
# E. D'Souza

# Converts the gff features to CSV
# for easy reading

library("data.table")
library("rtracklayer")
library("magrittr")
library("optparse")

option_list <- list(
    make_option(c("-m", "--mane_gff"),
        type = "character",
        help = "dataset file name (should be .gff)", 
        metavar = "character"
    ),
    make_option(c("-o", "--output_file"),
        type = "character",
        help = "Output file name (should end in .tsv)", 
        metavar = "character"
    ),
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Parse command line options
mane_file <- opt$mane_gff
output_file <- opt$output_file

# Read file
mane <- mane_file %>%
    readGFF() %>%
    as.data.table()

# Write data.table
output_file %>%
    fwrite(mane, ., sep = "\t")
