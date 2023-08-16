#!/usr/bin/env Rscript

# Processes the CADD conservation scores file to only include MANE transcripts
# and removes duplicate rows

library(data.table)
library(optparse)

# Define the options
option_list <- list(
    make_option(c("-i", "--input"), type = "character",
                help = "Path to UTR cons scores file", metavar = "character"),
    make_option(c("-m", "--mane"), type = "character",
                help = "Path to MANE summary file", metavar = "character"),
    make_option(c("-o", "--output"),
    type = "character", default = "utr_cons.txt",
                help = "Output file path", metavar = "character")
)

# Parse the options
args <- parse_args(OptionParser(option_list=option_list))
input_file <- args$input
mane_file <- args$mane
output_file <- args$output

# Original code, modified to use the file paths from CLI
b <- fread(input_file)
b <- unique(b)
mane <- fread(mane_file)
mane[, transcript_id := substr(Ensembl_nuc, 1, 15)]
transcript_ids <- mane$transcript_id
b <- b[FeatureID %in% transcript_ids]
setnames(b, c("chrom", "pos", "ref", "cadd", "transcript_id", "phast_cons",
"phylop", "gerp_rs", "gerp_rs_pval", "gerp_n", "gerp_s"))
fwrite(b, output_file, sep = "\t")
