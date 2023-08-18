#!/usr/bin/env Rscript

# Processes the CADD conservation scores file to only include MANE transcripts
# and removes duplicate rows

# Usage example:
# Rscript conservation_scores.R \
# -i /path/to/conservation_scores.tsv \
#  -m /path/to/mane_summary.tsv \
#  -g /path/to/g2t.tsv \
#  -o /path/to/output.tsv \

library(data.table)
library(magrittr)
library(optparse)

# Define the options
option_list <- list(
    make_option(c("-i", "--input"), type = "character",
                help = "Path to UTR cons scores file", metavar = "character"),
    make_option(c("-g", "--g2tfile"), type = "character",
                help = "Genome To Transcript coordinate file",
                metavar = "character"),
    make_option(c("-o", "--output"),
    type = "character", default = "utr_cons.txt",
                help = "Output file path", metavar = "character")
)

# Parse the options
args <- parse_args(OptionParser(option_list = option_list))
input_file <- args$input
g2tfile <- args$g2tfile
output_file <- args$output


# Reading in the files
dt <- fread(input_file)
g2t <- fread(g2tfile)

# Rename headers
setnames(dt, c("chr", "gpos", "ref", "cons_score",
    "ensembl_transcript_id", "phastcons",
    "phylop", "gerp_rs", "gerp_rs_pval", "gerp_n",
    "gerp_s", "raw_cadd", "phred_cadd"))

setnames(g2t,
    c("chr", "transcript_id", "strand", "exon_number",
    "gpos", "tpos", "ensembl_transcript_id")
)

dt[, chr:= paste0("chr", chr)]

# Filter to relevant cols
dt %<>% .[, .(
    chr, gpos, phastcons, phylop,
    gerp_rs, gerp_s, raw_cadd, phred_cadd
)]

g2t %<>% .[, .(
    chr, gpos, tpos, ensembl_transcript_id
)]


# Filter to unique sites
dt <- unique(dt)

# Set the keys
setkey(dt, chr, gpos)
setkey(g2t, chr, gpos)

# Merge the data.table
dt <- g2t[dt]

# Write to file
fwrite(dt, output_file, sep = "\t")
