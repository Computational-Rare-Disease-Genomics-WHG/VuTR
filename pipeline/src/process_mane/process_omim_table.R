# Description: This script processes the omim table to only include gene entries and remove entries without ensembl gene ids
# Input: omim table
# Output: processed omim table

library(optparse)
library(data.table)
library(magrittr)

# define options
options_list <- list(
  make_option(
    c("-i", "--input"),
    dest = "input",
    type = "character",
    help = "Path to the omim downloaded file from OMIM"
  ),
  make_option(
    c("-o", "--output"),
    dest = "output",
    type = "character",
    help = "Path to the output processed file"
  )
)
# Pass cmd line args
parser <- OptionParser(option_list = options_list)
options <- parse_args(parser)
input_file_path <- options$input
output_file_path <- options$output

b <- fread(input_file_path)

names(b) <- c("omim_entry", "type", "entrez", "hgnc", "ensembl_gene_id")
b %<>% .[type == "gene"]
b %<>% .[!(is.na(ensembl_gene_id)|ensembl_gene_id == "")]



fwrite(b, output_file_path, sep = "\t")
