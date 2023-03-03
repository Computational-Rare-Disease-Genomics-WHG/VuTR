

library(data.table)
library(magrittr)


# Pass cmd line args
parser <- OptionParser()

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
parser <- add_options(parser, options_list)

input_file_path <- options$input
output_file_path <- options$output

b <- fread(input_file_path)

names(b) <- c("omim_entry", "type", "entrez", "hgnc", "ensembl_gene_id")
b %<>% .[type == "gene"]
b %<>% .[!(is.na(ensembl_gene_id)|ensembl_gene_id == "")]



fwrite(b, output_file_path, sep = "\t")
