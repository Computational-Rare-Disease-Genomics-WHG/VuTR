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


option_list <- list(
  make_option(
    c("-m", "--mane_tsv"),
    type = "character",
    help = "Mane GFF file (as TSV)",
    metavar = "character"
  ),
  make_option(
    c("-p", "--output-no-prefix"),
    dest = "output_no_prefix",
    type = "character",
    help = "Output file with no prefix",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output_file"),
    dest = "output_file",
    type = "character",
    help = "Output file name (should end in .tsv)",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Get args from CLI
mane_tsv_file_path <- opt$mane_tsv
output_no_prefix_file_path <- opt$output_no_prefix
output_file_path <- opt$output_file


mane_features <- fread(mane_tsv_file_path, sep = "\t")

mane_features %<>% .[type == "five_prime_UTR"]
mane_features %<>% .[, .(seqid, start, end)]
names(mane_features) <- c("CHROM", "POS", "POS_TO")
mane_features[, CHROM := as.character(CHROM)]

# For use with gnomAD
fwrite(mane_features, output_no_prefix_file_path, sep = "\t", col.names = F)
mane_features[, CHROM := substr(CHROM, 4, nchar(CHROM))]

# For use with clinvar
fwrite(mane_features, output_file_path, sep = "\t", col.names = F)
