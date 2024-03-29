
library("data.table")
library("magrittr")
library("stringi")
library("stringr")
library("optparse")

# define options
options_list <- list(
  make_option(
    c("-m", "--mane-tsv"),
    dest = "mane_file",
    type = "character",
    help = "Path to the mane file"
  ),
  make_option(
    c("-o", "--output-file-path"),
    dest = "output_file_path",
    type = "character",
    help = "Path to the output intron file"
  ),
  make_option(
    c("-b", "--output-bed"),
    dest = "output_bed",
    type = "character",
    help = "Path to the output Bed file"
  )
)
parser <- OptionParser(option_list = options_list)
options <- parse_args(parser)
mane_file_path <- options$mane_file
output_file_path <- options$output_file_path
output_bed <- options$output_bed

mane_dt <- fread(mane_file_path)
mane <- mane_dt [
    type == "exon"
]
setorder(mane, transcript_id, exon_number)

flanking_size <- 20

# Filter to forward strand
flank_dt <- data.table(
    genome_start = c(mane$start - flanking_size, mane$end + 1),
    genome_end = c(mane$start - 1, mane$end + flanking_size),
    exon_id = rep(mane$exon_id, 2),
    exon_number = rep(mane$exon_number, 2),
    transcript_id = rep(mane$transcript_id),
    relative_position = rep(c("LHS", "RHS"), nrow(mane)),
    chr = rep(mane$seqid, 2),
    strand = rep(mane$strand, 2)
)

# Create intron feature_id
flank_dt[, intron_id := paste0(
    relative_position, "_", exon_number, "_", exon_id
)]

# Filter out the codons that are further
# than the transcription start site
flank_dt %<>%  .[
    !((exon_number == 1 & relative_position == "LHS" & strand == "+") | (
        exon_number == 1 & relative_position == "RHS" & strand == "-"
    ))
]

# Create a intron stichting flat table
fwrite(flank_dt, output_file_path, sep = "\t")

# Create bed files
fwrite(
    flank_dt[,
        .(
            chr,
            genome_start,
            genome_end,
            intron_id,
            ".",
            strand
        )
    ], output_bed, sep = "\t", col.names = FALSE
)
