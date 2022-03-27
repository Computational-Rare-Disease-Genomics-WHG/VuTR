# create_gpos_lookup.R
# E. D'Souza

# Generates a table mapping each transcript relative coordinate
# to its genomic coordinate in MANE (GRCh38) for easy swapping
# between the two entities

library(data.table)
library(magrittr)
library(rtracklayer)
library(optparser)

setwd("../../")
option_list <- list(
    make_option(c("-m", "--mane_version"),
        type = "character", default = "1.0",
        help = "dataset file name", metavar = "character"
    ),
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
mane_version <- opt$mane_version


# Read genomic feature file
genomic_mane <- "data/pipeline/MANE/%s/MANE.GRCh38.v%s.select_ensembl_genomic.tsv" %>% # nolint
    sprintf(., mane_version) %>%
    fread()
# Filter to exons
genomic_mane %<>% .[type == "five_prime_UTR"]

# Choose relevant cols
genomic_mane %<>% .[, .(seqid, start, end, strand, exon_number, transcript_id)]

# Rename cols and order by exon num
setnames(genomic_mane, "start", "gstart")
setnames(genomic_mane, "end", "gend")
setorderv(genomic_mane, cols = c("transcript_id", "exon_number"))


# Add relative transcript coordinates
genomic_mane[, exon_width := gend - gstart + 1]
genomic_mane[, tstart := cumsum(
    data.table::shift(exon_width, fill = 1)
), by = transcript_id]
genomic_mane[, tend := cumsum(
    exon_width
), by = transcript_id]

# Create lookup up table by reading along each exon sequentially
gt_lookup_dt <- genomic_mane %>%
    .[, .(
        gpos = if (strand == "+") seq(gstart, gend) else seq(gend, gstart),
        tpos = seq(tstart, tend)
    ), by = .(
        seqid,
        transcript_id,
        strand,
        exon_number
    )]

# Setkey for transcript and pos
setkey(gt_lookup_dt, transcript_id, tpos)

gt_lookup_dt[, ensembl_transcript_id := substr(transcript_id, 1, 15)]

# Write to file
fwrite(gt_lookup_dt,
    "data/pipeline/UTR_Genome_Transcript_Coordinates.tsv",
    sep = "\t"
)
