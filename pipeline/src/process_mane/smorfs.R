
library("data.table")
library("magrittr")
library("stringi")
library("stringr")
library("optparse")

option_list <- list(
    make_option(c("-m", "--mane_tsv"),
        type = "character",
        help = "Mane GFF file (as TSV)",
        metavar = "character"
    ),
    make_option(c("-s", "--smorf"),
        type = "character",
        help = "Smorf file with the bed location of the smorfs",
        metavar = "character"
    ),

    make_option(c("-g", "--g2tcoord"),
        type = "character",
        help = "Genome 2 Transcript output file",
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
mane_file <- opt$mane_tsv
smorf_file <- opt$smorf
genome_transcript_fn <- opt$g2tcoord
output_file_fn <- opt$output_file


smorfs <- fread(smorf_file)
names(smorfs) <- c("chr", "source", "orf_type",
    "start", "end", "V6","strand", "V8", "iORF_id") # nolint
smorfs[, chr := paste0("chr", chr)]
smorfs[chr == "chr24", chr := "chrX"]
smorfs[chr == "chr25", chr := "chrY"]
smorfs[, id := .I]
setkey(smorfs, id)


mane <- fread(mane_file)
mane %<>% .[type == "five_prime_UTR"]
mane %<>% .[tag == "MANE_Select"]
mane %<>% .[, .(seqid, start, end, strand, transcript_id)]
setkey(mane, start, end)


# Find all smorfs that lie within
# MANE exons
output <- smorfs[, {
    mane[seqid == .SD$chr &
    strand == .SD$strand &
    start <= .SD$start  &
    .SD$end <= end
    ][, .(transcript_id)]
}, by = id]

setkey(output, id)

annotated_output <- smorfs[, .(source, orf_type, iORF_id, strand,
    genome_start = start, genome_end = end, id)][output]


# Convert IDs to tpos using transcript id and tpos
tmap <- fread(genome_transcript_fn)


setkey(tmap, transcript_id, gpos)
setkey(annotated_output, transcript_id, genome_start)

# Convert to transcript relative coordinates
annotated_output[,transcript_start :=
    tmap[.(transcript_id=annotated_output$transcript_id, 
        gpos=annotated_output$genome_start), .(tpos)]]

annotated_output[,transcript_end :=
    tmap[.(transcript_id =annotated_output$transcript_id,
        gpos = annotated_output$genome_end), .(tpos)]]

# Swap for reverse strand
annotated_output[
    strand == "-", `:=` (
        transcript_start = transcript_end,
        transcript_end = transcript_start
    )
]

fwrite(annotated_output,
    output_file_fn,
    sep = "\t"
)
