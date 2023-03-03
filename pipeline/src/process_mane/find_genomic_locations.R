# find_genomic_locations.R
# E.D'Souza
# This looks up the transcript relative coordinates in ORFS
# and converts them to genomic locations


library("data.table")
library("magrittr")
library("optparse")


option_list <- list(
    make_option(c("-f", "--orf-features"),
        dest = "orf_features",
        type = "character",
        help = "Orf features file (produced by find_orfs.R)",
        metavar = "character"
    ),
    make_option(c("-m", "--mane-file"),
        dest = "mane_file",
        type = "character",
        help = "Mane file produced by convert_mane_features_to_tsv.R",
        metavar = "character"
    ),
    make_option(c("-o", "--output-file"),
        dest = "output_file",
        type = "character",
        help = "Output file name",
        metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
orf_features_file_path <- opt$orf_features
mane_file_path <- opt$mane_file
output_file_path <- opt$output_file


uorfs <- fread(orf_features_file_path)
mane_gff <- fread(mane_file_path) # nolint

# Filter to exons only
mane_gff %<>% .[type == "exon"]

# set order by the transcript id
setorder(mane_gff, "transcript_id", "exon_number")

# Create the width of the exon
mane_gff[, exon_length := end - start + 1]

# Create the relative transcript start and end widths
mane_gff[, transcript_end := cumsum(exon_length), by = transcript_id]
mane_gff[, transcript_start := cumsum(data.table::shift(exon_length,
    fill = 1
)), by = transcript_id]

# Get the start position by looking up mane_gff
setkey(mane_gff, transcript_id, transcript_start, transcript_end)

uorfs_genomic_positions <- data.table()

for (i in seq(nrow(uorfs))) {
    # Find values for this iteration
    orf_id <- uorfs[i]$orf_id
    this_transcript_id <- uorfs[i]$ensembl_transcript_id
    orf_start <- uorfs[i]$orf_start_codon
    orf_end <- uorfs[i]$orf_stop_codon

    # Filter the gff entries and strand corresponding to this transcript
    transcript_gff <- mane_gff[transcript_id == this_transcript_id]
    strand <- transcript_gff[1]$strand

    if (strand == "+") {
        # find the genomic positions for the ATG
        # and the start of all succeeding exons

        genomic_starts <- c(
            transcript_gff[(transcript_start <= orf_start) &
                (transcript_end >= orf_start)]$start +

                (orf_start - transcript_gff[(transcript_start <= orf_start) &
                    (transcript_end >= orf_start)]$transcript_start), # nolint

            # Find the genomic start positions
            # of all of the exons until the exon where the stop is located
            transcript_gff[transcript_start >= orf_start &
                transcript_end <= orf_end]$start
        )

        # find the genomic location of the end of all of the exons
        genomic_ends <- c(
            transcript_gff[(transcript_start >= orf_start) &
                (transcript_end <= orf_end)]$end,

            # Find the genomic position of the TGA/TAA/TAG
            transcript_gff[(transcript_start <= orf_end) &
                (transcript_end >= orf_end)]$start +
                (orf_end - transcript_gff[(transcript_start <= orf_end) &
                    (transcript_end >= orf_end)]$transcript_start)
        )
    } else {
        # For the reverse strand

        # Find the genomic position of the stop codon
        genomic_starts <- c(
            transcript_gff[(transcript_start <= orf_end) &
                (transcript_end >= orf_end)]$start +
                (orf_end - transcript_gff[(transcript_start <= orf_end) &
                    (transcript_end >= orf_end)]$transcript_end), # nolint
            transcript_gff[transcript_end <= orf_end &
                transcript_start >= orf_start]$start
        )

        # Find the genomic position of the ATG start codon
        genomic_ends <- c(
            transcript_gff[transcript_end <= orf_end &
                transcript_start >= orf_start]$end,
            transcript_gff[(transcript_start <= orf_start) &
                (transcript_end >= orf_start)]$start +
                (orf_start - transcript_gff[(transcript_start <= orf_start) &
                    (transcript_end >= orf_start)]$transcript_end)
        )
    }
    # Create a dataframe to row bind this to
    n <- length(genomic_ends)
    uorfs_genomic_positions <- rbind(
        uorfs_genomic_positions,
        data.table(
            orf_id = rep(orf_id, n),
            strand = rep(strand, n),
            transcript_id = rep(this_transcript_id, n),
            start = genomic_starts,
            end = genomic_ends
        )
    )
}

fwrite(uorfs_genomic_positions, output_file_path, sep = "\t")
