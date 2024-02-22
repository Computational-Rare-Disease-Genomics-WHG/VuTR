# Filters smorfs to mane locations
# and annotations their relative transcript position
# Input is the smORF tsv file curated by the CRDG Novo
# and is available upon request

library("data.table")
library("magrittr")
library("stringi")
library("stringr")
library("optparse")
library("logger")


#' Returns the TIS -6,+4
#' @param pos
#' @param cdna_seq
#' @returns context : substring containing context TIS (pos-6,pos+4)
get_context <- function(cdna_seq, pos) {
    context <- substr(cdna_seq, pos - 6, pos + 4)
    return(if (nchar(context) == 11 & !is.na(context)) context else NA)
}
#' Get the kozak context
#' @param seq the transcript sequence
#' @param start_site the location of the a in the start codon
#' @return the 7 base kozak context, "-" if doesn't apply
get_kozak_context <- function(seq, start_site) {
    context <- substr(seq, start_site - 3, start_site + 3)
    return(if (nchar(context) == 7) context else "-")
}

#' Get the kozak strength
#' @param kozak_context, the kozak context of the sequence
#' @return the strength of the kozak as per rules ("Strong", "Moderate", "Weak")
get_kozak_consensus_strength <- function(kozak_context) {
    if (nchar(kozak_context) != 7) {
        return("None")
    }
    # Get the first and last bases that determine the kozak strength
    fb <- substr(kozak_context, 1, 1)
    lb <- substr(kozak_context, 7, 7)

    # Check cases for the Kozak Strength
    if ((fb == "a" || fb == "g") && lb == "g") {
        return("Strong")
    } else if ((fb == "a" || fb == "g") || lb == "g") {
        return("Moderate")
    } else {
        return("Weak")
    }
}


option_list <- list(
  make_option(
    c("-m", "--mane-file"),
    dest = "mane_file",
    type = "character",
    help = "Mane GFF file (as TSV)",
    metavar = "character"
  ),
  make_option(
    c("-r", "--rna-file"),
    dest = "rna_file",
    type = "character",
    help = "File with MANE RNA sequences",
    metavar = "character"
  ),

  make_option(
    c("-e", "--translational-efficiency"),
    dest = "te_file",
    type = "character",
    help = "File with translational efficiency",
    metavar = "character"
  ),

  make_option(
    c("-s", "--smorf"),
    type = "character",
    help = "Smorf file with the bed location of the smorfs",
    metavar = "character"
  ),
  make_option(
    c("-g", "--g2tcoord"),
    type = "character",
    help = "Genome 2 Transcript output file",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output-file"),
    dest = "output_file",
    type = "character",
    help = "Output file name (should end in .tsv)",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Parse command line options
mane_file <- opt$mane_file
smorf_file <- opt$smorf
genome_transcript_fn <- opt$g2tcoord
output_file_fn <- opt$output_file
te_file <- opt$te_file
rna_file <- opt$rna_file

log_info("Reading files")
smorfs <- fread(smorf_file, sep="\t")
mane <- fread(mane_file, sep="\t")
tmap <- fread(genome_transcript_fn, sep="\t")
te <- fread(te_file, sep="\t")
rna <- fread(rna_file, sep="\t")

# Renaming headers
names(smorfs) <- c(
  "chr",
  "start",
  "end",
  "smorf_id",
  "score",
  "strand",
  "thick_start",
  "thick_end",
  "item_rgb",
  "block_count",
  "block_sizes",
  "block_starts",
  "aa_seq",
  "start_codon",
  "smorf_names",
  "smorf_datasets",
  "dataset_count",
  "name",
  "cluster",
  "filtering",
  "confidence",
  "type",
  "alternate_types",
  "tx_id",
  "gene_name",
  "gene_id",
  "alternate_transcripts",
  "alternate_gene_ids",
  "alternate_gene_transcripts",
  "smorf_length",
  "isoform_count"
) # nolint

all_chrs <- paste0("chr", c(1:22, "X", "Y"))
# Create a unique "id" for each smorf
smorfs[, id := .I]
mane[, stable_transcript_id := substring(transcript_id, 1, 15)]
setkey(smorfs, id)

# Convert smorf 0-based to 1-based coordinates bedFile
# to match the mane coordinates
# which are 1-based
smorfs[, start :=  start + 1]


log_info("Filtering")

# Filter down to five prime UTRs
mane %<>% .[type == "exon"]
mane %<>% .[, .(transcript_id, stable_transcript_id)]
mane %<>% unique()

# Filter to necessary chromosomes
# Filter smorfs to those that are mapped to MANE coordinates
# and are uORFs or uoORFs
smorfs %<>% .[chr %in% all_chrs &
  tx_id %in% mane$stable_transcript_id & 
  (type == "uORF" | type == "uoORF")
]

# Merge with MANE to get transcript_id
setkey(mane, stable_transcript_id)
setkey(smorfs, tx_id)
smorfs <- mane[smorfs]

# Rename a few columns
setnames(smorfs, "start", "genome_start")
setnames(smorfs, "end", "genome_end")

annotated_output <- smorfs

log_info("Annotating transcript coordinates")

# Convert IDs to tpos using transcript id and tpos
setkey(tmap, transcript_id, gpos)
setkey(annotated_output, transcript_id, genome_start)

# Convert to transcript relative coordinates
annotated_output[, transcript_start :=
  tmap[.(transcript_id = annotated_output$transcript_id,
  gpos = annotated_output$genome_start), .(tpos)]
]
setkey(tmap, transcript_id, gpos)
setkey(annotated_output, transcript_id, genome_end)
annotated_output[, transcript_end :=
  tmap[.(transcript_id = annotated_output$transcript_id,
  gpos = annotated_output$genome_end), .(tpos)]
]

# Swap for reverse strand
annotated_output[strand == "-", `:=` (
  transcript_start = transcript_end,
  transcript_end = transcript_start
)]

# # Add translational efficiency
# log_info("Adding translational efficiency")
rna %<>% .[, .(transcript_id, seq)]
setkey(rna, transcript_id)
setkey(annotated_output, transcript_id)
# Merge cdna_seq with annotated_output
annotated_output <- rna[annotated_output]
annotated_output[start_codon == "ATG", context := get_context(seq, transcript_start), by = id]
annotated_output[start_codon == "ATG", kozak_context := get_kozak_context(seq, transcript_start), by = id]
annotated_output[start_codon == "ATG", kozak_consensus_strength := get_kozak_consensus_strength(kozak_context), by = id]

# # Add translational efficiency
te %<>% .[, .(context, efficiency)]
setkey(te, context)
setkey(annotated_output, context)
annotated_output <- te[annotated_output]

# # Remove cdna_seq
annotated_output[, seq := NULL]

# # For non-ATG start codons, set context and efficiency to NA
annotated_output[
  start_codon != "ATG",
  `:=` (
    context = NA,
    efficiency = NA
  )
]


log_info("Writing output")
fwrite(
  annotated_output,
  output_file_fn,
  sep = "\t"
)
