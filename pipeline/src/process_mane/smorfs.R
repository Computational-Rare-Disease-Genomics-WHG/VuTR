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

option_list <- list(
  make_option(
    c("-m", "--mane-file"),
    dest = "mane_file",
    type = "character",
    help = "Mane GFF file (as TSV)",
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


log_info("Reading files")
smorfs <- fread(smorf_file, sep="\t")
mane <- fread(mane_file, sep="\t")
tmap <- fread(genome_transcript_fn, sep="\t")

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

# Filter to necessary chromosomes
smorfs %<>% .[chr %in% all_chrs]

# Create a unique "id" for each smorf
smorfs[, id := .I]
setkey(smorfs, id)

log_info("Filtering")
# Filter down to five prime UTRs
mane %<>% .[type == "five_prime_UTR"]
mane %<>% .[tag == "MANE_Select"]
mane %<>% .[, .(seqid, start, end, strand, transcript_id)]
setkey(mane, start, end)

# Convert smorf 0-based to 1-based coordinates bedFile
# to match the mane coordinates
# which are 1-based

smorfs[, start :=  start + 1]

# Find all smorfs that lie within the mane coordinates
output <- smorfs[, {
  mane[
    seqid == .SD$chr &
    strand == .SD$strand &
    start <= .SD$start  &
    .SD$end <= end][, .(transcript_id)]
}, by = id]


setkey(output, id)

# Rename a few columns
setnames(smorfs, "start", "genome_start")
setnames(smorfs, "end", "genome_end")

annotated_output <- smorfs[output]

# Convert IDs to tpos using transcript id and tpos
setkey(tmap, transcript_id, gpos)
setkey(annotated_output, transcript_id, genome_start)

# Convert to transcript relative coordinates
annotated_output[, transcript_start :=
  tmap[.(transcript_id = annotated_output$transcript_id,
  gpos = annotated_output$genome_start), .(tpos)]
]

annotated_output[, transcript_end :=
  tmap[.(transcript_id = annotated_output$transcript_id,
  gpos = annotated_output$genome_end), .(tpos)]
]

# Swap for reverse strand
annotated_output[strand == "-", `:=` (
  transcript_start = transcript_end,
  transcript_end = transcript_start
)]

log_info("Writing output")
fwrite(
  annotated_output,
  output_file_fn,
  sep = "\t"
)
