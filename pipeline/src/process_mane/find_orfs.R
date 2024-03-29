# find_uORFs.R

# finds all of the uORFS and their relavant features
# (frame / kozak / transcript start and end points)
# for a given transcript

library("data.table")
library("magrittr")
library("stringi")
library("stringr")
library("optparse")


#' Returns the TIS -6,+4
#' @param pos
#' @param cdna_seq
#' @returns context : substring containing context TIS (pos-6,pos+4)
get_context <- function(cdna_seq, pos) {
    context <- substr(cdna_seq, pos - 6, pos + 4)
    return(if (nchar(context) == 11 & !is.na(context)) context else "-")
}

#' Process the strings (to cdna) and
#' effsiciencies of the translational efficiency table
#' @param dt : The data.table of the read translation efficiency
#' @return data.table : the necessary columns and the efficiency per context
process_te_table <- function(dt) {
    dt[, tis_context := gsub("u", "t", tolower(sequence))] # nolint
    setkey(dt, tis_context) # nolint
    # Add quantiles
    eff_quantiles <- quantile(dt$efficiency, prob = seq(0, 1, 0.1))
    dt[, efficiency_decile := as.numeric(cut(efficiency, # nolint
        eff_quantiles,
        include.lowest = T,
        label = F
    ))]

    return(dt[, .(
        tis_context, # nolint
        efficiency, # nolint
        efficiency_decile # nolint
    )])
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


#' Reads through the sequence in triplets from the orf start site
#' and gives the sequence and position of the stop codon
#' @param seq the full transcript sequence
#' @param orf_start start_site
#' @return list of the orf_seq and the orf stop pos
find_stop <- function(seq, orf_start_codon) {
    # stop codon patterns
    stop_grep_pattern <- paste0(c("taa", "tag", "tga"), collapse = "|")
    # split the sequence into triplets looking for stop codons
    orf_codons <- strsplit(substr(
        seq,
        orf_start_codon, nchar(seq)
    ), "(?<=.{3})", perl = TRUE)[[1]]

    if (any(grepl(stop_grep_pattern, orf_codons))) {
        stop_codon_index <- grep(stop_grep_pattern, orf_codons)[1]
        # collapse codons until the stop
        orf_seq <- paste0(orf_codons[1:stop_codon_index], collapse = "")
        return(list(
            orf_seq = orf_seq,
            orf_stop_codon = orf_start_codon + 3 * stop_codon_index
        ))
    } else {
        return(list(orf_seq = "", orf_stop_codon = NA))
    }
}


categorize_frame <- function(start_pos, ref_pos) {
    if ((ref_pos - start_pos) %% 3 == 0) {
        return("Inframe")
    } else if ((ref_pos - start_pos) %% 3 == 1) {
        return("Out-of-Frame (1bp)")
    } else if ((ref_pos - start_pos) %% 3 == 2) {
        return("Out-of-Frame (2bp)")
    }
}

# define options
options_list <- list(
  make_option(
    c("-m", "--rna-file"),
    dest = "mane_rna_file",
    type = "character",
    help = "Path to the mane file"
  ),
  make_option(
    c("-g", "--g2tcoord"),
    dest = "g2tcoord",
    type = "character",
    help = "Genome 2 Transcript output file",
    metavar = "character"
  ),
  make_option(
    c("-t", "--te-file"),
    dest = "te_file",
    type = "character",
    help = "Path to the output RNA file"
  ),
  make_option(
    c("-o", "--output-orf-features"),
    dest = "output_orf_features",
    type = "character",
    help = "Path to the output RNA feature file"
  )
)


parser <- OptionParser(option_list = options_list)
options <- parse_args(parser)

mane_file_path <- options$mane_rna_file
te_file_path <- options$te_file
g2t_file_path <- options$g2tcoord
output_orf_features_file_path <- options$output_orf_features

transcripts <- fread(mane_file_path, sep = "\t")

# Apply for each transcript id
print("Finding orfs...")
transcripts[, orfs := {
    if (start_site_pos > 1) {

        # find all uORFs, defined by
        # the prescence of a uAUG
        orfs <- data.table(str_locate_all(substr(
            seq,
            1,
            start_site_pos
        ), "atg")[[1]])[, .(start)]

        # find the stop codon
        colnames(orfs) <- c("orf_start_codon")

        if (nrow(orfs) >= 1) {
            # Find the ORF sequence and stop sites
            orfs[, (c(
                "orf_seq",
                "orf_stop_codon"
            )) := find_stop(seq, orf_start_codon), by = orf_start_codon]

            # Characterize which are uORFS and which are oORFs
            orfs[orf_stop_codon > start_site_pos, orf_type := "oORF"]
            orfs[orf_stop_codon <= start_site_pos, orf_type := "uORF"]

            # find the frame
            orfs[, frame := categorize_frame(
                start_site_pos,
                orf_start_codon
            ), by = orf_start_codon]

            # find the kozak context
            orfs[,
                kozak_context := get_kozak_context(
                    seq,
                    orf_start_codon
                ),
                by = orf_start_codon
            ]

            # Add tis contexts
            orfs[,
                context := get_context(
                    seq,
                    orf_start_codon
                ),
                by = orf_start_codon
            ]

            orfs[,
                kozak_consensus_strength := get_kozak_consensus_strength(
                    kozak_context
                ),
                by = orf_start_codon
            ] # nolint

            # name each of the orfs
            orfs[, orf_id := paste(
                orf_type,
                .I,
                ensembl_transcript_id,
                sep = ":"
            ), by = orf_type]

            # Return the data.table
            .(.(orfs))
        }
    }
}, by = ensembl_transcript_id]
print("...Completed")

# Combine the data tables together
orfs <- transcripts[
    ,
    rbindlist(setNames(
        orfs,
        ensembl_transcript_id
    ), idcol = "ensembl_transcript_id")
]
print("Finding Genomic Coordinates...")

# Add genomic coordinates
genome_mapper <- fread(g2t_file_path, sep = "\t")
setkey(genome_mapper, transcript_id, tpos)
g <- genome_mapper[, .(transcript_id, gpos, tpos)]
setkey(g, transcript_id, tpos)
orfs[, orf_start_codon_genome := g[.(
    orfs$ensembl_transcript_id,
    orfs$orf_start_codon), gpos]]
orfs[, orf_stop_codon_genome := g[.(
    orfs$ensembl_transcript_id,
orfs$orf_stop_codon), gpos]]

print("Completed")

# Add translational efficiency
te <- fread(te_file_path)
setkey(te, context)
setkey(orfs, context)
orfs <- te[orfs]

setkey(orfs, orf_id)
# Write to file
fwrite(orfs,output_orf_features_file_path,sep = "\t")
