# find_uORFs.R

# finds all of the uORFS and their relavant features
# (frame / kozak / transcript start and end points)
# for a given transcript

library("data.table")
library("magrittr")
library("stringi")
library("stringr")

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


transcripts <- fread(
    "../../data/pipeline/MANE/0.93/MANE_transcripts_v0.93.tsv",
    sep = "\t"
)


# Apply for each transcript id
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

        if (nrow(orfs) > 1) {
            # Find the ORF sequence and stop sites
            orfs[, (c(
                "orf_seq",
                "orf_stop_codon"
            )) := find_stop(seq, orf_start_codon), by = orf_start_codon]

            # Characterize which are uORFS and which are oORFs
            orfs[orf_stop_codon > start_site_pos, orf_type := "oORF"]
            orfs[orf_stop_codon < start_site_pos, orf_type := "uORF"]

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


# Combine the data tables together
orfs <- transcripts[
    ,
    rbindlist(setNames(
        orfs,
        ensembl_transcript_id
    ), idcol = "ensembl_transcript_id")
]

# Write to file
fwrite(orfs,
    "../../data/pipeline/ORFS_Features_0.93.tsv",
    sep = "\t"
)