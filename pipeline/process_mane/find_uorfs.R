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


feature_file <- fread("data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.select_ensembl_genomic.tsv", sep = "\t") # nolint
transcripts <- fread("data/pipeline/MANE/0.93/MANE_transcripts_v0.93.tsv", sep = "\t") # nolint

feature_file[, width := end - start + 1]
transcript_feats <- feature_file[,
    .(
        five_prime_utr_length = sum(.SD[type == "five_prime_UTR"]$width),
        three_prime_utr_length = sum(.SD[type == "three_prime_UTR"]$width),
        num_five_prime_utr_exons = nrow(.SD[type == "five_prime_UTR"]),
        start_site_pos = sum(.SD[type == "five_prime_UTR"]$width) + 1,
        cds_start = nrow(.SD[type == "five_prime_UTR"]) + 1,
        cds_end = sum(.SD[type == "CDS"]$width) + nrow(.SD[type == "five_prime_UTR"]), # nolint
        cds_length = sum(.SD[type == "CDS"]$width)
    ),
    by = transcript_id
]
setkey(transcript_feats, transcript_id)
setkey(transcripts, ensembl_transcript_id)

b <- transcript_feats[transcripts]


transcripts %<>% .[1:10]

# calculate UTR stats
# such as UTR length / number of exons


# find the canonical start codon
transcripts[, .(
    if (start_site_pos > 1) {

        # find all uORFs, defined by
        # the prescence of a uAUG
        orfs <- data.table(str_locate_all(substr(
            seq,
            1, start_site_pos
        ), "atg")[[1]])[, .(start)]

        # find the stop codon
        colnames(ORFS) <- c("orf_start_codon")

        # Find the ORF sequence and stop sites
        orfs[, (c(
            "orf_seq",
            "orf_stop_codon"
        )) := find_stop(seq, orf_start_codon)]

        # Characterize which are uORFS and which are oORFs
        orfs[orf_stop_codon > start_site_pos, orf_type := "oORF"]
        orfs[orf_stop_codon > start_site_pos, orf_type := "uORF"]

        # find the frame
        orfs[, frame := categorize_frame(start_pos, orf_start_codon)]

        # find the kozak context
        orfs[, kozak_context := get_kozak_context(cdna_seq, start_site_pos)]
        orfs[, kozak_consensus_strength := get_kozak_consensus_strength(kozak_context)] # nolint

        # name each of the orfs
        orfs[, orf_id := paste(
            orf_type,
            .I,
            ensembl_transcript_id,
            sep = ":"
        ), by = ensembl_transcript_id]

        # Return the data.table
        .(.(orfs))
    }
), by = ensembl_transcript_id]