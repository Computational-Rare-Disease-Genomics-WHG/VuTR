# transcripts.R
# E. D'Souza

# Read and the MANE transcripts and their annotations and turn it into a
# a readable csv file

library("magrittr")
library("data.table")
library("seqinr")
library("stringr")

# Read through MANE
ensembl_rna <-
    read.fasta(
            "../../data/pipeline/MANE/0.93/MANE.GRCh38.v0.93.select_ensembl_rna.fna.gz",
        as.string = T
    )

# Look at structure
ensembl_rna[1] %>% str

# functions from https://davetang.org/muse/2013/05/09/using-the-r-seqinr-package/
rna_names <- getName(ensembl_rna) %>% unlist
rna_annotations <- getAnnot(ensembl_rna) %>% unlist
rna_seqs <-
    sapply(rna_names, function (i)
        ensembl_rna[[i]][1]) %>% as.vector

# Create a dt with these sequences
mane_rna_dt <-
    data.table(
        ensembl_transcript_id = rna_names,
        annotations = rna_annotations,
        seq = rna_seqs
    )
setkey(mane_rna_dt, ensembl_transcript_id)

# A function that take the fasta header and change it to a named list
get_field_val <- function(ann) {
    # Convert into list separated by spaces
    p <- str_split(ann, " ") %>% unlist

    # Annotate the second field as type:cdna
    p[2] <- sprintf("type:%s", p[2])

    # paste the description from first occurance till the end of the string
    description_pos <- grep("description", p)

    # Concatenate the last couple of fields
    p <- c(p[1:description_pos - 1], paste(p[description_pos:length(p)], collapse = " "))

    # Remove first ID as it is not necessary
    p <- p[2:length(p)]

    # Split the annotation based into a fixed width matrix with
    # first_column being the field name and the second one the value
    key_val <- str_split_fixed(p, ":", 2)

    # Convert matrix into named list
    transcript_annotation <- as.vector(key_val[, 2])
    transcript_annotation <-
        as.list(setNames(as.vector(key_val[, 2]), as.vector(key_val[, 1])))
    return(transcript_annotation)
}

# Find the name of the new columns
new_col_names <- names(get_field_val(mane_rna_dt[1]$annotations))
# Add annotation fields into the data.table
mane_rna_dt[, (new_col_names) := (get_field_val(annotations)), by = ensembl_transcript_id]

# parse the chromosome column
mane_rna_dt[, (c("build", "chr", "start", "end", "strand")) := data.table(str_split_fixed(chromosome, ":", 5))]
mane_rna_dt[, chromosome := NULL]

# Save relevant features in file
fwrite(mane_rna_dt, "../../data/pipeline/MANE/0.93/MANE_transcripts_v0.93.tsv", sep="\t")
