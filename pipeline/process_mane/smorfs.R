
library("data.table")
library("magrittr")
library("stringi")
library("stringr")
library("optparse")


smorfs <- fread("SMORFS/all_final_orfCDS.txt")
names(smorfs) <- c('chr', 'source', 'orf_type', 
    'start', 'end', 'V6','strand', 'V8', 'iORF_id') # nolint
smorfs[, chr := paste0("chr", chr)]
smorfs[chr == 'chr24', chr := 'chrX' ]
smorfs[chr == 'chr25', chr := 'chrY' ]
smorfs[, id := .I]
setkey(smorfs, id)


mane <- fread("MANE/1.0/MANE.GRCh38.v1.0.ensembl_genomic.tsv")
mane %<>% .[type == 'five_prime_UTR']
mane %<>% .[tag == 'MANE_Select']
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

annotated_output <- smorfs[, .(source, orf_type, iORF_id,
    genome_start = start, genome_end = end, id)][output]


# Convert IDs to tpos using transcript id and tpos
tmap <- fread("UTR_Genome_Transcript_Coordinates.tsv")


setkey(tmap, transcript_id, gpos)
setkey(annotated_output, transcript_id, genome_start)

# Convert to transcript relative coordinates
annotated_output[,transcript_start :=
    tmap[.(transcript_id=annotated_output$transcript_id, 
        gpos=annotated_output$genome_start), .(tpos)]]

annotated_output[,transcript_end :=
    tmap[.(transcript_id=annotated_output$transcript_id, 
        gpos=annotated_output$genome_end), .(tpos)]]


fwrite(annotated_output,
    "SMORFS/smorfs_locations.tsv",
    sep="\t")
