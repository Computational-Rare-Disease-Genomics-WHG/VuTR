b = fread("mim2gene.txt")
names(b) <- c("omim_entry", "type", "entrez", "hgnc", "ensembl_gene_id")
b %<>% .[type == "gene"]
b %<>% .[!(is.na(ensembl_gene_id)|ensembl_gene_id =="")]
fwrite(b, "processed_omim.txt", sep="\t")
