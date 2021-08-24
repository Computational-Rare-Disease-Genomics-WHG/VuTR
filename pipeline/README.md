# Pipeline 

This is a description of the upstream pipeline to annotate all of the ClinVar and gnomAD UTR variants for the UTR Visualization apps with the UTRannotator. 

## Environment 
Use a conda environment named pipeline 

```sh
conda env create --file ../envs/pipeline.yml pipeline
conda activate
```

Other dependencies, need to install `tabix`, `docker`, `snakemake`. 

## Install Docker 

```sh 
docker pull ensemblorg/ensembl-vep

# To run vep standalone 
docker run -t -i ensemblorg/ensembl-vep ./vep

# To view the image interactively running vep 
docker run -t -i -v $(pwd)/vep_data/:/opt/vep/.vep/ ensemblorg/ensembl-vep bash
```


## Sorfs DB 

For reproducibility, here are the options and attributes for the sORFS.download

```
<!DOCTYPE Query><Query client="true" processor="TSV" limit="-1" header="1"><Dataset name="BioMart" config="Human"><Filter name="human__annotation_104" value="5UTR" filter_list=""/><Attribute name="human__sorf_id_104"/><Attribute name="human__chr_104"/><Attribute name="human__sorf_end_104"/><Attribute name="human__strand_104"/><Attribute name="human__sorf_begin_104"/><Attribute name="human__upstream_gene_distance_104"/><Attribute name="human__downstream_gene_distance_104"/><Attribute name="human__tr_seq_104"/><Attribute name="human__id_104"/><Attribute name="human__start_codon_104"/></Dataset></Query>
```





## Key scripts for the pre-processing pipeline. 

Create a snakemake pipeline for the following scripts.

1. `download_clinvar.sh` : Downloads the weekly update of clinvar
2. `download_gnomAD.sh` : Downloads gnomAD (Need to only do this once)
3. `find_UTR_regions.py` : Filter to MANE UTR regions 
4. `filter_to_UTRs.sh` : Runs tabix to filter UTR regions to the necessary.
5. `run_docker_vep.sh` : Runs docker on the whole set + UTR annotator
6. `model.py` : Describing the sql alch
7. `create_sqlite_db.py` : Creates / updates a sqlite3 database from gnomAD and VEP for use in the web_server. 


Commandline Argument Options : 

`gnomad version`, `latest clinvar`, `assembly`, `output_db_sqlite3`, `mane_version`, `ENSEMBL_V` (Ensure that Ensembl version matches the mane version). 


Aim afterwards is to turn this into a cron job that can do this regularly. 



