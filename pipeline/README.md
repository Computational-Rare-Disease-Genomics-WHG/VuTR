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

## Key scripts for the pre-processing pipeline. 

Create a snakemake pipeline for the following scripts.

1. `download_clinvar.sh` : Downloads the weekly update of clinvar
2. `download_gnomAD.sh` : Downloads gnomAD (Need to only do this once)
3. `find_UTR_regions.py` : Filter to MANE UTR regions 
4. `filter_to_UTRs.sh` : Runs tabix to filter UTR regions to the necessary.
5. `run_docker_vep.sh` : Runs docker on the whole set + UTR annotator
6. `create_sqlite_db.py` : Creates / updates a sqlite3 database from gnomAD and VEP for use in the web_server. 

Commandline Argument Options : 

`gnomad version`, `latest clinvar`, `assembly`, `output_db_sqlite3`, `mane_version`, `ENSEMBL_V` (Ensure that Ensembl version matches the mane version). 


Aim afterwards is to turn this into a cron job that can do this regularly. 



