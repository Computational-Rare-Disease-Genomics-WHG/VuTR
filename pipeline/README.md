# Pipeline 


This is a description of the upstream pipeline to annotate all of the ClinVar and gnomAD UTR variants for the UTR Visualization apps with the UTRannotator. 


## Environment 

Use a conda environment named pipeline 

```sh
conda env create --file ../envs/pipeline_environment.yml pipeline

conda activate
```

## Key scripts 

Create a snakemake pipeline for the following scripts.

1. `download_clinvar.sh` : Downloads the weekly update of clinvar
2. `download_gnomAD.sh` : Downloads gnomAD (Need to only do this once)
3. `run_docker_vep.sh` : Runs docker on the whole set
4. `create_sqlite_db.py` : Creates a sqlite3 database for use in the web_server 

Commandline Argument Options : 
`gnomad version`, `latest clinvar`, `assembly`, `output_db_sqlite3`

## Install Docker 

```sh 
docker pull ensemblorg/ensembl-vep

# To run vep standalone 
docker run -t -i ensemblorg/ensembl-vep ./vep

# To view the image interactively running vep 
docker run -t -i -v $(pwd)/UTRannotator/:/opt/vep/.vep/Plugins ensemblorg/ensembl-vep bash
```
