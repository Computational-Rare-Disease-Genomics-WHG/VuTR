# Pipeline

This the processing scripts to annotate all of the ClinVar and gnomAD UTR variants for the UTR Visualization apps with the UTRannotator to a SQLite3 database.

## Install Docker

We can interactivately look at the Ensembl docker image using.

```sh
docker pull ensemblorg/ensembl-vep

# To run vep standalone
docker run -t -i ensemblorg/ensembl-vep ./vep

# To view the image interactively running vep
docker run -t -i -v $(realpath ../../data/pipeline/vep_data/)/vep_data/:/opt/vep/.vep/ ensemblorg/ensembl-vep bash
```

## Key scripts for the pre-processing pipeline

Create a snakemake pipeline for the following scripts.


1. `convert_MANE_to_csv.R` : Converts the gff file to a CSV file.
2. `filter_to_UTRs.R` : Filter to MANE UTR exons
3. `run_tabix.sh` : Runs tabix to filter variants to UTR based on the bed file of the previous script
4. `run_docker_vep.sh` : Runs docker on the whole set + UTR annotator
5. `parse_vep_output.py` : Cleans up the UTR annotator output file from VE!P to a nicer tsv file.
6. `create_sqlite_db.py` : Reads all of the processed data and converts to a SQLite3 database that will be used in the server.

Commandline Argument Options :

`gnomad version`, `latest clinvar`, `assembly`, `output_db_sqlite3`, `mane_version`, `ENSEMBL_V` (Ensure that Ensembl version matches the mane version).

Aim afterwards is to turn this into a cron job that can do this regularly.
