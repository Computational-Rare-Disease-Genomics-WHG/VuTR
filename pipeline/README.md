## UTR Visualization App Pipeline 

This directory  contains the code to query and process the data as a pre-requisite to the pipline
Outputs 

First, is the `downloads` folders that queries all the tools used in the UTR-Visualization Application

Second, is the process 

To run the UTR-Visualization application, we need to download and process all of the data required by the webserver application. 

This document is brief outline of how to run (before moving to a Snakemake pipeline.)

### Pre-requisites 

Ensure that coreutils are install if you are using macOS

```bash
brew install coreutils
```

Other dependencies include `wget`, `docker (Optional if not using vep)`, `vep`. 

Ensure that you also have ENSEMBL's *homo_sapiens* cache installed for the version of Ensembl being used.

All other python dependencies are packaged in conda. 

Go to the root directory and install the requirements from the conda environment

```bash
conda env create -f pipeline-env.yml
conda activate utr-app-pipeline
```



## Running the pipline 

The aim is to make this automated through Snakemake / Nextflow.

Firstly, we need to create the output directory structure where all of the output files will be held. To do this, run the following command. 

```bash 
bash create_data_dirs.sh
```

Now, if we look at the root of this github directory, a new folder called `data` should be created. 

Start with downloads, pass a mane version (tested on 0.93).

```bash 
cd ~/downloads
bash download_all.sh --mane_v MANE_V # Optional argument, default mane version is 0.93
```

Run scripts within each of the directories in the following order. 

- downloads
- process_mane
- variants
- database

Commands to run each of the scripts are in the directories.


