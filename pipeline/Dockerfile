# Create a dockerfile with the following dependencies 
# to run the pipeline reproducibly

# To build and install and send to dockerhub

# docker login (Need to do this once)

# Build for both arm64 and amd64
# docker buildx build --platform linux/amd64,linux/arm64 -t elstonndsouza/docker-vutr-pipeline:latest --push .

# Development notes: 

# Was building this for on the ARM M1 however
# we run into issues with numpy for arm64. 
# This is okay because ARM M1 macs can 
# run docker images.

# Installing hail was also a goal however, we run into 
# issues due to the fact that the python version is a 
# bit too high at python3.9, hail currently only works
# stably for python3.7. Perhaps


###################################################
# Stage 1 - docker container to build ensembl-vep #
###################################################
FROM ensemblorg/ensembl-vep:release_103

# Install  R
# Install bedtools, htslib and samtools, bedtools 
# and other dependencies

USER root
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
    apt-get install -y r-base \
     r-base-dev \
     samtools \
     tabix \
     bcftools \ 
     bedtools \
     openjdk-8-jre-headless \
     g++ \
     libopenblas-base \
     liblapack3 



# Install R packages from CRAN and Bioconductor
RUN R -e "install.packages(c('magrittr', 'optparse', 'logger' ,'data.table', 'BiocManager'), repos='https://cran.rstudio.com/')" && \
    R -e "BiocManager::install(c('rtracklayer', 'seqinr'))"

WORKDIR /root

# Install Python 
RUN apt-get install -y python3-pip python3-dev python3-venv

# Install python packages
RUN pip3 install sqlalchemy numpy pandas scipy matplotlib
