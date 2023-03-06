# Not working
# Create a dockerfile with the following dependencies 
# to run the pipeline reproducibly

FROM alpine:latest

# Install system dependencies for VEP, R and Bioconductor
RUN apk update && \
    apk upgrade && \
    apk add --no-cache perl wget unzip R R-dev make gcc g++ libxml2 libxml2-dev curl perl-dbd-mysql
# Install VEP
RUN wget https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/104.0.tar.gz && \
    tar -zxvf 104.0.tar.gz && \
    mv ensembl-vep-release-104.0 vep && \
    rm 104.0.tar.gz && \
    cd vep && \
    perl INSTALL.pl --AUTO cpanminus && \
    perl INSTALL.pl --AUTO Tools && \
    cd ..

# Install Conda
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -O miniconda.sh && \
    /bin/sh miniconda.sh -bfp /usr/local && \
    rm miniconda.sh && \
    /usr/local/bin/conda clean -tipsy && \
    ln -s /usr/local/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /usr/local/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Copy the environment.yml file into the container
COPY pipeline-env.yaml .

# Create a new Conda environment based on the environment file
RUN /opt/conda/bin/conda env create -f pipeline-env.yaml && \
    /opt/conda/bin/conda clean -afy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.profile && \
    echo "conda activate vep-r" >> ~/.profile

# Install R packages from CRAN and Bioconductor
RUN R -e "install.packages(c('magrittr', 'data.table'), repos='https://cran.rstudio.com/')" && \
    R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')" && \
    R -e "BiocManager::install(c('rtracklayer'))"

ENV PATH /opt/conda/bin:$PATH
