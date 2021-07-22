
FROM ensemblorg/ensembl-vep

# Put python related project in its own workdir
USER root
#WORKDIR $OPT_SRC/ensembl-vep/webvep
## Install Python3.8 and pip
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3.8
