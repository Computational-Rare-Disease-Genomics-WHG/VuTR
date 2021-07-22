
# Taking inspiration from the approach taken in https://github.com/IanVermes/vep_api 
FROM ensemblorg/ensembl-vep

# Put python related project in its own workdir
USER root
#WORKDIR $OPT_SRC/ensembl-vep/webvep
## Install Python3.8 and pip
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3.8
    
# Run Flask API as an entry point for VEP

# TODO EDIT THE FOLLOWING TO CREATE THE REST API

# Create Workdir
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m virtualenv --python=/usr/bin/python3 $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# copy the dependencies file to the working directory
COPY requirements.txt .

# install dependencies
RUN pip3 install -r requirements.txt

# copy the content of the local src directory to the working directory
COPY ./app /app
COPY wsgi.py .
COPY setup_gunicorn.sh . 

# To run as gunicorn 
ENTRYPOINT ["./setup_gunicorn.sh"]
