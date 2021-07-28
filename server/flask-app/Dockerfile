# Dockerfile 
# Elston N. D'Souza

# set base image (host OS)
FROM ubuntu:18.04
RUN apt-get update && apt-get install \
  -y --no-install-recommends python3 python3-virtualenv

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

