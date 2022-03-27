# Dockerfile
# Elston N. D'Souza

# set base image (host OS)
# Taking inspiration from the approach taken in https://github.com/IanVermes/vep_api
FROM python:3.8-slim-buster


RUN apt-get update && apt-get install \
  -y --no-install-recommends python3 python3-virtualenv

# Create Workdir
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# copy the dependencies file to the working directory
COPY requirements.txt .


# copy the content of the local src directory to the working directory
COPY ./app /app
COPY wsgi.py .
COPY setup_gunicorn.sh .
ENV FLASK_ENV='production'

# To run as gunicorn
ENTRYPOINT ["./setup_gunicorn.sh"]
