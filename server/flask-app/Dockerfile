# Dockerfile
# Elston N. D'Souza

# set base image (host OS)
# set base image (host OS)
FROM python:3.8-slim-buster

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
CMD ["gunicorn"  , "-b", "0.0.0.0:8080", "wsgi:app"]
