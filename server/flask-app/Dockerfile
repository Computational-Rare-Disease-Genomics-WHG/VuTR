# Dockerfile
# Elston N. D'Souza

# set base image (host OS)
# set base image (host OS)
FROM python:3.9.7-slim-buster

# Create Workdir
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# copy the dependencies file to the working directory
COPY requirements.txt .
RUN pip install -r requirements.txt

# copy the content of the local src directory to the working directory
COPY ./app /app
COPY wsgi.py .
COPY setup_gunicorn.sh .
ENV FLASK_ENV='production'

# To run as gunicorn
ENTRYPOINT ["./setup_gunicorn.sh"]
