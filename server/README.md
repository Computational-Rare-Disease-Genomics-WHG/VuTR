# Running locally for testing and development

## Installation of dependencies 

- Docker
- Docker-VEP 
- UTR annotator through docker vep
- Download Ensembl Cache
- Create and source an virtual env and then `cd flask-app; pip3 install -r requirements.txt`


## Running the flask application server 

```bash
# Once cd'd into flask-app, then 
export FLASK_ENV=development 
export FLASK_APP=app 
flask run 
```
## Running through docker 

```bash
sudo docker build -t utrapp . 
sudo docker run --rm -d -p 8080:8080 utrapp

# To kill the container 

# find the name of the container 
docker ps 
docker stop "<CONTAINER_NAME>"
```

# Running in a production environment (VPS / Web server)

## Installing docker on the machine

```bash
docker-compose up -d 
```
