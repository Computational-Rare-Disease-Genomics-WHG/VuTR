# Running locally for testing and development

## Installation of dependencies 

- Docker
- Docker-VEP 
- UTR annotator through docker vep
```{bash}

```
- Python Env and requirements. 

```{bash}
cd server # Change to this exact directory

export FLASK_ENV=development 
export FLASK_APP=app 
flask run 
```
# Running in a production environment (VPS / Web server)

## Installing docker on the machine

```{bash}
docker-compose up -d 
```