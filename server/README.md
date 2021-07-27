# Running locally for testing and development

## Installation of dependencies 

- Docker
- Docker-VEP 
- UTR annotator through docker vep
- Create a virtual env and `cd flask-app; pip3 install -r requirements.txt`

- Running the flask application server 

```bash
cd flask-app/
export FLASK_ENV=development 
export FLASK_APP=app 
flask run 
```
# Running in a production environment (VPS / Web server)

## Installing docker on the machine

```bash
docker-compose up -d 
```
