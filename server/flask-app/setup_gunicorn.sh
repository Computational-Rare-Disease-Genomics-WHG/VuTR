#!/bin/sh
source /opt/venv/bin/activate
gunicorn wsgi:app -w 3 --threads 2 -b 0.0.0.0:8080
