
# Dockerfile
# Elston N. D'Souza

# Multi-stage build for distroless
FROM python:3.9.7-slim-buster AS builder

# Set working directory
WORKDIR /app

# Install dependencies with all required packages
COPY requirements.txt .
RUN pip install --no-cache-dir --target=/app/packages "setuptools<81" wheel pip && \
    pip install --no-cache-dir --target=/app/packages -r requirements.txt

# Runtime stage using distroless
FROM gcr.io/distroless/python3-debian11

# Copy Python packages from builder
COPY --from=builder /app/packages /app/packages

# Set environment variables
ENV PYTHONPATH="/app/packages"
ENV FLASK_DEBUG=false
ENV FLASK_ENV=production
ENV PYTHONUNBUFFERED=1

# Set working directory
WORKDIR /app

# Copy application files
COPY ./app ./app
COPY wsgi.py .

# To run as gunicorn
ENTRYPOINT ["python", "-m", "gunicorn", "wsgi:app", "-w", "3", "--threads", "1", "-b", "0.0.0.0:8080"]
