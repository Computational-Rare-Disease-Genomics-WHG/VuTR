#!/bin/bash

# The following script deploys the VuTR container using podman.
# It requires two arguments: --features-db and --variant-store-db,
# which specify the paths to the respective database files on the host system.


# Check if podman is installed
if ! command -v podman &> /dev/null; then
    echo "Error: podman is not installed or not in PATH."
    exit 1
fi

# Default values
FEATURES_DB=""
VARIANT_STORE_DB=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --features-db)
            FEATURES_DB="$2"
            shift 2
            ;;
        --variant-store-db)
            VARIANT_STORE_DB="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 --features-db FEATURES_DB_PATH --variant-store-db VARIANT_STORE_DB_PATH"
            echo "Deploy the VuTR container."
            echo ""
            echo "Options:"
            echo "  --features-db FEATURES_DB_PATH        Path to the features.db file (required)"
            echo "  --variant-store-db VARIANT_STORE_DB_PATH  Path to the variant_store.db file (required)"
            echo "  --help, -h                            Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information."
            exit 1
            ;;
    esac
done

# Check if required arguments are provided
if [ -z "$FEATURES_DB" ] || [ -z "$VARIANT_STORE_DB" ]; then
    echo "Error: Both --features-db and --variant-store-db are required."
    echo "Use --help for usage information."
    exit 1
fi

# Check if the files exist
if [ ! -f "$FEATURES_DB" ]; then
    echo "Error: '$FEATURES_DB' does not exist or is not a file."
    exit 1
fi

if [ ! -f "$VARIANT_STORE_DB" ]; then
    echo "Error: '$VARIANT_STORE_DB' does not exist or is not a file."
    exit 1
fi

# Pull the latest image from GitHub Container Registry
if ! podman pull ghcr.io/computational-rare-disease-genomics-whg/vutr:latest; then
    echo "Error: Failed to pull the image."
    exit 1
fi

# Run the container
if ! podman run -d \
  --replace \
  --name vutr \
  --restart unless-stopped \
  --hostname vutr \
  -p 8080:8080 \
  -v "$FEATURES_DB:/db/features.db" \
  -v "$VARIANT_STORE_DB:/db/variant_store.db" \
  --health-cmd "curl -f http://localhost:8080/health" \
  --health-interval 1m30s \
  --health-timeout 30s \
  --health-retries 5 \
  --health-start-period 30s \
  -e PYTHONUNBUFFERED=1 \
  ghcr.io/computational-rare-disease-genomics-whg/vutr:latest; then
    echo "Error: Failed to run the container."
    exit 1
fi

echo "Container 'VuTR' started successfully."