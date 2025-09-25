#!/bin/bash
set -euo pipefail

# Deploy VuTR using rootless Podman + systemd user unit (no quadlets, no compose).

# Configuration
IMAGE="ghcr.io/computational-rare-disease-genomics-whg/vutr:latest"
CTR_NAME="vutr"
PORT=8080
UNIT_FILE="container-${CTR_NAME}.service"
SYSTEMD_USER_DIR="${HOME}/.config/systemd/user"

# Parse arguments
FEATURES_DB=""
VARIANT_STORE_DB=""

usage(){
  echo "Usage: $0 --features-db PATH --variant-store-db PATH"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --features-db)
      FEATURES_DB="${2:-}"
      shift 2
      ;;
    --variant-store-db)
      VARIANT_STORE_DB="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown argument: $1"
      usage
      ;;
  esac
done

# Validate
command -v podman >/dev/null || { echo "Error: podman not found"; exit 1; }
if [[ -z "$FEATURES_DB" || ! -f "$FEATURES_DB" ]]; then
  echo "Error: invalid features-db: $FEATURES_DB"
  exit 1
fi
if [[ -z "$VARIANT_STORE_DB" || ! -f "$VARIANT_STORE_DB" ]]; then
  echo "Error: invalid variant-store-db: $VARIANT_STORE_DB"
  exit 1
fi

mkdir -p "$SYSTEMD_USER_DIR"

echo ">> Pull image $IMAGE"
podman pull "$IMAGE"

# Stop old unit if exists
echo ">> Stopping old systemd unit, if any"
systemctl --user stop "$UNIT_FILE" 2>/dev/null || true

# Remove old container
echo ">> Removing existing container (if any)"
podman rm -f "$CTR_NAME" 2>/dev/null || true

# Run new container
echo ">> Running container $CTR_NAME"
podman run -d \
  --name "$CTR_NAME" \
  --hostname "$CTR_NAME" \
  -p "$PORT:$PORT" \
  -v "${FEATURES_DB}:/db/features.db:Z" \
  -v "${VARIANT_STORE_DB}:/db/variant_store.db:Z" \
  -e PYTHONUNBUFFERED=1 \
  "$IMAGE"

# Generate systemd unit
echo ">> Generating systemd unit"
podman generate systemd \
  --name "$CTR_NAME" \
  --new \
  --files \
  --restart-policy always

# Move the generated unit file to user units directory
echo ">> Installing systemd user unit"
mv -f "$UNIT_FILE" "$SYSTEMD_USER_DIR/$UNIT_FILE"

# Reload systemd user daemon
echo ">> Reloading user systemd"
systemctl --user daemon-reload

# Enable and start
echo ">> Enabling + starting $UNIT_FILE"
systemctl --user enable --now "$UNIT_FILE"

cat <<EOF

✅ Deployed VuTR under systemd --user

Service unit: ${UNIT_FILE}
Container:    ${CTR_NAME}
Image:        ${IMAGE}

Commands:
  systemctl --user status ${UNIT_FILE}
  journalctl --user -u ${UNIT_FILE} -f
  podman logs -f ${CTR_NAME}
  podman pull ${IMAGE} && systemctl --user restart ${UNIT_FILE}

Make sure your user has lingering enabled so the service restarts across logouts/reboots:
  loginctl enable-linger $USER

EOF