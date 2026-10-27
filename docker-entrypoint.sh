#!/bin/sh
# Shiny Server launches each app as the `shiny` user and does not carry the
# container environment through to the R process, so `docker run -e ECMX_MODE=web`
# alone silently leaves the app in unrestricted local mode. R does read
# .Renviron from the application directory at startup, so translate the
# container environment into that file before starting the server.
set -eu

APP_DIR=/srv/shiny-server/ecmexplainer

{
  echo "ECMX_MODE=${ECMX_MODE:-web}"
  echo "ECMX_MAX_UPLOAD_MB=${ECMX_MAX_UPLOAD_MB:-50}"
  echo "ECMX_MAX_SAMPLES=${ECMX_MAX_SAMPLES:-500}"
  echo "ECMX_MAX_TRAIN=${ECMX_MAX_TRAIN:-200}"
  echo "ECMX_MAX_PER_PHENO=${ECMX_MAX_PER_PHENO:-60}"
  echo "ECMX_WORKERS=${ECMX_WORKERS:-2}"
  echo "ECMX_CORES=${ECMX_CORES:-2}"
} > "$APP_DIR/.Renviron"

echo "ECMExplainer starting in ${ECMX_MODE:-web} mode:"
sed 's/^/  /' "$APP_DIR/.Renviron"

exec /usr/bin/shiny-server
