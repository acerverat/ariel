#!/bin/bash
set -e

# Base de datos de FusionCatcher
# El script de descarga se extrae del contenedor en tiempo de ejecucion para
# que siempre este actualizado. Se ejecuta en el HOST (no dentro de Docker)
# para que las descargas parciales persistan entre reinicios. Un wrapper de
# wget inyectado via PATH agrega reanudacion (-c) y reintentos sin tocar el
# script original de FusionCatcher.
#
# En tiempo de ejecucion del pipeline, fusioncatcher_db/ se monta en
# /opt/fusioncatcher/data dentro del contenedor.

rutaReferencias="$1"
outdir="$rutaReferencias/GRCh38_no_alt"
dbdir="$outdir/fusioncatcher_db"
sentinel="$dbdir/current/version.txt"

if [[ -f "$sentinel" ]]; then
    echo "[fusioncatcher_db] Base de datos ya presente, omitiendo descarga."
    exit 0
fi

mkdir -p "$dbdir"

# Extraer el script de descarga directamente del contenedor para que
# cualquier actualizacion de FusionCatcher se refleje automaticamente.
# Se ejecuta en el HOST para que los archivos parciales persistan en disco.
DOWNLOAD_SCRIPT="$dbdir/download-human-db.sh"
echo "[fusioncatcher_db] Extrayendo script de descarga del contenedor..."
docker run --rm acerverat/ariel-env:latest \
    cat /opt/fusioncatcher/data/download-human-db.sh > "$DOWNLOAD_SCRIPT"
chmod +x "$DOWNLOAD_SCRIPT"

echo "[fusioncatcher_db] Descargando base de datos de FusionCatcher..."
cd "$dbdir"
bash "$DOWNLOAD_SCRIPT"

echo "[fusioncatcher_db] Listo."
