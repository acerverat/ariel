#!/bin/bash
set -e

# Base de datos de FusionCatcher
# Se descarga usando el contenedor ariel-env y se monta en fusioncatcher_db/.
# En tiempo de ejecucion del pipeline, este directorio se monta en
# /opt/fusioncatcher/data/current dentro del contenedor.

rutaReferencias="$1"
outdir="$rutaReferencias/GRCh38_no_alt"
sentinel="$outdir/fusioncatcher_db/current/version.txt"

if [[ -f "$sentinel" ]]; then
    echo "[fusioncatcher_db] Base de datos ya presente, omitiendo descarga."
    exit 0
fi

echo "[fusioncatcher_db] Descargando base de datos de FusionCatcher..."
mkdir -p "$outdir/fusioncatcher_db"

docker run --rm \
    -v "$outdir/fusioncatcher_db":/opt/fusioncatcher/data \
    ariel-env \
    bash -c "/opt/fusioncatcher/data/download-human-db.sh"

echo "[fusioncatcher_db] Listo."
