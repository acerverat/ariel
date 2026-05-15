#!/bin/bash
set -e

# Indice de Salmon (GRCh38, partial SA index) desde refgenie

rutaReferencias="$1"
outdir="$rutaReferencias/GRCh38_no_alt"
sentinel="$outdir/salmon_index/info.json"

if [[ -f "$sentinel" ]]; then
    echo "[salmon_index] Indice de Salmon ya presente, omitiendo descarga."
    exit 0
fi

echo "[salmon_index] Descargando indice de Salmon desde refgenie..."
mkdir -p "$outdir"
cd "$outdir"

wget -c "http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/salmon_partial_sa_index?tag=default" \
    -O salmon_index.tar.gz

tar -xzf salmon_index.tar.gz

# refgenie extrae a salmon_partial_sa_index/default/
mv salmon_partial_sa_index/default salmon_index
rm -rf salmon_partial_sa_index salmon_index.tar.gz

echo "[salmon_index] Listo."
