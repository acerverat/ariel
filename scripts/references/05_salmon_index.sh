#!/bin/bash
set -e

# Genera el indice de Salmon para GRCh38 usando el transcriptoma de Gencode v42.
# El transcriptoma se descarga, se indexa y luego se elimina para ahorrar espacio.

rutaReferencias="$1"
outdir="$rutaReferencias/GRCh38_no_alt"
sentinel="$outdir/salmon_index/info.json"

if [[ -f "$sentinel" ]]; then
    echo "[salmon_index] Indice de Salmon ya presente, omitiendo generacion."
    exit 0
fi

echo "[salmon_index] Descargando transcriptoma Gencode v42..."
mkdir -p "$outdir"
cd "$outdir"

wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.transcripts.fa.gz
gunzip -f gencode.v42.transcripts.fa.gz

echo "[salmon_index] Generando indice de Salmon..."
docker run -u $(id -u):$(id -g) --rm \
    -v "$outdir":"$outdir" \
    -w "$outdir" \
    ariel-env \
    salmon index \
        -t gencode.v42.transcripts.fa \
        -i salmon_index \
        --threads 4

rm gencode.v42.transcripts.fa

echo "[salmon_index] Listo."
