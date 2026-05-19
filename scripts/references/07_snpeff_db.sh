#!/bin/bash
set -e

# Descarga la base de datos de SnpEff para el genoma indicado.
# El genoma por defecto es GRCh38.99; puede sobreescribirse con el segundo argumento.
# La base de datos se guarda en GRCh38_no_alt/snpeff_db/ y se monta en
# /snpeff_data dentro del contenedor en tiempo de ejecucion del pipeline.

rutaReferencias="$1"
genome="${2:-GRCh38.p14}"
outdir="$rutaReferencias/GRCh38_no_alt"
snpeffdb="$outdir/snpeff_db"
sentinel="$snpeffdb/${genome}/snpEffectPredictor.bin"

if [[ -f "$sentinel" ]]; then
    echo "[snpeff_db] Base de datos $genome ya presente, omitiendo descarga."
    exit 0
fi

echo "[snpeff_db] Descargando base de datos $genome..."
mkdir -p "$snpeffdb"

docker run --rm \
    -v "$snpeffdb":/snpeff_data \
    ariel-env \
    java -jar /opt/snpEff/snpEff.jar download \
        -dataDir /snpeff_data \
        -v "$genome"

echo "[snpeff_db] Listo."
