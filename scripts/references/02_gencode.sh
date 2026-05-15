#!/bin/bash
set -e

# Anotaciones Gencode v42 (GTF y GFF3) para GRCh38

rutaReferencias="$1"
outdir="$rutaReferencias/GRCh38_no_alt"
sentinel="$outdir/gencode.v42.annotation.gtf"

if [[ -f "$sentinel" ]]; then
    echo "[gencode] Anotaciones ya presentes, omitiendo descarga."
    exit 0
fi

echo "[gencode] Descargando anotaciones Gencode v42..."
mkdir -p "$outdir"
cd "$outdir"

wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gff3.gz

gunzip -f gencode.v42.annotation.gtf.gz
gunzip -f gencode.v42.annotation.gff3.gz

echo "[gencode] Listo."
