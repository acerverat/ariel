#!/bin/bash
set -e

# Indice de STAR 2.7.10b para GRCh38
# Requiere que las referencias de Cicero (01_cicero.sh) y Gencode (02_gencode.sh) ya esten instaladas.

rutaReferencias="$1"
outdir="$rutaReferencias/GRCh38_no_alt"
sentinel="$outdir/STAR_2.7.10b_index/SAindex"

if [[ -f "$sentinel" ]]; then
    echo "[star_index] Indice de STAR ya presente, omitiendo generacion."
    exit 0
fi

fasta="$outdir/cicero_references/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa"
gtf="$outdir/gencode.v42.annotation.gtf"

if [[ ! -f "$fasta" ]]; then
    echo "[star_index] ERROR: FASTA no encontrado. Ejecuta primero 01_cicero.sh."
    exit 1
fi

if [[ ! -f "$gtf" ]]; then
    echo "[star_index] ERROR: GTF no encontrado. Ejecuta primero 02_gencode.sh."
    exit 1
fi

echo "[star_index] Generando indice de STAR..."
mkdir -p "$outdir/STAR_2.7.10b_index"

docker run -u $(id -u):$(id -g) --rm \
    -v "$outdir":"$outdir" \
    -w "$outdir" \
    acerverat/ariel-env:latest \
    STAR \
    --runMode genomeGenerate \
    --genomeDir STAR_2.7.10b_index \
    --genomeFastaFiles "$fasta" \
    --sjdbGTFfile "$gtf" \
    --runThreadN 4

echo "[star_index] Listo."
