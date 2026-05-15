#!/bin/bash
set -e

# Referencias de Cicero y RNApeg
# Descarga desde: https://github.com/stjude/CICERO#reference

rutaReferencias="$1"
outdir="$rutaReferencias/GRCh38_no_alt"
sentinel="$outdir/cicero_references/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa"

if [[ -f "$sentinel" ]]; then
    echo "[cicero] Referencias ya presentes, omitiendo descarga."
    exit 0
fi

echo "[cicero] Descargando referencias de Cicero/RNApeg..."
mkdir -p "$outdir"
cd "$outdir"

wget -c https://zenodo.org/records/5088371/files/reference.tar.gz
tar -xvzf reference.tar.gz
mv reference cicero_references/
rm reference.tar.gz

echo "[cicero] Listo."
