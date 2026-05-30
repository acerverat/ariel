#!/bin/bash
set -e

# Descarga el VCF de ClinVar (GRCh38) para anotacion con SnpSift.
# Se guarda junto a la base de datos de SnpEff en snpeff_db/
# y se monta en /snpeff_data dentro del contenedor.

rutaReferencias="$1"
outdir="$rutaReferencias/GRCh38_no_alt/snpeff_db"
mkdir -p "$outdir"

if [[ -f "$outdir/clinvar.vcf.gz.tbi" ]]; then
    echo "[clinvar] ClinVar ya presente, omitiendo descarga."
else
    echo "[clinvar] Descargando ClinVar VCF (GRCh38)..."
    wget -c "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" \
         -O "$outdir/clinvar.vcf.gz"
    wget -c "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi" \
         -O "$outdir/clinvar.vcf.gz.tbi"
    echo "[clinvar] Listo."
fi

if [[ -f "$outdir/MANE_select.tsv" ]]; then
    echo "[mane] MANE Select ya presente, omitiendo descarga."
else
    echo "[mane] Descargando MANE Select (GRCh38)..."
    # Actualizar la URL si hay una nueva version en:
    # https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
    wget -c "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz" \
         -O "$outdir/MANE_select.tsv.gz"
    gunzip "$outdir/MANE_select.tsv.gz"
    echo "[mane] Listo."
fi
