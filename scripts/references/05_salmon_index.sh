#!/bin/bash
set -e

# Genera el indice de Salmon para GRCh38 usando el transcriptoma de Gencode v42.
# El transcriptoma se descarga, se indexa y luego se elimina para ahorrar espacio.
# Tambien genera ensg_enst.tsv (gene_id, transcript_id, gene_name) desde la anotacion GTF.

rutaReferencias="$1"
outdir="$rutaReferencias/GRCh38_no_alt"

# ── Indice de Salmon ──────────────────────────────────────────────────────────

if [[ -f "$outdir/salmon_index/info.json" ]]; then
    echo "[salmon_index] Indice de Salmon ya presente, omitiendo generacion."
else
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
fi

# ── Tabla ENSG / ENST / gene_name ─────────────────────────────────────────────

if [[ -f "$outdir/geneId_transcriptId_geneName.tsv" ]]; then
    echo "[ensg_enst] Tabla ya presente, omitiendo generacion."
else
    gtf="$outdir/gencode.v42.annotation.gtf"
    if [[ ! -f "$gtf" ]]; then
        echo "[ensg_enst] ERROR: $gtf no encontrado. Ejecuta primero 02_gencode.sh." >&2
        exit 1
    fi

    echo "[ensg_enst] Generando geneId_transcriptId_geneName.tsv desde la anotacion Gencode v42..."
    gawk 'BEGIN{OFS="\t"; print "gene_id\ttranscript_id\tgene_name"}
          $3=="transcript" {
            match($0, /gene_id "([^"]+)"/, g)
            match($0, /transcript_id "([^"]+)"/, t)
            match($0, /gene_name "([^"]+)"/, n)
            print g[1], t[1], n[1]
          }' "$gtf" > "$outdir/geneId_transcriptId_geneName.tsv"

    echo "[ensg_enst] Listo. $(( $(wc -l < "$outdir/geneId_transcriptId_geneName.tsv") - 1 )) transcriptos escritos."
fi
