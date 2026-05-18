#!/bin/bash
set -e

# Descarga el indice de Salmon para GRCh38 desde refgenie (salmon_partial_sa_index).
# El archivo incluye el GTF de anotacion que se usa para generar geneId_transcriptId_geneName.tsv.

rutaReferencias="$1"
outdir="$rutaReferencias/GRCh38_no_alt"

REFGENIE_URL="http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/salmon_partial_sa_index?tag=default"

# ── Indice de Salmon ──────────────────────────────────────────────────────────

if [[ -f "$outdir/salmon_index/info.json" ]]; then
    echo "[salmon_index] Indice de Salmon ya presente, omitiendo descarga."
else
    mkdir -p "$outdir"
    cd "$outdir"

    echo "[salmon_index] Descargando indice de Salmon desde refgenie..."
    wget -c "$REFGENIE_URL" -O salmon_partial_sa_index.tgz

    echo "[salmon_index] Extrayendo..."
    tar -xzf salmon_partial_sa_index.tgz

    # El archivo refgenie extrae en salmon_partial_sa_index/default/salmon_partial_sa_index/
    idx_src=$(find . -name "info.json" | head -1 | xargs dirname)
    mv "$idx_src" "$outdir/salmon_index"

    rm -f salmon_partial_sa_index.tgz
    rm -rf salmon_partial_sa_index/

    echo "[salmon_index] Listo."
fi

# ── Tabla ENSG / ENST / gene_name ─────────────────────────────────────────────

if [[ -f "$outdir/geneId_transcriptId_geneName.tsv" ]]; then
    echo "[ensg_enst] Tabla ya presente, omitiendo generacion."
else
    gtf="$(ls "$outdir/salmon_index/"*.gtf 2>/dev/null | head -1)"
    if [[ -z "$gtf" ]]; then
        echo "[ensg_enst] ERROR: no se encontro GTF en $outdir/salmon_index/. Ejecuta primero este mismo script para generar el indice." >&2
        exit 1
    fi

    echo "[ensg_enst] Generando geneId_transcriptId_geneName.tsv desde $gtf ..."
    gawk 'BEGIN{OFS="\t"; print "gene_id\ttranscript_id\tgene_name"}
          $3=="transcript" {
            match($0, /gene_id "([^"]+)"/, g)
            match($0, /transcript_id "([^"]+)"/, t)
            match($0, /gene_name "([^"]+)"/, n)
            print g[1], t[1], n[1]
          }' "$gtf" > "$outdir/geneId_transcriptId_geneName.tsv"

    echo "[ensg_enst] Listo. $(( $(wc -l < "$outdir/geneId_transcriptId_geneName.tsv") - 1 )) transcriptos escritos."
fi
