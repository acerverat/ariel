#!/usr/bin/env bash
# run_local.sh — for local development and testing.
# Copy to your run folder, fill in the values marked <...>, then execute.
# Do not commit your filled-in copy to the repo.

set -euo pipefail

# ── Fill these in before running ───────────────────────────────────────────────
ARIEL_DIR="$HOME/ARIEL"          # e.g. $HOME/ARIEL
REFERENCE_DIR="$HOME/Prueba_Ariel/referencias/GRCh38_no_alt"   # absolute path to GRCh38_no_alt/
TPM_PANEL="$HOME/Prueba_Ariel/datos/tpmPanel_fixed.tsv"         # absolute path to tpmPanel.tsv
SAMPLE_SHEET="$HOME/Prueba_Ariel/datos/sample_sheet.tsv"          # absolute path to sample_sheet.tsv
BASE_DIR="$HOME/Prueba_Ariel/RESULTADOS"             # parent directory where run folders will be created
WORK_DIR=""                                # leave empty to default to <RUN_DIR>/work
THREADS_STAR=1
THREADS_RASCALL=1
THREADS_FUNGI=15

# ── Validation ─────────────────────────────────────────────────────────────────
for var in ARIEL_DIR REFERENCE_DIR TPM_PANEL SAMPLE_SHEET BASE_DIR; do
    [[ "${!var}" == "<"* ]] && { echo "Error: fill in $var before running." >&2; exit 1; }
done

if [[ "${1:-}" == "--name" ]]; then
    RUN_NAME="${2:?'--name requires a value'}"
else
    read -rp "Nombre de la corrida: " RUN_NAME
    [[ -z "$RUN_NAME" ]] && { echo "Error: el nombre de la corrida no puede estar vacío." >&2; exit 1; }
fi

# ── Derived paths ──────────────────────────────────────────────────────────────
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${BASE_DIR}/${RUN_NAME}"
[[ -z "$WORK_DIR" ]] && WORK_DIR="${BASE_DIR}/work"

# ── Create run directory structure ────────────────────────────────────────────
mkdir -p "${RUN_DIR}/resultados" "${RUN_DIR}/reportesQC" "${RUN_DIR}/logs"

# ── Write params.yaml ──────────────────────────────────────────────────────────
cat > "${RUN_DIR}/params.yaml" <<EOF
# Rutas de entrada
runSampleSheet: "$SAMPLE_SHEET"
referenceDir:   "$REFERENCE_DIR"
tpmPanel:       "$TPM_PANEL"

# Rutas de salida
resultsDir: "${RUN_DIR}/resultados"

# Hilos
threadsSTAR:    $THREADS_STAR
threadsRascall: $THREADS_RASCALL
threadsFungi:   $THREADS_FUNGI
EOF

# ── Run pipeline ───────────────────────────────────────────────────────────────
nextflow run "${ARIEL_DIR}/main.nf" \
    -name        "${RUN_NAME}_${TIMESTAMP}" \
    -w           "$WORK_DIR" \
    -params-file "${RUN_DIR}/params.yaml" \
    -with-report "${RUN_DIR}/logs/report.html" \
    -with-trace  "${RUN_DIR}/logs/trace.txt" \
    -resume

echo "Done. Results in ${RUN_DIR}/resultados/, logs in ${RUN_DIR}/logs/"
