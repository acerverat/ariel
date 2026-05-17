#!/bin/bash
set -euo pipefail

# Genera config/params.yaml a partir de rutas introducidas interactivamente.
# Ejecutar desde el directorio raiz del proyecto:
#   bash scripts/generaParams.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$SCRIPT_DIR/../config/params.yaml"

ask() {
    local prompt="$1"
    local default="${2:-}"
    local value

    if [[ -n "$default" ]]; then
        read -rp "$prompt [$default]: " value
        echo "${value:-$default}"
    else
        read -rp "$prompt: " value
        if [[ -z "$value" ]]; then
            echo "Error: campo obligatorio." >&2
            exit 1
        fi
        echo "$value"
    fi
}

echo "============================================"
echo " Configuracion de parametros de ARIEL"
echo " Salida: $OUT"
echo "============================================"
echo ""

runSampleSheet=$(ask  "Ruta del SampleSheet TSV (Sample/R1/R2)")
referenceDir=$(ask    "Directorio de referencias")
tpmPanel=$(ask        "Tabla de TPMs de referencia (panel de genes)")
resultsDir=$(ask      "Directorio de resultados")
reportsDir=$(ask      "Directorio de reportes de QC")
threadsSTAR=$(ask     "Hilos para STAR" "4")
threadsRascall=$(ask  "Hilos para RaScALL" "4")

cat > "$OUT" <<EOF
# Parametros del flujo de trabajo ARIEL.
# Generado por scripts/generaParams.sh
# Pasar a Nextflow con: nextflow run main.nf -c config/nextflow.config -params-file config/params.yaml

# Rutas de entrada
runSampleSheet: "$runSampleSheet"
referenceDir:   "$referenceDir"
tpmPanel:       "$tpmPanel"

# Rutas de salida
resultsDir: "$resultsDir"
reportsDir: "$reportsDir"

# Hilos
threadsSTAR:    $threadsSTAR
threadsRascall: $threadsRascall
EOF

echo ""
echo "============================================"
echo " config/params.yaml generado."
echo "============================================"
