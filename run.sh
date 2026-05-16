#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/config/nextflow.config"
PARAMS="$SCRIPT_DIR/config/params.yaml"
MAIN="$SCRIPT_DIR/main.nf"

# Directorio de trabajo de Nextflow.
# Prioridad: variable de entorno NXF_WORK > valor por defecto ($PWD/work).
# Se puede sobreescribir antes de llamar al script:
#   NXF_WORK=/ruta/alternativa bash run.sh
WORK_DIR="${NXF_WORK:-$PWD/work}"

# Verificar dependencias
if ! command -v nextflow &> /dev/null; then
    echo "Error: nextflow no encontrado en PATH." >&2
    exit 1
fi

if ! command -v docker &> /dev/null; then
    echo "Error: docker no encontrado en PATH." >&2
    exit 1
fi

# Verificar archivos de configuracion
for f in "$CONFIG" "$PARAMS" "$MAIN"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: archivo no encontrado: $f" >&2
        exit 1
    fi
done

echo "Directorio de trabajo: $WORK_DIR"

nextflow run "$MAIN" \
    -c "$CONFIG" \
    -params-file "$PARAMS" \
    -w "$WORK_DIR" \
    "$@"
