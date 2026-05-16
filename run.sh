#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/config/nextflow.config"
PARAMS="$SCRIPT_DIR/config/params.yaml"
MAIN="$SCRIPT_DIR/main.nf"

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

nextflow run "$MAIN" \
    -c "$CONFIG" \
    -params-file "$PARAMS" \
    "$@"
