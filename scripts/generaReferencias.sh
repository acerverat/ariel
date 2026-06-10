#!/bin/bash
set -e

# Script maestro para la instalacion de referencias de ARIEL.
# Llama a cada subscript de forma independiente para que puedan
# reejecutarse individualmente si alguno falla.
#
# Uso: bash generaReferencias.sh /ruta/completa/referencias
#
# Subscripts disponibles (en scripts/references/):
#   01_cicero.sh          Referencias de Cicero y RNApeg
#   02_gencode.sh         Anotaciones Gencode v42
#   03_star_index.sh      Indice de STAR (requiere 01 y 02)
#   04_fusioncatcher_db.sh Base de datos de FusionCatcher
#   05_salmon_index.sh    Indice de Salmon
#   06_snpeff_db.sh       Base de datos de SnpEff
#   07_clinvar.sh         VCF de ClinVar y tabla MANE Select

if [[ -z "$1" ]]; then
    echo "No se introdujo una ruta."
    echo "Uso: bash generaReferencias.sh /la_ruta_completa"
    echo "Ejemplo: bash generaReferencias.sh /home/yun/Documentos/referencias"
    exit 1
fi

rutaReferencias="$1"
SCRIPTS_DIR="$(dirname "$(realpath "$0")")/references"

echo "=========================================="
echo " Instalacion de referencias de ARIEL"
echo " Directorio: $rutaReferencias"
echo "=========================================="
echo ""

bash "$SCRIPTS_DIR/01_cicero.sh"          "$rutaReferencias"
bash "$SCRIPTS_DIR/02_gencode.sh"         "$rutaReferencias"
bash "$SCRIPTS_DIR/03_star_index.sh"      "$rutaReferencias"
bash "$SCRIPTS_DIR/04_fusioncatcher_db.sh" "$rutaReferencias"
bash "$SCRIPTS_DIR/05_salmon_index.sh"    "$rutaReferencias"
bash "$SCRIPTS_DIR/06_snpeff_db.sh"       "$rutaReferencias"
bash "$SCRIPTS_DIR/07_clinvar.sh"         "$rutaReferencias"

echo ""
echo "=========================================="
echo " Todas las referencias instaladas."
echo "=========================================="
