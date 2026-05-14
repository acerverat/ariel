#!/bin/bash

# ==========================================
# Generador de archivo de configuración
# ==========================================

echo "=========================================="
echo " Configuración del pipeline de Nextflow"
echo "=========================================="
echo ""

# ==========================================
# Función para validar campos vacíos
# ==========================================

validate_input() {

    local value="$1"
    local field="$2"
    local example="$3"

    if [[ -z "$value" ]]; then
        echo ""
        echo "ERROR: El campo '${field}' es obligatorio."
        echo ""
        echo "Ejemplo:"
        echo "$example"
        echo ""
        exit 1
    fi
}

# =========================
# Parámetros obligatorios
# =========================

read -p "Directorio de resultados: " resultsDir

validate_input \
"$resultsDir" \
"resultsDir" \
"/home/user/resultsAAA"

# ------------------------------------------

read -p "Directorio de reportes de calidad: " reportsDir

validate_input \
"$reportsDir" \
"reportsDir" \
"/home/user/reportsAAA"

# ------------------------------------------

read -p "Ruta del SampleSheet TSV: " runSampleSheet

validate_input \
"$runSampleSheet" \
"runSampleSheet" \
"/home/user/sampleSheetAAA"

# ------------------------------------------

read -p "Directorio de referencias: " referenceDir

validate_input \
"$referenceDir" \
"referenceDir" \
"/home/user/reference"

# ------------------------------------------

read -p "Ruta de la tabla TPMs: " exprDir

validate_input \
"$exprDir" \
"exprDir" \
"/home/user/additionalTPMs.tsv"

# ------------------------------------------

read -p "Ruta de la tabla ENSG/ENST: " ensg_enst_table

validate_input \
"$ensg_enst_table" \
"ensg_enst_table" \
"/home/user/ensg_enst_conversion.tsv"

# -----------------------------------------

read -p "Ruta del directorio de trabajo de nextflow: " workDir

validate_input \
"$workDir" \
"workDir" \
"/home/user/directorio_de_trabajoAAA"


echo ""
echo "=========================================="
echo " Parametros default"
echo "=========================================="
echo ""

# =========================
# method_counts
# =========================

default_method_counts=10

read -p "¿Deseas usar method_counts por default (${default_method_counts})? [Y/n]: " use_default

if [[ "$use_default" =~ ^[Nn]$ ]]; then

    read -p "Nuevo valor para method_counts: " method_counts

    validate_input \
    "$method_counts" \
    "method_counts" \
    "10"

else
    method_counts=$default_method_counts
fi

# =========================
# supporting_reads
# =========================

default_supporting_reads=3

read -p "¿Deseas usar supporting_reads por default (${default_supporting_reads})? [Y/n]: " use_default

if [[ "$use_default" =~ ^[Nn]$ ]]; then

    read -p "Nuevo valor para supporting_reads: " supporting_reads

    validate_input \
    "$supporting_reads" \
    "supporting_reads" \
    "3"

else
    supporting_reads=$default_supporting_reads
fi

# =========================
# threadsSTAR
# =========================

default_threadsSTAR=4

read -p "¿Deseas usar threadsSTAR por default (${default_threadsSTAR})? [Y/n]: " use_default

if [[ "$use_default" =~ ^[Nn]$ ]]; then

    read -p "Nuevo valor para threadsSTAR: " threadsSTAR

    validate_input \
    "$threadsSTAR" \
    "threadsSTAR" \
    "4"

else
    threadsSTAR=$default_threadsSTAR
fi

# =========================
# threadsRascall
# =========================

default_threadsRascall=4

read -p "¿Deseas usar threadsRascall por default (${default_threadsRascall})? [Y/n]: " use_default

if [[ "$use_default" =~ ^[Nn]$ ]]; then

    read -p "Nuevo valor para threadsRascall: " threadsRascall

    validate_input \
    "$threadsRascall" \
    "threadsRascall" \
    "4"

else
    threadsRascall=$default_threadsRascall
fi

# ==========================================
# Genera archivo nextflow.config
# ==========================================

cat > nextflow.config <<EOF
params.resultsDir="${resultsDir}"
params.reportsDir="${reportsDir}"

// tsv con las columnas sample, r1 y r2
params.runSampleSheet="${runSampleSheet}"

// directorio con las referencias de Salmon, STAR, RNApeg, Cicero y Arriba
params.referenceDir="${referenceDir}"

// Tabla de TPMs de muestras de un panel de genes.
params.exprDir="${exprDir}"

// archivo tsv que tiene el gene id y ensembl id
params.ensg_enst_table="${ensg_enst_table}"

// directorio de trabajo de nextflow
workDir="${workDir}"

// parametros fussionSummary
params.method_counts='${method_counts}'
params.supporting_reads='${supporting_reads}'

// hilos
params.threadsSTAR='${threadsSTAR}'
params.threadsRascall='${threadsRascall}'

// Configuracion de los contenedores de Docker cuando Nextflow lo utiliza
docker {
 runOptions = '-u \$(id -u):\$(id -g)'
 enabled = true
 temp = 'auto'
 fixOwnership = true
}


// Configuracion de los modulos
process {
  errorStrategy = 'ignore'
}
EOF

echo ""
echo "=========================================="
echo " Archivo nextflow.config generado"
echo "=========================================="
