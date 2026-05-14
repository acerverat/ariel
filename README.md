# ARIEL
### Análisis de RNA-seq para IdEntificación de subtipos de Leucemia

Pipeline de Nextflow para la detección de fusiones génicas y variantes en muestras de RNA-seq de Leucemia Linfoblástica Aguda (LLA).

## Descripción

ARIEL integra múltiples herramientas para la detección y análisis de fusiones génicas:

| Herramienta | Función |
|-------------|---------|
| **STAR** | Alineamiento de lecturas contra genoma de referencia |
| **Arriba** | Detección de fusiones a partir de BAM |
| **RNApeg** | Generación de junctions para Cicero |
| **Cicero** | Detección de fusiones a partir de BAM y junctions |
| **FusionCatcher** | Detección de fusiones con base de datos propia |
| **RaScALL** | Detección de fusiones, SNVs, fusiones IGH, deleciones focales y DUX4 (especializado en LLA) |
| **Fungi** | Consenso de fusiones detectadas por los métodos anteriores |
| **Salmon** | Cuantificación de expresión génica |
| **ExprClusters** | Agrupamiento de muestras por expresión (CRLF2, DUX4) |
| **FusionSummary** | Reporte final con resultados integrados |

## Requisitos

- [Nextflow](https://www.nextflow.io/) >= 22.x
- [Docker](https://www.docker.com/)

## Instalación

### 1. Construir las imágenes de Docker

```bash
cd docker

# Imagen principal: STAR, Salmon, Arriba, FusionCatcher, Fungi
bash build_docker.sh

# Imagen de RaScALL
bash buildRascall.sh
```

### 2. Preparar las referencias

El script `generaReferencias.sh` descarga y construye los índices necesarios en el directorio indicado:

```bash
bash scripts/generaReferencias.sh /ruta/al/directorio/referencias
```

Esto crea `GRCh38_no_alt/` con:

```
GRCh38_no_alt/
├── cicero_references/        # Referencias de Cicero y RNApeg
├── gencode.v42.annotation.gtf
├── gencode.v42.annotation.gff3
└── STAR_2.7.10b_index/       # Índice de STAR
```

> El índice de Salmon (`salmon_index/`) debe generarse por separado y colocarse dentro de `GRCh38_no_alt/`.

### 3. Generar la configuración

```bash
bash scripts/generaNextflowConfig.sh
```

El script solicita de forma interactiva las rutas y parámetros, y genera un archivo `nextflow.config`.

### SampleSheet

Archivo TSV con las columnas `sample`, `r1`, `r2`:

```
sample	r1	r2
CA001	/ruta/CA001_R1.fastq.gz	/ruta/CA001_R2.fastq.gz
CA002	/ruta/CA002_R1.fastq.gz	/ruta/CA002_R2.fastq.gz
```

## Uso

```bash
nextflow run main.nf -c nextflow.config
```

## Parámetros

| Parámetro | Descripción | Default |
|-----------|-------------|---------|
| `resultsDir` | Directorio de resultados | — |
| `runSampleSheet` | Ruta del SampleSheet (TSV) | — |
| `referenceDir` | Directorio de referencias | — |
| `exprDir` | Tabla de TPMs de referencia (panel de genes) | — |
| `ensg_enst_table` | Tabla de conversión ENSG/ENST | — |
| `workDir` | Directorio de trabajo de Nextflow | — |
| `method_counts` | Métodos mínimos que deben detectar una fusión | `10` |
| `supporting_reads` | Lecturas de soporte mínimas | `3` |
| `threadsSTAR` | Hilos para STAR | `4` |
| `threadsRascall` | Hilos para RaScALL | `4` |

## Estructura de resultados

```
resultsDir/
├── alignments/          # BAM y BAI (STAR)
├── junction/            # Junctions (RNApeg)
├── fusions/
│   ├── arriba/          # Fusiones detectadas por Arriba
│   ├── fusioncatcher/   # Fusiones detectadas por FusionCatcher
│   └── cicero/          # Fusiones detectadas por Cicero
├── fungi/               # Consenso de fusiones (Fungi)
├── rascall/             # Resultados de RaScALL
├── quantification/      # Cuantificación de expresión (Salmon)
├── ExprClusters/        # Resultados de clustering de expresión
└── reports/
    ├── report_summary.tsv    # Fusiones, subtipos y breakpoints por muestra
    └── otros_hallazgos.tsv   # SNVs, deleciones focales, expresión CRLF2/DUX4
```

## Estructura del proyecto

```
ARIEL/
├── config/
│   └── nextflow.config          # Plantilla de configuración
├── docker/
│   ├── Dockerfile               # Imagen principal (ball_classifier_pruebas)
│   ├── Rascall_Dockerfile       # Imagen de RaScALL (rascall:1.0)
│   ├── build_docker.sh          # Construye ball_classifier_pruebas
│   └── buildRascall.sh          # Construye rascall:1.0
├── modules/
│   └── modules.nf               # Módulos del pipeline
└── scripts/
    ├── generaNextflowConfig.sh  # Generador interactivo de nextflow.config
    └── generaReferencias.sh     # Descarga y construye referencias
```
