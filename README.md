# ARIEL
### Análisis de RNA-seq para IdEntificación de subtipos de Leucemia

> [English version](README.en.md)

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

`generaReferencias.sh` es el script maestro que ejecuta en orden los subscripts de `scripts/references/`. Cada subscript verifica si sus referencias ya están instaladas antes de descargar, por lo que es seguro reejecutar sin reinstalar lo que ya existe. Si un paso falla, se puede reejecutar solo ese subscript.

```bash
# Instalar todas las referencias
bash scripts/generaReferencias.sh /ruta/al/directorio/referencias

# O ejecutar un paso individual, por ejemplo:
bash scripts/references/03_star_index.sh /ruta/al/directorio/referencias
```

| Script | Contenido |
|--------|-----------|
| `01_cicero.sh` | Referencias de Cicero y RNApeg |
| `02_gencode.sh` | Anotaciones Gencode v42 (GTF y GFF3) |
| `03_star_index.sh` | Índice de STAR (requiere 01 y 02) |
| `04_fusioncatcher_db.sh` | Base de datos de FusionCatcher |
| `05_salmon_index.sh` | Índice de Salmon (desde transcriptoma Gencode v42) y tabla `geneId_transcriptId_geneName.tsv` |

Esto crea `GRCh38_no_alt/` con:

```
GRCh38_no_alt/
├── cicero_references/               # Referencias de Cicero y RNApeg
├── fusioncatcher_db/                # Base de datos de FusionCatcher
├── salmon_index/                    # Índice de Salmon
├── gencode.v42.annotation.gtf
├── gencode.v42.annotation.gff3
├── geneId_transcriptId_geneName.tsv # Tabla gene_id / transcript_id / gene_name
└── STAR_2.7.10b_index/              # Índice de STAR
```

### 3. Generar la configuración

```bash
bash scripts/generaNextflowConfig.sh
```

El script solicita de forma interactiva las rutas y parámetros, y genera un archivo `nextflow.config`.

### SampleSheet

Archivo TSV con las columnas `Sample`, `R1`, `R2`:

```
Sample	R1	R2
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
| `ensg_enst_table` | Tabla `geneId_transcriptId_geneName.tsv` generada por `05_salmon_index.sh` | — |
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
    ├── hallazgos_principales.csv  # Fusiones principales, subtipos y breakpoints por muestra
    └── hallazgos_otros.csv        # Fusiones secundarias, SNVs, deleciones focales, expresión CRLF2/DUX4
```

## Autores

- Alejandra Cervera
- Yun Hernández
- Erubiel Castillo

## Estructura del proyecto

```
ARIEL/
├── bin/
│   ├── preExprCluster.R                    # Genera matriz TPM desde cuantificaciones de Salmon
│   ├── kmeans.R                            # Clustering k-means por gen y generación de boxplots
│   ├── generaReporteHallazgosPrincipales.R # Reporte de fusiones principales e integración con RaScALL
│   └── generaReporteHallazgosOtros.R       # Reporte de fusiones secundarias, SNVs, deleciones y expresión
├── config/
│   └── nextflow.config                     # Plantilla de configuración
├── docker/
│   ├── Dockerfile                          # Imagen principal (ariel-env)
│   ├── Rascall_Dockerfile                  # Imagen de RaScALL (rascall:1.0)
│   ├── build_docker.sh                     # Construye ariel-env
│   └── buildRascall.sh                     # Construye rascall:1.0
├── modules/
│   └── modules.nf                          # Módulos del pipeline
├── scripts/
│   ├── generaNextflowConfig.sh             # Generador interactivo de nextflow.config
│   ├── generaReferencias.sh                # Script maestro de referencias
│   └── references/                         # Subscripts por herramienta
│       ├── 01_cicero.sh
│       ├── 02_gencode.sh
│       ├── 03_star_index.sh
│       ├── 04_fusioncatcher_db.sh
│       └── 05_salmon_index.sh
└── main.nf                                 # Flujo de trabajo principal
```
