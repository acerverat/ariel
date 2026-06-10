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
| **FreeBayes** | Llamado de variantes (SNVs e indels) |
| **SnpEff + SnpSift** | Anotación funcional de variantes y filtrado con ClinVar |
| **FusionSummary** | Reporte final con resultados integrados |

## Requisitos

- [Nextflow](https://www.nextflow.io/) >= 22.x
- [Docker](https://www.docker.com/)

## Instalación

### 1. Obtener las imágenes de Docker

Las imágenes están disponibles en Docker Hub y se descargan automáticamente al ejecutar el pipeline. También se pueden descargar manualmente:

```bash
docker pull acerverat/ariel-env:latest
docker pull acerverat/rascall:1.0
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
| `06_snpeff_db.sh` | Base de datos de SnpEff (GRCh38.p14) |
| `07_clinvar.sh` | VCF de ClinVar (GRCh38) y tabla MANE Select para anotación de variantes |

Esto crea `GRCh38_no_alt/` con:

```
GRCh38_no_alt/
├── cicero_references/               # Referencias de Cicero y RNApeg
├── fusioncatcher_db/                # Base de datos de FusionCatcher
├── salmon_index/                    # Índice de Salmon
├── snpeff_db/                       # Base de datos de SnpEff, ClinVar y MANE Select
├── gencode.v42.annotation.gtf
├── gencode.v42.annotation.gff3
├── geneId_transcriptId_geneName.tsv # Tabla gene_id / transcript_id / gene_name
└── STAR_2.7.10b_index/              # Índice de STAR
```

### 3. Configurar los parámetros

Copia `run_local.sh` o `run_remote.sh` a tu directorio de trabajo, completa las rutas marcadas con `<...>` y ejecútalo. El script genera automáticamente el archivo de parámetros y lanza el pipeline.

### SampleSheet

Archivo TSV con las columnas `Sample`, `R1`, `R2`:

```
Sample	R1	R2
CA001	/ruta/CA001_R1.fastq.gz	/ruta/CA001_R2.fastq.gz
CA002	/ruta/CA002_R1.fastq.gz	/ruta/CA002_R2.fastq.gz
```

## Uso

```bash
# Ejecución local (desarrollo y pruebas)
bash run_local.sh

# Ejecución desde una versión publicada en GitHub
bash run_remote.sh
```

## Parámetros

| Parámetro | Descripción | Default |
|-----------|-------------|---------|
| `resultsDir` | Directorio de resultados | — |
| `runSampleSheet` | Ruta del SampleSheet (TSV) | — |
| `referenceDir` | Directorio de referencias (`GRCh38_no_alt/`) | — |
| `threadsSTAR` | Hilos para STAR | `4` |
| `threadsRascall` | Hilos para RaScALL | `4` |
| `threadsFungi` | Hilos para Fungi | `15` |
| `fungiFilterEnsembl` | Filtros Ensembl de Fungi | `"same_gene,homologs"` |
| `fungiFilterDb` | Filtros de base de datos de Fungi | `"banned,paralog"` |
| `fungiFilterMinCount` | Conteo mínimo de lecturas para Fungi | `0` |

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
├── variants/
│   ├── freebayes/       # VCFs de FreeBayes
│   └── snpeff/          # VCFs anotados por SnpEff+SnpSift
├── qc/
│   ├── beforeTrimm/     # MultiQC antes del filtrado
│   ├── afterTrimm/      # FastQC + Fastp + MultiQC después del filtrado
│   └── trimmed/         # Lecturas filtradas (Fastp)
└── reportes/
    ├── hallazgos_principales.csv  # Fusiones principales, subtipos y breakpoints por muestra
    └── variantes_NM.tsv           # Variantes MANE Select con clasificación ClinVar
```

## Autores

- Alejandra Cervera
- Yun Hernández
- Erubiel Castillo

## Estructura del proyecto

```
ARIEL/
├── bin/
│   ├── generaReporteHallazgosPrincipales.R # Reporte de fusiones principales e integración con RaScALL
│   ├── generaReporteHallazgosOtros.R       # Reporte de fusiones secundarias, SNVs y deleciones
│   └── parse_vcf_freebayes.R               # Reporte de variantes MANE Select con clasificación ClinVar
├── config/
│   └── params.yaml                         # Referencia de parámetros disponibles
├── docker/
│   ├── Dockerfile                          # Imagen principal (ariel-env)
│   └── build_ariel.sh                      # Construye ariel-env
├── modules/
│   └── modules.nf                          # Módulos del pipeline
├── scripts/
│   ├── generaReferencias.sh                # Script maestro de referencias (pasos 01–07)
│   └── references/
│       ├── 01_cicero.sh
│       ├── 02_gencode.sh
│       ├── 03_star_index.sh
│       ├── 04_fusioncatcher_db.sh
│       ├── 05_salmon_index.sh
│       ├── 06_snpeff_db.sh
│       └── 07_clinvar.sh
├── main.nf                                 # Flujo de trabajo principal
├── run_local.sh                            # Ejecución local (desarrollo y pruebas)
└── run_remote.sh                           # Ejecución desde versión publicada en GitHub
```
