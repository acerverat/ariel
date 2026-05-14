# ARIEL
### Artificial intelligence Rnaseq IdEntification of Leukemia subtypes

> [Versión en español](README.md)

Nextflow pipeline for gene fusion detection and variant analysis in RNA-seq samples from Acute Lymphoblastic Leukemia (ALL).

## Description

ARIEL integrates multiple bioinformatics tools for the detection and analysis of gene fusions:

| Tool | Function |
|------|----------|
| **STAR** | Read alignment against a reference genome |
| **Arriba** | Fusion detection from BAM files |
| **RNApeg** | Splice junction generation for Cicero |
| **Cicero** | Fusion detection from BAM and junctions |
| **FusionCatcher** | Fusion detection using its own database |
| **RaScALL** | Detection of fusions, SNVs, IGH fusions, focal deletions, and DUX4 (ALL-specialized) |
| **Fungi** | Consensus of fusions detected by the methods above |
| **Salmon** | Gene expression quantification |
| **ExprClusters** | Sample clustering by expression (CRLF2, DUX4) |
| **FusionSummary** | Final integrated results report |

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 22.x
- [Docker](https://www.docker.com/)

## Installation

### 1. Build the Docker images

```bash
cd docker

# Main image: STAR, Salmon, Arriba, FusionCatcher, Fungi
bash build_docker.sh

# RaScALL image
bash buildRascall.sh
```

### 2. Prepare reference files

The `generaReferencias.sh` script downloads and builds the required indices in the specified directory:

```bash
bash scripts/generaReferencias.sh /path/to/references
```

This creates `GRCh38_no_alt/` with:

```
GRCh38_no_alt/
├── cicero_references/        # Cicero and RNApeg references
├── fusioncatcher_db/         # FusionCatcher database
├── salmon_index/             # Salmon index (refgenie, partial SA)
├── gencode.v42.annotation.gtf
├── gencode.v42.annotation.gff3
└── STAR_2.7.10b_index/       # STAR index
```

### 3. Generate the configuration

```bash
bash scripts/generaNextflowConfig.sh
```

The script interactively prompts for paths and parameters, and generates a `nextflow.config` file.

### SampleSheet

A TSV file with columns `sample`, `r1`, `r2`:

```
sample	r1	r2
CA001	/path/to/CA001_R1.fastq.gz	/path/to/CA001_R2.fastq.gz
CA002	/path/to/CA002_R1.fastq.gz	/path/to/CA002_R2.fastq.gz
```

## Usage

```bash
nextflow run main.nf -c nextflow.config
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `resultsDir` | Output directory | — |
| `runSampleSheet` | Path to the SampleSheet (TSV) | — |
| `referenceDir` | References directory | — |
| `exprDir` | Reference TPM table (gene panel) | — |
| `ensg_enst_table` | ENSG/ENST conversion table | — |
| `workDir` | Nextflow work directory | — |
| `method_counts` | Minimum number of methods required to call a fusion | `10` |
| `supporting_reads` | Minimum supporting reads for a fusion | `3` |
| `threadsSTAR` | Threads for STAR | `4` |
| `threadsRascall` | Threads for RaScALL | `4` |

## Output structure

```
resultsDir/
├── alignments/          # BAM and BAI files (STAR)
├── junction/            # Junction files (RNApeg)
├── fusions/
│   ├── arriba/          # Fusions detected by Arriba
│   ├── fusioncatcher/   # Fusions detected by FusionCatcher
│   └── cicero/          # Fusions detected by Cicero
├── fungi/               # Fusion consensus (Fungi)
├── rascall/             # RaScALL results
├── quantification/      # Expression quantification (Salmon)
├── ExprClusters/        # Expression clustering results
└── reports/
    ├── report_summary.tsv    # Fusions, subtypes, and breakpoints per sample
    └── otros_hallazgos.tsv   # SNVs, focal deletions, CRLF2/DUX4 expression
```

## Authors

- Alejandra Cervera
- Yun Hernández
- Erubiel Castillo

## Project structure

```
ARIEL/
├── config/
│   └── nextflow.config          # Configuration template
├── docker/
│   ├── Dockerfile               # Main image (ariel-env)
│   ├── Rascall_Dockerfile       # RaScALL image (rascall:1.0)
│   ├── build_docker.sh          # Builds ariel-env
│   └── buildRascall.sh          # Builds rascall:1.0
├── modules/
│   └── modules.nf               # Pipeline modules
└── scripts/
    ├── generaNextflowConfig.sh  # Interactive nextflow.config generator
    └── generaReferencias.sh     # Downloads and builds reference files
```
