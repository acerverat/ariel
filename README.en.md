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
| **FreeBayes** | Variant calling (SNVs and indels) |
| **SnpEff + SnpSift** | Functional variant annotation and ClinVar filtering |
| **FusionSummary** | Final integrated results report |

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 22.x
- [Docker](https://www.docker.com/)

## Installation

### 1. Get the Docker images

Images are available on Docker Hub and are pulled automatically when the pipeline runs. They can also be pulled manually:

```bash
docker pull acerverat/ariel-env:latest
docker pull acerverat/rascall:1.0
```

### 2. Prepare reference files

`generaReferencias.sh` is the master script that runs the subscripts in `scripts/references/` in order. Each subscript checks whether its references are already installed before downloading, so it is safe to rerun without reinstalling existing data. If a step fails, only that subscript needs to be rerun.

```bash
# Install all references
bash scripts/generaReferencias.sh /path/to/references


# Or run a single step, for example:
bash scripts/references/03_star_index.sh /path/to/references
```

| Script | Content |
|--------|---------|
| `01_cicero.sh` | Cicero and RNApeg references |
| `02_gencode.sh` | Gencode v42 annotations (GTF and GFF3) |
| `03_star_index.sh` | STAR index (requires 01 and 02) |
| `04_fusioncatcher_db.sh` | FusionCatcher database |
| `05_salmon_index.sh` | Salmon index (built from Gencode v42 transcriptome) and `geneId_transcriptId_geneName.tsv` table |
| `06_snpeff_db.sh` | SnpEff database (GRCh38.p14) |
| `07_clinvar.sh` | ClinVar VCF (GRCh38) and MANE Select table for variant annotation |

This creates `GRCh38_no_alt/` with:

```
GRCh38_no_alt/
├── cicero_references/               # Cicero and RNApeg references
├── fusioncatcher_db/                # FusionCatcher database
├── salmon_index/                    # Salmon index
├── snpeff_db/                       # SnpEff database, ClinVar, and MANE Select
├── gencode.v42.annotation.gtf
├── gencode.v42.annotation.gff3
├── geneId_transcriptId_geneName.tsv # gene_id / transcript_id / gene_name table
└── STAR_2.7.10b_index/              # STAR index
```

### 3. Configure parameters

Copy `run_local.sh` or `run_remote.sh` to your working directory, fill in the paths marked with `<...>`, and run it. The script automatically generates the parameters file and launches the pipeline.

### SampleSheet

A TSV file with columns `Sample`, `R1`, `R2`:

```
Sample	R1	R2
CA001	/path/to/CA001_R1.fastq.gz	/path/to/CA001_R2.fastq.gz
CA002	/path/to/CA002_R1.fastq.gz	/path/to/CA002_R2.fastq.gz
```

## Usage

```bash
# Local execution (development and testing)
bash run_local.sh

# Execution from a published GitHub release
bash run_remote.sh
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `resultsDir` | Output directory | — |
| `runSampleSheet` | Path to the SampleSheet (TSV) | — |
| `referenceDir` | References directory (`GRCh38_no_alt/`) | — |
| `threadsSTAR` | Threads for STAR | `4` |
| `threadsRascall` | Threads for RaScALL | `4` |
| `threadsFungi` | Threads for Fungi | `15` |
| `fungiFilterEnsembl` | Fungi Ensembl filters | `"same_gene,homologs"` |
| `fungiFilterDb` | Fungi database filters | `"banned,paralog"` |
| `fungiFilterMinCount` | Minimum read count for Fungi | `0` |

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
├── variants/
│   ├── freebayes/       # FreeBayes VCFs
│   └── snpeff/          # SnpEff+SnpSift annotated VCFs
├── qc/
│   ├── beforeTrimm/     # MultiQC before trimming
│   ├── afterTrimm/      # FastQC + Fastp + MultiQC after trimming
│   └── trimmed/         # Trimmed reads (Fastp)
└── reportes/
    ├── hallazgos_principales.csv  # Main fusions, subtypes, and breakpoints per sample
    └── variantes_NM.tsv           # MANE Select variants with ClinVar classification
```

## Authors

- Alejandra Cervera
- Yun Hernández
- Erubiel Castillo

## Project structure

```
ARIEL/
├── bin/
│   ├── generaReporteHallazgosPrincipales.R # Main fusion report and RaScALL integration
│   ├── generaReporteHallazgosOtros.R       # Secondary fusions, SNVs, and deletions report
│   └── parse_vcf_freebayes.R               # MANE Select variants report with ClinVar classification
├── config/
│   └── params.yaml                         # Reference of available parameters
├── docker/
│   ├── Dockerfile                          # Main image (ariel-env)
│   └── build_ariel.sh                      # Builds ariel-env
├── modules/
│   └── modules.nf                          # Pipeline modules
├── scripts/
│   ├── generaReferencias.sh                # Master reference installation script (steps 01–05)
│   └── references/
│       ├── 01_cicero.sh
│       ├── 02_gencode.sh
│       ├── 03_star_index.sh
│       ├── 04_fusioncatcher_db.sh
│       ├── 05_salmon_index.sh
│       ├── 06_snpeff_db.sh
│       └── 07_clinvar.sh
├── main.nf                                 # Main workflow
├── run_local.sh                            # Local execution (development and testing)
└── run_remote.sh                           # Execution from a published GitHub release
```
