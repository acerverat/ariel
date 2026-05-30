#!/usr/bin/env Rscript
# parse_vcf_freebayes.R
# Parsea VCFv4.2 (FreeBayes v1.3.10 + SnpEff ANN / GRCh38)
# Extrae transcriptos MANE Select con gen, mutacion (NM_:HGVS.c), proteina (HGVS.p), VAF%
# Solo variantes con clasificacion ClinVar reconocida
# Uso: parse_vcf_freebayes.R <vcf_file> <mane_select_tsv> <output_tsv>

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Uso: parse_vcf_freebayes.R <vcf> <mane_select_tsv> <output_tsv>")

vcf_file    <- args[1]
mane_file   <- args[2]
output_file <- args[3]


# Carga los IDs MANE Select desde el archivo de referencia pre-descargado
load_mane_select <- function(mane_file) {
  mane_raw <- read.table(
    mane_file, sep = "\t", header = TRUE,
    quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE
  )
  colnames(mane_raw)[1] <- sub("^#", "", colnames(mane_raw)[1])
  nm_select <- mane_raw$RefSeq_nuc[mane_raw$MANE_status == "MANE Select"]
  sub("\\.\\d+$", "", nm_select)
}


# Parsea el campo ANN y devuelve un data.frame con una fila por transcripto MANE Select
parse_ann <- function(ann_string, alt_alleles, ao_vals, ro, mane_set) {

  allele_to_ao <- setNames(ao_vals, alt_alleles)

  entries <- strsplit(ann_string, ",")[[1]]

  rows <- lapply(entries, function(entry) {
    f <- strsplit(entry, "\\|")[[1]]
    if (length(f) < 11) return(NULL)

    transcript <- f[7]
    if (!startsWith(transcript, "NM_")) return(NULL)

    transcript_base <- sub("\\.\\d+$", "", transcript)
    if (!(transcript_base %in% mane_set)) return(NULL)

    allele <- f[1]
    gene   <- f[4]
    hgvs_c <- f[10]
    hgvs_p <- f[11]

    ao <- allele_to_ao[[allele]]
    if (is.na(ao)) ao <- ao_vals[1]

    vaf <- if ((ao + ro) > 0) round(ao / (ao + ro) * 100, 2) else NA_real_

    data.frame(
      gen      = gene,
      mutacion = paste0(transcript, ":", hgvs_c),
      proteina = hgvs_p,
      VAF_pct  = vaf,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, Filter(Negate(is.null), rows))
}


# Parsea un VCF completo y devuelve un data.frame
parse_vcf <- function(vcf_file, mane_set) {

  sample_name <- sub("_annotated\\.vcf$", "", basename(vcf_file))
  message("Procesando: ", sample_name)

  lines <- readLines(vcf_file, warn = FALSE)

  header_idx <- grep("^#CHROM", lines)[1]
  col_names  <- strsplit(sub("^#", "", lines[header_idx]), "\t")[[1]]
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[nzchar(data_lines)]

  if (!length(data_lines)) {
    message("  Sin variantes.")
    return(NULL)
  }

  vcf <- read.table(
    text             = paste(data_lines, collapse = "\n"),
    sep              = "\t",
    header           = FALSE,
    quote            = "",
    comment.char     = "",
    stringsAsFactors = FALSE,
    col.names        = col_names
  )

  sample_col <- tail(col_names, 1)

  clnsig_allowed <- c(
    "Conflicting_classifications_of_pathogenicity",
    "Likely_pathogenic",
    "Pathogenic",
    "Uncertain_significance",
    "association",
    "drug_response",
    "risk_factor"
  )

  all_rows <- lapply(seq_len(nrow(vcf)), function(i) {

    row <- vcf[i, ]

    fmt_keys   <- strsplit(row$FORMAT, ":")[[1]]
    fmt_vals   <- strsplit(row[[sample_col]], ":")[[1]]
    fmt        <- setNames(fmt_vals, fmt_keys)
    alleles_gt <- strsplit(fmt["GT"], "[/|]")[[1]]
    if (all(alleles_gt %in% c("0", "."))) return(NULL)

    clnsig_match <- regmatches(
      row$INFO,
      regexpr("(?<=CLNSIG=)[^;]+", row$INFO, perl = TRUE)
    )
    if (!length(clnsig_match) || !(clnsig_match %in% clnsig_allowed)) return(NULL)
    clnsig_val <- clnsig_match

    ro     <- suppressWarnings(as.numeric(fmt["RO"]))
    ao_raw <- fmt["AO"]
    if (is.na(ro) || is.na(ao_raw)) return(NULL)

    ao_vals <- suppressWarnings(as.numeric(strsplit(ao_raw, ",")[[1]]))
    if (any(is.na(ao_vals))) return(NULL)

    ann_match <- regmatches(
      row$INFO,
      regexpr("(?<=ANN=)[^;]+", row$INFO, perl = TRUE)
    )
    if (!length(ann_match) || !nzchar(ann_match)) return(NULL)

    alt_alleles <- strsplit(as.character(row$ALT), ",")[[1]]

    df <- parse_ann(ann_match, alt_alleles, ao_vals, ro, mane_set)
    if (is.null(df) || !nrow(df)) return(NULL)

    df$muestra       <- sample_name
    df$clasificacion <- clnsig_val
    df
  })

  do.call(rbind, Filter(Negate(is.null), all_rows))
}


# ── Procesamiento ──────────────────────────────────────────────────────────────
mane_set <- load_mane_select(mane_file)
result   <- parse_vcf(vcf_file, mane_set)

empty_df <- data.frame(
  muestra = character(), gen = character(), mutacion = character(),
  proteina = character(), VAF_pct = numeric(), clasificacion = character()
)

if (is.null(result) || !nrow(result)) {
  message("No se encontraron variantes MANE Select con clasificacion ClinVar.")
  write.table(empty_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  result <- unique(result[, c("muestra", "gen", "mutacion", "proteina", "VAF_pct", "clasificacion")])
  write.table(result, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message("\nListo. ", nrow(result), " filas exportadas a:\n  ", output_file)
}
