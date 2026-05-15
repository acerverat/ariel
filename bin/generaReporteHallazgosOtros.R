#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
# args[1]: fusiones_otras.csv    (salida de generaReporteHallazgosPrincipales.R)
# args[2]: rascall_data.csv      (salida de generaReporteHallazgosPrincipales.R)
# args[3]: directorio de Cicero  (con archivos _all_ITD.txt)
# args[4]: sample sheet (TSV)
# args[5]: log2tpm_CRLF2_3_clusters.tsv

###### Read inputs

message("Reading ", args[1])
fusiones_otras <- read.table(args[1], header = TRUE, sep = "\t", na.strings = "NA")

message("Reading ", args[2])
rascall_data <- read.table(args[2], header = TRUE, sep = "\t", na.strings = "NA")

message("Reading ", args[4])
all_samples <- read.table(args[4], header = TRUE, sep = "\t") |>
  select(Muestra = Sample)

message("Reading ", args[5])
crlf2_expr <- read.table(args[5], header = TRUE, sep = "\t") |>
  select(Muestra = sample, CRLF2_expr = cluster) |>
  mutate(CRLF2_expr = recode(CRLF2_expr, c1 = "Bajo-c1", c2 = "Medio-c2", c3 = "Alto-c3"))

message("Reading cicero ITD files from: ", args[3])
cicero_itd_files <- list.files(args[3], pattern = "_all_ITD\\.txt$", full.names = TRUE)

hallazgos_crlf2 <- lapply(cicero_itd_files, \(f) {
  df <- read.table(f, header = TRUE, sep = "\t", quote = "", fill = TRUE)
  if (any(grepl("CRLF2", df$geneA))) {
    data.frame(Muestra = sub("_all_ITD\\.txt$", "", basename(f)),
               Duplicacion_CRLF2 = "CRLF2 dup")
  }
}) |>
  bind_rows()

###### Filter fusions
# Criteria:
#   1. Supporting_reads >= 10
#   2. If called by cicero, medal > 2
#   3. No NAs in either gene of the fusion name

hallazgos_fusiones <- fusiones_otras |>
  filter(
    # No NA gene names (string "NA" or actual NA)
    !is.na(Fusion),
    !grepl("(^NA--|--NA$|--NA--|^NA$)", Fusion),
    # Supporting reads threshold
    Supporting_reads >= 10,
    # Cicero medal filter: only applies when cicero called the fusion
    !grepl("cicero", Metodos, ignore.case = TRUE) | (Cicero_medal > 2),
    # Remove fusions where any database annotates them as fully overlapping genes
    !grepl("fully_overlapping", coalesce(annotations, ""), ignore.case = TRUE)
  ) |>
  select(-Supporting_reads, -Cicero_medal, -annotations) |>
  mutate(
    # Merge subtype columns: combine when both are informative, otherwise keep whichever has a value
    Subtipo = case_when(
      Subtipo != "-" & Subtipo_Emergente != "-" ~ paste(Subtipo, Subtipo_Emergente, sep = ", "),
      Subtipo != "-"                            ~ Subtipo,
      Subtipo_Emergente != "-"                  ~ Subtipo_Emergente,
      TRUE                                      ~ "-"
    ),
    # Single SR column: max across all three callers, ignoring missing ("-")
    SR = pmax(
      as.integer(na_if(SR_Arriba,        "-")),
      as.integer(na_if(SR_Cicero,        "-")),
      as.integer(na_if(SR_Fusioncatcher, "-")),
      na.rm = TRUE
    )
  ) |>
  select(-Subtipo_Emergente, -Punto_de_corte, -SR_Arriba, -SR_Cicero, -SR_Fusioncatcher)

###### RaScALL non-fusion rows (one row per sample, alterations collapsed)

hallazgos_snv <- rascall_data |>
  filter(Target_Type == "SNV", !is.na(Alteration)) |>
  group_by(Muestra = File) |>
  summarise(SNV_RaScALL = paste(Alteration, collapse = "; "), .groups = "drop")

hallazgos_deleciones <- rascall_data |>
  filter(Target_Type == "focalDeletion", !is.na(Alteration)) |>
  group_by(Muestra = File) |>
  summarise(Deleciones_Focales_Rascall = paste(Alteration, collapse = "; "), .groups = "drop")

hallazgos_dux4 <- rascall_data |>
  filter(Target_Type == "DUX4", !is.na(Alteration)) |>
  group_by(Muestra = File) |>
  summarise(DUX4r_Rascall = paste(Alteration, collapse = "; "), .groups = "drop")

###### Combine all hallazgos
# For samples already present in fusions: join annotation columns onto existing rows.
# For samples with no fusions at all: add a single new row with all available annotations.

rascall_annotations <- list(hallazgos_snv, hallazgos_deleciones, hallazgos_dux4, hallazgos_crlf2) |>
  reduce(full_join, by = "Muestra")

# All-sample annotations: CRLF2_expr for every sample in the cohort
all_sample_annotations <- all_samples |>
  left_join(crlf2_expr, by = "Muestra") |>
  left_join(rascall_annotations, by = "Muestra")

hallazgos_otros <- hallazgos_fusiones |>
  left_join(all_sample_annotations, by = "Muestra") |>
  bind_rows(
    anti_join(all_sample_annotations, hallazgos_fusiones, by = "Muestra")
  )

###### Write output

hallazgos_otros <- hallazgos_otros |>
  mutate(across(everything(), ~replace_na(as.character(.), "-")))

write.table(hallazgos_otros, "hallazgos_otros.csv", sep = "\t", row.names = FALSE, quote = FALSE)

