#!/usr/bin/env Rscript
library(vroom)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
# args[1]: sample sheet (TSV)
# args[2]: bp_consensus_report.tsv
# args[3]: directorio de resultados de RaScALL

make_canonical <- function(fusion) {
  g1 <- sub("--.*", "", fusion)
  g2 <- sub(".*--", "", fusion)
  paste(pmin(g1, g2), pmax(g1, g2), sep = "--")
}

###### Read inputs

# Sample sheet
message("Reading ", args[1])
samples <- read.table(args[1], header = TRUE, sep = "\t")

# bp_consensus
message("Reading ", args[2])
bp <- read.table(args[2], header = TRUE, sep = "\t")

# RaScALL
message("Reading from folder: ", args[3])
Rascall_Folder <- args[3]

Rascall_file_names <- list.files(Rascall_Folder,
                                 pattern = "_final_variants.csv",
                                 recursive = TRUE,
                                 full.names = FALSE)

Rascall.data <- vroom(
  file.path(Rascall_Folder, Rascall_file_names),
  col_select = c("File", "Target_Type", "Alteration", "Type")
)

# Create fusion table from RaScALL
Fusion.Rascall <- Rascall.data |>
  filter(grepl("Fusion", Target_Type), !is.na(Alteration)) |>
  filter(!grepl("-NA$|^NA-|--NA|NA--", Alteration)) |>
  select(-Target_Type) |>
  rename(Sample = File, FusionName = Alteration) |>
  mutate(
    FusionName = gsub("-", "--", FusionName),
    Lista_Metodos_Fusiones = "RaScALL"
  ) |>
  select(FusionName, Sample, Lista_Metodos_Fusiones)

# Collapse reversed fusion pairs in bp, 
# keeping orientation with highest Supporting_reads
bp_collapsed <- bp |>
  mutate(canonical = make_canonical(FusionName)) |>
  group_by(canonical, Sample) |>
  mutate(reversed_also_present = n() > 1) |>
  slice_max(order_by = Supporting_reads, n = 1, with_ties = FALSE) |>
  slice_max(order_by = Methods_count, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(-canonical)

# Orientation lookup from bp_collapsed
bp_orientation <- bp_collapsed |>
  mutate(canonical = make_canonical(FusionName)) |>
  select(canonical, Sample, FusionName_bp = FusionName)

# Deduplicate Fusion.Rascall, using bp_collapsed orientation where available
rascall_dedup <- Fusion.Rascall |>
  mutate(canonical = make_canonical(FusionName)) |>
  group_by(canonical, Sample) |>
  summarise(
    reversed_also_present_rascall = n() > 1,
    Lista_Metodos_Fusiones = first(Lista_Metodos_Fusiones),
    FusionName_fallback = first(FusionName),
    .groups = "drop"
  ) |>
  left_join(bp_orientation, by = c("canonical", "Sample")) |>
  mutate(FusionName = coalesce(FusionName_bp, FusionName_fallback)) |>
  select(-FusionName_bp, -FusionName_fallback)

emerging_patterns <- c(
  "DUX4"        = "DUX4",
  "MEF2D"       = "MEF2D",
  "ZNF384"      = "ZNF384",
  "ZNF384-like" = "SMARCA2--ZNF362|ZNF362--SMARCA2|TAF15--ZNF362|ZNF362--TAF15",
  "BCL/MYC"     = "BCL6|BCL2|MYC",
  "NUTM1"       = "NUTM1",
  "IKZF1"       = "IKZF1",
  "CRLF2r"      = "DDX3X|USP9X",
  "PAX5r"       = "PAX5",
  "IGH--CEBPE"  = "IGH--CEBPE|CEBPE--IGH"
)

# Integrate bp_collapsed and rascall_dedup
integrated <- full_join(
  rascall_dedup |> rename(FusionName_rascall = FusionName),
  bp_collapsed |>
    mutate(canonical = make_canonical(FusionName)) |>
    rename(FusionName_bp = FusionName),
  by = c("canonical", "Sample")
) |>
  mutate(
    FusionName = coalesce(FusionName_rascall, FusionName_bp),
    reversed_also_present = coalesce(reversed_also_present, FALSE) | coalesce(reversed_also_present_rascall, FALSE),
    Methods_count = case_when(
      !is.na(FusionName_rascall) & !is.na(FusionName_bp) ~ Methods_count + 1L,
      !is.na(FusionName_rascall) & is.na(FusionName_bp)  ~ 1L,
      TRUE ~ Methods_count
    ),
    Methods_list = case_when(
      !is.na(FusionName_rascall) & !is.na(FusionName_bp) ~ paste(Methods_list, Lista_Metodos_Fusiones, sep = ","),
      !is.na(FusionName_rascall) & is.na(FusionName_bp)  ~ Lista_Metodos_Fusiones,
      TRUE ~ Methods_list
    )
  ) |>
  select(-canonical, -FusionName_rascall, -FusionName_bp,
         -Lista_Metodos_Fusiones, -reversed_also_present_rascall) |>
  relocate(FusionName, Sample)

###### Build output table

# Annotate with subtype and SR columns before splitting
hallazgos_annotated <- integrated |>
  mutate(
    Subtipo = case_when(
      grepl("BCR--ABL1|ABL1--BCR", FusionName)                      ~ "Ph",
      grepl("ETV6--RUNX1|RUNX1--ETV6", FusionName)                  ~ "ETV6--RUNX1",
      grepl("TCF3--PBX1|PBX1--TCF3", FusionName)                    ~ "TCF3--PBX1",
      grepl("TCF3--HLF|HLF--TCF3|TCF4--HLF|HLF--TCF4", FusionName) ~ "TCF3--HLF",
      grepl("IGH--IL3|IL3--IGH", FusionName)                        ~ "IGH--IL3",
      grepl("KMT2A", FusionName)                                    ~ "KMT2Ar",
      grepl("CRLF2|JAK2|ABL|EPOR|PDGFRB", FusionName)              ~ "Ph-like",
      grepl("ETV6", FusionName)                                     ~ "ETV6-like",
      .default = "-"
    ),
    Subtipo_Emergente = map_chr(FusionName, \(fn) {
      matched <- names(emerging_patterns)[sapply(emerging_patterns, grepl, x = fn)]
      if (length(matched) == 0) "-" else paste(matched, collapse = ", ")
    }),
    SR_Arriba        = str_extract(additional_info, "(?<=arriba:SR=)\\d+") |> replace_na("-"),
    SR_Cicero        = str_extract(additional_info, "(?<=cicero:SR=)\\d+") |> replace_na("-"),
    SR_Fusioncatcher = str_extract(additional_info, "(?<=fusioncatcher:SR=)\\d+") |> replace_na("-"),
    es_principal = grepl("RaScALL", coalesce(Methods_list, "")) |
                   #Subtipo != "-" |
                   (coalesce(Methods_count, 0L) >= 3 & coalesce(Supporting_reads, 0L) >= 10)
  )

# Helper to apply final column selection and Spanish names
final_cols <- function(df) {
  df |> select(
    Muestra          = Sample,
    Fusion           = FusionName,
    Metodos          = Methods_list,
    Subtipo,
    Subtipo_Emergente,
    Punto_de_corte   = best_bp,
    SR_Arriba,
    SR_Cicero,
    SR_Fusioncatcher
  )
}

hallazgos_principales <- hallazgos_annotated |> filter(es_principal) |> final_cols()

# fusiones_otras keeps extra filtering columns needed by generaReporteHallazgosOtros.R
fusiones_otras <- hallazgos_annotated |>
  filter(!es_principal) |>
  mutate(
    Cicero_medal = str_extract(additional_info, "(?<=medal=)\\d+") |> as.integer()
  ) |>
  select(
    Muestra          = Sample,
    Fusion           = FusionName,
    Metodos          = Methods_list,
    Subtipo,
    Subtipo_Emergente,
    Punto_de_corte   = best_bp,
    SR_Arriba,
    SR_Cicero,
    SR_Fusioncatcher,
    Supporting_reads,
    Cicero_medal,
    annotations
  )

write.table(hallazgos_principales, "hallazgos_principales.csv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(fusiones_otras,        "fusiones_otras.csv",        sep = "\t", row.names = FALSE, quote = FALSE)
write.table(Rascall.data,          "rascall_data.csv",          sep = "\t", row.names = FALSE, quote = FALSE)

