#!/usr/bin/env Rscript

# Recolecta la cuantificacion de Salmon y genera una matriz de TPM
# combinada con la matriz de referencia (panel de genes).
#
# Args:
#   1. Rutas de archivos .quant.sf separadas por coma
#   2. geneId_transcriptId_geneName.tsv  (gene_id, transcript_id, gene_name)
#   3. Tabla de TPM de referencia (panel de genes, con columna ensembl_id y gene_names)

library(tximport)
library(readr)
library(dplyr)
library(vroom)

args <- commandArgs(trailingOnly = TRUE)

files  <- unlist(strsplit(args[1], ","))
nombres <- gsub(".*/(CA[0-9]+)\\.quant\\.sf", "\\1", files)
names(files) <- nombres

archivoRef <- read.table(args[2], sep = "\t", header = TRUE)
tx2gene    <- archivoRef[, c("transcript_id", "gene_id")]

txi <- tximport(files,
                type          = "salmon",
                tx2gene       = tx2gene,
                ignoreTxVersion = FALSE,
                dropInfReps   = TRUE)

table.out <- txi$abundance

write.table(table.out,
            file      = "txiAbundance.tsv",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = TRUE)

geneId_geneName  <- unique(archivoRef[, c("gene_id", "gene_name")])
table.out.names  <- merge(geneId_geneName, table.out, by.x = "gene_id", by.y = 0)

write.table(table.out.names,
            file      = "exprTable.tsv",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE)

tablaGenes          <- table.out.names
tablaGenes$gene_id  <- sub("\\..*", "", tablaGenes$gene_id)
tablaGenes$gene_name <- NULL
colnames(tablaGenes) <- sub("\\..*", "", colnames(tablaGenes))

tabla92 <- read.table(args[3], sep = "\t", header = TRUE)

tablasMerged <- merge(tabla92, tablaGenes, by.x = "ensembl_id", by.y = "gene_id", all.x = TRUE)

write.table(tablasMerged,
            file      = "merged_exprTable.tsv",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE)

tablasMergedSum <- tablasMerged %>%
  group_by(ensembl_id) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
            across(where(~ !is.numeric(.)), ~ first(.x)))

tablasMergedSum <- tablasMergedSum[, c("ensembl_id", "gene_names",
                                       setdiff(colnames(tablasMergedSum),
                                               c("ensembl_id", "gene_names")))]

# resolve merge duplicates: keep .y columns, drop .x/.y suffixes
cols_y    <- grep("\\.y$", colnames(tablasMergedSum), value = TRUE)
cols_base <- sub("\\.y$", "", cols_y)
for (col in cols_base) {
  tablasMergedSum[[col]] <- tablasMergedSum[[paste0(col, ".y")]]
}
tablasMerged <- tablasMergedSum[, !grepl("\\.x$|\\.y$", colnames(tablasMergedSum))]
tablasMerged <- na.omit(tablasMerged)
tablasMerged$ensembl_id <- NULL

write.table(tablasMerged,
            file      = "tpm.tsv",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE)
