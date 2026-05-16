#!/usr/bin/env Rscript

# Agrupa muestras por expresion de un gen usando k-means y genera boxplots.
#
# Args:
#   1. Sample sheet TSV (columna Sample)
#   2. Tabla de TPM (tpm.tsv, salida de preExprCluster.R)
#   3. Gen a clusterizar (e.g. CRLF2)
#   4. Numero de clusters k

library(dplyr)
library(ggplot2)
library(vroom)

args <- commandArgs(trailingOnly = TRUE)

sample_sheet <- read.table(args[1], header = TRUE)
gene <- args[3]
k    <- as.integer(args[4])

expr_dir <- read.table(args[2], sep = "\t", header = TRUE, row.names = 1)
tpm <- expr_dir
tpm[is.na(tpm)] <- 0

log2tpm <- as.data.frame(log(as.matrix(tpm) + 1))
write.table(data.frame(Genes = rownames(log2tpm), log2tpm),
            "expr_log2tpm.tsv", sep = "\t", row.names = FALSE)
write.table(data.frame(Genes = rownames(tpm), tpm),
            "expr_tpm.tsv", sep = "\t", row.names = FALSE)

mykmeans <- function(mat, gene, k) {
  mydf   <- t(mat[gene, ])
  if (sum(mydf) < 1) return(NULL)
  km.res <- kmeans(mydf, k, nstart = 25)

  cl <- data.frame()
  for (i in 1:k) {
    cl <- rbind(cl, data.frame(
      sample  = names(km.res$cluster[which(km.res$cluster == order(km.res$centers)[i])]),
      cluster = paste0("c", i)
    ))
  }

  new_df <- merge(mydf, cl, by.x = 0, by.y = 1)
  colnames(new_df)[1] <- "sample"
  write.table(new_df,
              paste0("log2tpm_", gene, "_", k, "_clusters.tsv"),
              sep = "\t", row.names = FALSE)
  return(new_df)
}

kmeans_df <- mykmeans(log2tpm, gene, k)

if (!is.null(kmeans_df)) {
  plot_df <- kmeans_df
  plot_df$plotname <- kmeans_df$sample
  plot_df <- plot_df %>%
    mutate(plotname = ifelse(plotname %in% sample_sheet$Sample, plotname, ""))

  p <- ggplot(plot_df,
              aes(x = cluster, y = !!sym(gene), color = cluster)) +
    geom_boxplot(outlier.shape = NA)

  p1 <- p +
    ggtitle(paste0(gene, " Expression groups in B-ALL patients")) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank()) +
    ylab("Log2 TPM Expression") + xlab("Expression groups") +
    scale_color_discrete(name = "Events")

  p2 <- p1 + geom_jitter(shape = 16, position = position_jitter(0.2))
  ggsave(paste0("boxplot_", gene, ".png"), plot = p2, device = png)

  p3 <- p2 + geom_text(aes(label = plotname), size = 4)
  ggsave(paste0("boxplot_nombres_", gene, ".png"), plot = p3, device = png)
}
