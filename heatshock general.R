library(tidyverse)
library(corrr)
library(ggcorrplot)
library(factoextra)
library(edgeR)
library(rtracklayer)

data <- readr::read_delim("complex_yeast_heatshock.tsv", delim="\t", col_names = T)
gtf <- rtracklayer::import("Saccharomyces_cerevisiae.R64-1-1.75.gtf") %>%
  as.data.frame() %>%
  filter(type == "gene")
gene_length <- gtf[order(match(gtf$gene_id, data$gene_id)), ] %>%
  select(gene_id, width) %>%
  rename(width = "Length")

rpkm <- data %>%
  select(-gene_id) %>%
  DGEList(genes=gene_length) %>%
  calcNormFactors() %>%
  rpkm()

# Cluster hierarchically according to columns
cor <- cor(rpkm)
as.dist(1 - cor) %>%
  hclust(method="complete") %>%
  plot()

