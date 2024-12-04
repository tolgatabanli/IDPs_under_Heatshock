library(DESeq2)
library(tidyverse)
library(SummarizedExperiment)
library(edgeR)

# You should transform to FPKM or RPKM

colData <- read_tsv("sample_mapping.tsv", col_names = T) %>%
  mutate(across(everything(), as.factor)) %>%
  filter(temperature != 25) %>%
  column_to_rownames(var = "sample")

heatshock_counts <- read_tsv(here::here("complex_yeast_heatshock.tsv"), col_names = T) %>%
  column_to_rownames(var = "gene_id") %>%
  relocate(rownames(colData)) %>%
  select(-starts_with("25")) %>%
  as.matrix()

gtf <- rtracklayer::import('Saccharomyces_cerevisiae.R64-1-1.75.gtf') %>%
  as.data.frame() %>%
  filter(type == "transcript") %>%
  column_to_rownames(var = "gene_id") %>%
  '['(rownames(heatshock_counts), ) %>%
  GRanges()

# RPKM
rpkm <- heatshock_counts %>%
  DGEList(genes=data.frame(Length=gtf@ranges@width)) %>%
  calcNormFactors() %>%
  rpkm()
# But now they are not integers

se <- SummarizedExperiment(assays=list(counts=heatshock_counts),
                     rowRanges=gtf, colData=colData)

# DESeq
dds <- DESeqDataSet(se, design = ~ temperature + time + knockout)
dds$temperature = relevel(deseq$temperature, ref = "37")
deseq <- DESeq(deseq)

# Results, Contrastin and Plotting
res <- results(deseq)
res <- results(res, contrast = c("temperature", "37", "42"))
summary(res)

DESeq2::plotMA(res)

# LFC ?
lfc <- lfcShrink(res, coef = "temperature_42_vs_37", type = "apeglm")
plotMA(lfc)
