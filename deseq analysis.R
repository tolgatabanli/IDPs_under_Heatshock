library(DESeq2)
library(tidyverse)
library(SummarizedExperiment)
library(edgeR)

# You should transform to FPKM or RPKM

colData <- read_tsv("sample_mapping.tsv", col_names = T) %>%
  mutate(temperature = factor(temperature),
         replicate = factor(replicate)) %>%
  column_to_rownames(var = "sample")

heatshock_counts <- read_tsv(here::here("complex_yeast_heatshock.tsv"), col_names = T) %>%
  column_to_rownames(var = "gene_id") %>%
  relocate(rownames(colData)) %>%
  as.matrix()

gtf <- rtracklayer::import('Saccharomyces_cerevisiae.R64-1-1.75.gtf') %>%
  as.data.frame() %>%
  filter(type == "transcript") %>%
  GRanges()

# TODO
rpkm <- heatshock_counts %>%
  DGEList(genes=data.frame(Length=gtf@ranges@width)) %>%
  calcNormFactors() %>%
  rpkm()
# But now they are not integers, which DESeq requires

se <- SummarizedExperiment(assays=list(counts=rpkm),
                     rowRanges=gtf, colData=colData)

deseq <- DESeqDataSet(se, design = ~ replicate + temperature)
deseq$temperature = relevel(deseq$temperature, ref = "25")
res <- DESeq(deseq)
res <- results(res, contrast = c("temperature", "37", "42"))
summary(res)


lfc <- lfcShrink(deseq, coef = "replicate_2_vs_1", type = "apeglm")











deseq2 <- DESeqDataSetFromMatrix(countData = heatshock_counts,
                                 colData = colData,
                                 design = ~ replicate + temperature) %>%
  DESeq()

results(deseq2)
