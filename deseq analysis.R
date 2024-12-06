library(DESeq2)
library(tidyverse)
library(SummarizedExperiment)
library(edgeR)

# You should transform to FPKM or RPKM

colData <- read_tsv("sample_mapping.tsv", col_names = T) %>%
  mutate(across(everything(), as.factor)) %>%
  column_to_rownames(var = "sample")

heatshock_counts <- read_tsv(here::here("complex_yeast_heatshock.tsv"), col_names = T) %>%
  column_to_rownames(var = "gene_id") %>%
  relocate(rownames(colData)) %>%
  as.matrix()

se <- SummarizedExperiment(assays=list(counts=heatshock_counts),
                           colData=colData)

temp25_t0 <- colData %>%
  filter(temperature == 25, time == 0)

knockouts <- colData %>%
  select(knockout) %>%
  distinct() %>%
  pull()

ref = "Wildtype"

for (knock in knockouts) {
  subsample <- temp25_t0 %>% filter(knockout %in% c(ref, as.character(knock)))
  
  subcounts <- heatshock_counts %>% select(rownames(subsample)) # TODO
  
  dds <- DESeqDataSet(countData = subcounts,
                      colData = subsample,
                      design = ~ knockout)
  dds$temperature = relevel(deseq$temperature, ref = "37")
  deseq <- DESeq(deseq)
}

# DESeq


# Results, Contrastin and Plotting
res <- results(deseq)
res <- results(res, contrast = c("temperature", "37", "42"))
summary(res)

DESeq2::plotMA(res)

# LFC ?
lfc <- lfcShrink(res, coef = "temperature_42_vs_37", type = "apeglm")
plotMA(lfc)
