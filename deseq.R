library(DESeq2)
library(tidyverse)
library(edgeR)
library(ggpubr)
library(plotly)
library(grid)

colData <- read_tsv("sample_mapping.tsv", col_names = T) %>%
  mutate(across(everything(), as.factor)) %>%
  column_to_rownames(var = "sample")

heatshock_counts <- read_tsv(here::here("complex_yeast_heatshock.tsv"), col_names = T) %>%
  column_to_rownames(var = "gene_id") %>%
  relocate(rownames(colData))

knockouts <- colData %>%
  select(knockout) %>%
  distinct() %>%
  pull() %>%
  relevel("Wildtype")

idps <- readRDS("IDP decisions/commons_modes.Rds")
deseq_results <- list()
deseq_results <- readRDS("deseq_results.Rds") # Read already done deseq results

# returns two plots comparing with temp=25 and t=0 of the respective strain
multiple_deseq <- function(countData, colData, knock, temp, t) {
  subColData <- colData %>%
    filter(knockout == knock) %>%
    filter((temperature == 25 & time == 0) | temperature == temp)
  
  subsample <- subColData %>%
    filter(time %in% c(0, t))
  subcounts <- heatshock_counts %>% select(rownames(subsample))
  
  dds <- DESeqDataSetFromMatrix(countData = subcounts,
                                colData = subsample,
                                design = ~ time)
  deseq <- DESeq(dds, quiet = TRUE)
  res <- results(deseq)
  
  deseq_results[[paste0(knock,"_",temp,"_",t)]] <<- res
  
  res <- res %>%
    as.data.frame() %>%
    drop_na()
  significant_idps <- res$padj < 0.05 & rownames(res) %in% idps
  grob <- grobTree(textGrob(
    paste("Significant IDP ratio:\n",
          round(sum(significant_idps) / length(idps), digits = 2)),
    x=0.5,  y=0.6))
  
  p <- ggplot(data = res,
              aes(x = log2FoldChange, y = -log10(pvalue),
                  colour = significant_idps)) +
    scale_color_manual(values = c("#0400ff", "#747880"),
                       labels = c("Significant IDPs", "Others"),
                       breaks = c(TRUE, FALSE)) +
    geom_point() +
    annotation_custom(grob) +
    labs(title = paste(knock, temp, "C, t =", t),
         color = "Significance")
  return(p) # TODO: ggplot has two fields that it returns
}


# DESeq
times <- c(10, 30)
temps <- c(37, 42)
knockouts <- colData %>%
  select(knockout) %>%
  unique() %>%
  pull()

# First changing time, then temp and lastly knockouts.
combs <- expand.grid(time = times, temperature = temps, knockout = knockouts)
plots <- list()
for (row in 1:nrow(combs)) {
  plots[[row]] <- multiple_deseq(heatshock_counts,
                                 colData,
                                 combs[row,"knockout"],
                                 combs[row,"temperature"],
                                 combs[row, "time"])
}

# TODO: the color aesthetic throws length errors!
ggarrange(plotlist = plots, ncol = 2, nrow = 8, common.legend = T)
ggsave("volcanoes.png", device = "png", height = 25, width = 10)
