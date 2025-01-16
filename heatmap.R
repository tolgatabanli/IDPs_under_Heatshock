library(pheatmap)
library(tidyverse)

# PREPARE
biomart_mapping <- read_tsv("IDP decisions/biomart_mapping.tsv")

commons <- readRDS("IDP decisions/commons_modes.Rds")
disprot <- read_tsv("DisProt/DisProt release_2024_12.tsv") %>%
  select(acc) %>%
  distinct() %>%
  inner_join(biomart_mapping, by = join_by(acc == uniprotswissprot)) %>%
  pull(ensembl_gene_id)
heat_matrix_control <- data.frame(rowname = disprot)

heat_matrix <- data.frame(rowname = commons)

# Make heatmap matrix for analyticals
df_names <- names(deseq)
deseq <- readRDS("deseq_results.Rds") %>%
  map(data.frame) %>%
  map(rownames_to_column) %>%
  map(select, rowname, log2FoldChange) %>%
  imap(~ setNames(.x, gsub("log2FoldChange", .y, names(.x)))) %>%
  map(right_join, heat_matrix) %>%
  purrr::reduce(function(x, y) {
    full_join(x, y)
  }) %>%
  drop_na() %>% # 498 elements become 491 elements
  column_to_rownames() %>%
  as.matrix()

pheatmap(deseq, show_rownames = F,
         main = "Heatmap from Commons")

# Heatmap with DisProt as control (?)
df_names <- names(deseq)
deseq <- readRDS("deseq_results.Rds") %>%
  map(data.frame) %>%
  map(rownames_to_column) %>%
  map(select, rowname, log2FoldChange) %>%
  imap(~ setNames(.x, gsub("log2FoldChange", .y, names(.x)))) %>%
  map(right_join, heat_matrix_control) %>%
  purrr::reduce(function(x, y) {
    full_join(x, y)
  }) %>%
  drop_na() %>%
  column_to_rownames() %>%
  as.matrix()

pheatmap(deseq, show_rownames = F,
         main = "Heatmap from DisProt")

# Heatmap with All Proteins
df_names <- names(deseq)
deseq <- readRDS("deseq_results.Rds") %>%
  map(data.frame) %>%
  map(rownames_to_column) %>%
  map(select, rowname, log2FoldChange) %>%
  imap(~ setNames(.x, gsub("log2FoldChange", .y, names(.x)))) %>%
  purrr::reduce(function(x, y) {
    full_join(x, y)
  }) %>%
  drop_na() %>%
  column_to_rownames() %>%
  as.matrix()

pheatmap(deseq, show_rownames = F,
         main = "Heatmap with All Genes")
