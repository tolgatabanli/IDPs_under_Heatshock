library(pheatmap)
library(tidyverse)

# PREPARE
biomart_mapping <- read_tsv("IDP decisions/biomart_mapping.tsv")

idps <- readRDS("IDP decisions/commons_modes.Rds")
disprot <- read_tsv("DisProt/DisProt release_2024_12.tsv") %>%
  select(acc) %>%
  distinct() %>%
  inner_join(biomart_mapping, by = join_by(acc == uniprotswissprot)) %>%
  pull(ensembl_gene_id)
heat_matrix_control <- data.frame(rowname = disprot)

heat_matrix <- data.frame(rowname = idps)
deseq <- readRDS("deseq_reference37.Rds") %>%
  map(data.frame) %>%
  map(rownames_to_column) %>%
  map(select, rowname, log2FoldChange)

# Heatmap of IDPs
joined_matrix <- deseq %>%
  imap(~ setNames(.x, gsub("log2FoldChange", .y, names(.x)))) %>%
  map(right_join, heat_matrix) %>% # joins here with IDPs
  purrr::reduce(full_join) %>%
  drop_na() %>%
  column_to_rownames() %>%
  as.matrix()

dev.off()
png("significant plots/heatmap_ref37_idps.png",
    width = 960, height = 1080, res = 150)
pheatmap(joined_matrix, show_rownames = F,
         main = "Heatmap of IDPs")
dev.off()

# Heatmap with DisProt as control (?)
joined_matrix <- deseq %>%
  imap(~ setNames(.x, gsub("log2FoldChange", .y, names(.x)))) %>%
  map(right_join, heat_matrix_control) %>%
  purrr::reduce(full_join) %>%
  drop_na() %>%
  column_to_rownames() %>%
  as.matrix()

dev.off()
png("significant plots/heatmap_ref37_disprot.png",
    width = 960, height = 1080, res = 150)
pheatmap(joined_matrix, show_rownames = F,
         main = "Heatmap from DisProt")
dev.off()

# Heatmap with Non-IDPs
joined_matrix <- deseq %>%
  map(filter, !(rowname %in% idps)) %>%
  imap(~ setNames(.x, gsub("log2FoldChange", .y, names(.x)))) %>%
  purrr::reduce(full_join) %>%
  drop_na() %>%
  column_to_rownames() %>%
  as.matrix()

dev.off()
png("significant plots/heatmap_ref37_non-idps.png",
    width = 960, height = 1080, res = 150)
pheatmap(joined_matrix, show_rownames = F,
         main = "Heatmap of Non-IDPs")
dev.off()

