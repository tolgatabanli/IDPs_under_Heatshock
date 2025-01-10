library(tidyverse)
library(FactoMineR)
library(factoextra)

deseq_results <- readRDS("deseq_results.Rds")
gene_names <- rownames(deseq_results[[1]])

dir <- "IDP decisions/"
files <- list.files(
  path = dir,
  pattern = "^all_proteins",
  full.names = TRUE
)
files <- files[!(grepl("alpha80.tsv$", files) | (grepl("ensembl.tsv$", files)))]

idp_scores <- files %>%
  map(read_tsv) %>%
  purrr::reduce(inner_join, by = "ensembl_gene_id")


deseq <- deseq_results %>%
  map(as.data.frame) %>%
  map(dplyr::select, log2FoldChange) %>%
  list_cbind() %>%
  as.data.frame() %>%
  add_column(gene_names, .before = 1)

d <- deseq %>%
  inner_join(idp_scores, by = join_by(gene_names == ensembl_gene_id)) %>%
  drop_na()

pca <- scale(d[-1]) %>%
  princomp()

fviz_eig(pca)
ggsave("scree_plot.jpeg")

fviz_pca_var(pca, repel = TRUE)
ggsave("biplot.jpeg")

fviz_cos2(pca, choice = "var", axes = 1:2)
fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)


