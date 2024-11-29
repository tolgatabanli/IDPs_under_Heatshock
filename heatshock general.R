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

pca <- prcomp(t(rpkm)) # PCs are LCs of genes
fviz_eig(pca, addlabels = TRUE) # for scree
fviz_pca_var(pca, col.var = "black") # biplot?
plot(pca$x[, 1], pca$x[, 2]) # position of experiments of pc plot

# to find the best number of clusters of genes based on pca
pca$x[, c(1,2)] %>%
  fviz_nbclust(FUNcluster = kmeans, k.max=8)
pca$x[, c(1,2)] %>%
  eclust("kmeans", hc_metric = "euclidean", k=3)




# Now do the pc-s as experiments

pca <- prcomp(rpkm) # PCs are LCs of experiment columns
fviz_eig(pca, addlabels = TRUE) # for scree
fviz_pca_var(pca, col.var = "black") # biplot?
plot(pca$x[, 1], pca$x[, 2]) # position of experiments of pc plot



# to find the best number of clusters of genes based on pca
pca_scores <- pca$x %>% as.data.frame()

plddt_scores_cif <- readRDS("plddt_scores_cif.Rds")
idps <- which(plddt_scores_cif < 50)
names(idps) # CONVERT TO YEAST gene_id

highlight_index <- which(data$gene_id %in% names(idps))
highlight_index <- which(data$gene_id == "YPL085W")

pca_scores$highlight <- ifelse(1:nrow(pca_scores) == highlight_index, "highlight", "normal")

ggplot(pca_scores, aes(x = PC1, y = PC2, color = highlight)) +
  geom_point() +
  scale_color_manual(values = c("normal" = "black", "highlight" = "red")) +
  labs(title = "PCA with Highlighted Point") +
  theme_minimal()

pca$x[, c(1,2)] %>%
  fviz_nbclust(FUNcluster = kmeans, k.max=8)
pca$x[, c(1,2)] %>%
  log10() %>%
  drop_na() %>%
  eclust("kmeans", hc_metric = "euclidean", k=2)



