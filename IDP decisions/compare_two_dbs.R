library(tidyverse)
library(ggvenn)
library(gprofiler2)
library(ggpubr)

alphafold_all_proteins_local <- read_tsv("IDP decisions/all_proteins_from_alphafold_with_mean.tsv")
alphafold_all_proteins_global <- read_tsv("IDP decisions/all_proteins_from_alphafold_global_ensembl.tsv")
alphafold_all_proteins_mode <- read_tsv("IDP decisions/all_proteins_from_alphafold_with_mode.tsv")
alphafold_all_proteins_perc_alpha60 <- read_tsv("IDP decisions/all_proteins_from_alphafold_with_perc_alpha60.tsv")
aiupred_all_proteins_mean <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_mean.tsv") %>%
  mutate(mean_disorder = mean_disorder*100)
aiupred_all_proteins_mode <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_mode.tsv") %>%
  mutate(mode_disorder = mode_disorder*100)
aiupred_all_proteins_perc_alpha50 <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_perc_alpha50.tsv") %>%
  mutate(aiupred_perc_50 = aiupred_perc_50*100)
aiupred_all_proteins_perc_alpha60 <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_perc_alpha60.tsv") %>%
  mutate(aiupred_perc_60 = aiupred_perc_60*100)
aiupred_all_proteins_perc_alpha80 <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_perc_alpha80.tsv") %>%
  mutate(aiupred_perc_80 = aiupred_perc_80*100)


# Comparison plots

## Within database

### AF local vs global
inner_join(alphafold_all_proteins_local, alphafold_all_proteins_global) %>%
  ggplot(aes(x = mean_plddt, y = lddt)) +
  geom_point() +
  stat_cor(aes(label = paste(after_stat(rr.label),
                             after_stat(p.value), sep = "~`, p =`~")),
           method="pearson", p.accuracy = 0.001) +
  labs(title = "Global pLDDT scores against mean of local pLDDTs") +
  xlab("Mean pLDDT scores") +
  ylab("Global pLDDT scores")


### AIUPred mean vs perc50
inner_join(aiupred_all_proteins_mean, aiupred_all_proteins_perc_alpha50) %>%
  ggplot(aes(x = mean_disorder, y = aiupred_perc_50)) +
  geom_point() +
  stat_cor(aes(label = paste(after_stat(rr.label),
                             after_stat(p.value), sep = "~`, p =`~")),
           method="pearson", p.accuracy = 0.001) +
  labs(title = expression(
    paste("AIUPred Mean vs ", Percent[alpha==50], " Disorder Scores"))) +
  xlab("Mean disorder") +
  ylab("Percent disorder")

### AIUPred perc50 vs mode
inner_join(aiupred_all_proteins_mode, aiupred_all_proteins_perc_alpha50) %>%
  ggplot(aes(x = mode_disorder, y = aiupred_perc_50)) +
  geom_point() +
  stat_cor(aes(label = paste(after_stat(rr.label),
                             after_stat(p.value), sep = "~`, p =`~")),
           method="pearson", p.accuracy = 0.001) +
  labs(title = expression(
    paste("AIUPred Mode vs ", Percent[alpha==50], " Disorder Scores"))) +
  xlab("Mode disorder") +
  ylab("Percent disorder")

### AIUPred perc80 vs mode
inner_join(aiupred_all_proteins_mode, aiupred_all_proteins_perc_alpha80) %>%
  ggplot(aes(x = mode_disorder, y = aiupred_perc_80)) +
  geom_point() +
  stat_cor(aes(label = paste(after_stat(rr.label),
                             after_stat(p.value), sep = "~`, p =`~")),
           method="pearson", p.accuracy = 0.001) +
  labs(title = expression(
    paste("AIUPred Mean vs ", Percent[alpha==80], " Disorder Scores"))) +
  xlab("Mode disorder") +
  ylab("Percent disorder")

## Between DBs
### AF global vs AIUPred Mean
inner_join(alphafold_all_proteins_global, aiupred_all_proteins_mean) %>%
  ggplot(aes(x = lddt, y = mean_disorder)) +
  geom_point() +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.value), sep = "~`,`~")),
           method="pearson") +
  labs(title = "Mean AIUPred scores against global pLDDTs") +
  xlab("Global pLDDT") +
  ylab("Mean disorder score")

### AF global vs AIUPred Mode - SCHLECHT MAYBE DELETE
inner_join(alphafold_all_proteins_global, aiupred_all_proteins_mode) %>%
  ggplot(aes(x = lddt, y = mode_disorder)) +
  geom_point() +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.value), sep = "~`,`~")),
           method="pearson") +
  labs(title = "Mode AIUPred scores against global pLDDTs") +
  xlab("Global pLDDT") +
  ylab("Mode disorder score")

### AF Mode vs AIUPred Mode
inner_join(alphafold_all_proteins_mode, aiupred_all_proteins_mode) %>%
  ggplot(aes(x = mode_disorder, y = mode_plddt)) +
  geom_point() +
  labs(title = "AIUPred Mode vs AlphaFold Mode") +
  xlab("Mode disorder") +
  ylab("Mode pLDDT")

### AF perc60 vs AIUPred perc60
inner_join(alphafold_all_proteins_perc_alpha60, aiupred_all_proteins_perc_alpha60) %>%
  ggplot(aes(x = alphafold_perc_60, y = aiupred_perc_60)) +
  geom_point() +
  labs(title = "AIUPred Perc60 vs AlphaFold Perc60") +
  xlab("Percent Disorder AlphaFold") +
  ylab("Percent Disorder AIUPred")


aiupred_idps_mean <- read_tsv(here::here("IDP decisions/idps_from_aiupred_mean_ensembl.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()
aiupred_idps_perc <- read_tsv(here::here("IDP decisions/idps_from_aiupred_perc_alpha80_ensembl.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()
aiupred_idps_mode <- read_tsv(here::here("IDP decisions/idps_from_aiupred_mode_ensembl.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()
alphafold_idps_global <- read_tsv(here::here("IDP decisions/idps_from_global_plddt_scores.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()
alphafold_idps_local <- read_tsv(here::here("IDP decisions/idps_from_local_plddt_scores.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()
alphafold_idps_mode <- read_tsv(here::here("IDP decisions/idps_from_alphafold_mode.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()


# AIUPred means vs AlphaFold global
ggvenn(list(aiupred = aiupred_idps_mean, alphafold = alphafold_idps_global))
commons <- dplyr::intersect(aiupred_idps_mean, alphafold_idps_global)
gost(commons,
     organism = "scerevisiae",
     correction_method = "bonferroni")$result$term_name

# AIUPred Percentages vs AlphaFold
ggvenn(list(aiupred = aiupred_idps_perc, alphafold = alphafold_idps_global))
commons <- dplyr::intersect(aiupred_idps_perc, alphafold_idps_global)
gost(commons,
     organism = "scerevisiae",
     correction_method = "bonferroni")$result[c("term_name", "p_value")]

# AIUPred Percentages vs Mode
ggvenn(list(aiupred_perc = aiupred_idps_perc, aiupred_mode = aiupred_idps_mode))
commons <- dplyr::intersect(aiupred_idps_perc, aiupred_idps_mode)
gost(commons,
     organism = "scerevisiae",
     correction_method = "bonferroni")$result[c("term_name", "p_value")]

# AIUPred Mode and Alphafold Mode
ggvenn(list(aiupred_mode = aiupred_idps_mode, alphafold_mode = alphafold_idps_mode))
commons <- dplyr::intersect(aiupred_idps_mode, alphafold_idps_mode)
commons %>% saveRDS("IDP decisions/commons_modes.Rds")
gost(commons,
     organism = "scerevisiae",
     correction_method = "bonferroni") %>%
  gostplot()






# PCA
library(factoextra)
library(corrr)

central_tendencies <- list(alphafold_all_proteins_global,
                           alphafold_all_proteins_mode,
                           alphafold_all_proteins_perc_alpha60,
                           aiupred_all_proteins_mean,
                           aiupred_all_proteins_mode,
                           aiupred_all_proteins_perc_alpha60)

central_matrix <- central_tendencies %>%
  reduce(inner_join, multiple = "any")
protein_names <- central_matrix$ensembl_gene_id

pca <- central_matrix %>%
  select(-ensembl_gene_id) %>%
  scale() %>%
  prcomp()

fviz_eig(pca, addlabels = TRUE)
fviz_pca_var(pca, repel = T)
fviz_pca_ind(pca)
fviz_pca_biplot(pca)

pca_2d <- pca$x[,1:2] %>% as.data.frame()
fviz_nbclust(pca_2d, kmeans)
pca_cluster <- kmeans(pca_2d, centers = 3)
fviz_cluster(pca_cluster, data = pca_2d)

disorder_clusters <- data.frame(ensembl_gene_id = protein_names,
                                cluster = pca_cluster$cluster)
disorder_clusters %>%
  filter(cluster == 1) %>%
  pull(ensembl_gene_id) %>%
  saveRDS("idps_from_pca_from_central_tendencies.Rds")
