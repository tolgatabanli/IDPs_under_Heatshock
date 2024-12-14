library(tidyverse)
library(ggvenn)
library(gprofiler2)
library(ggpubr)

alphafold_all_proteins_local <- read_tsv("IDP decisions/all_proteins_from_alphafold_local_ensembl.tsv")
alphafold_all_proteins_global <- read_tsv("IDP decisions/all_proteins_from_alphafold_global_ensembl.tsv")
alphafold_all_proteins_mode <- read_tsv("IDP decisions/all_proteins_from_alphafold_with_mode.tsv")
aiupred_all_proteins_mean <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_mean.tsv") %>%
  mutate(mean_disorder = mean_disorder*100)
aiupred_all_proteins_mode <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_mode.tsv") %>%
  mutate(mode_disorder = mode_disorder*100)
aiupred_all_proteins_perc_alpha50 <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_perc_alpha50.tsv") %>%
  mutate(perc = perc*100) %>%
  rename(perc50 = perc)
aiupred_all_proteins_perc_alpha80 <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_perc_alpha80.tsv") %>%
  mutate(perc = perc*100) %>%
  rename(perc80 = perc)

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
  ggplot(aes(x = mean_disorder, y = perc50)) +
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
  ggplot(aes(x = mode_disorder, y = perc50)) +
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
  ggplot(aes(x = mode_disorder, y = perc80)) +
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

# AIUPred Move and Alphafold Mode
ggvenn(list(aiupred_mode = aiupred_idps_mode, alphafold_mode = alphafold_idps_mode))
commons <- dplyr::intersect(aiupred_idps_mode, alphafold_idps_mode)
commons %>% saveRDS("IDP decisions/commons_modes.Rds")
gost(commons,
     organism = "scerevisiae",
     correction_method = "bonferroni")$result[c("term_name", "p_value")]
