library(tidyverse)
library(ggvenn)
library(gprofiler2)
library(ggpubr)

alphafold_all_proteins_local <- read_tsv("IDP decisions/all_proteins_from_alphafold_local_ensembl.tsv")
alphafold_all_proteins_global <- read_tsv("IDP decisions/all_proteins_from_alphafold_global_ensembl.tsv")
aiupred_all_proteins_mean <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_mean.tsv") %>%
  mutate(mean_disorder = mean_disorder*100)
aiupred_all_proteins_mode <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_mode.tsv") %>%
  mutate(mode_disorder = mode_disorder*100)

# plot all proteins from two DBs against each other
inner_join(alphafold_all_proteins_local, aiupred_all_proteins_mean) %>%
  ggplot(aes(x = mean_plddt, y = mean_disorder)) +
  geom_point() +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.value), sep = "~`,`~")),
           method="pearson", p.accuracy = 0.001) +
  labs(title = "Mean AIUPred scores against mean of local pLDDTs") +
  xlab("Mean pLDDT") +
  ylab("Mean disorder score")
inner_join(alphafold_all_proteins_global, aiupred_all_proteins_mean) %>%
  ggplot(aes(x = lddt, y = mean_disorder)) +
  geom_point() +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.value), sep = "~`,`~")),
           method="pearson") +
  labs(title = "Mean AIUPred scores against global pLDDTs") +
  xlab("Global pLDDT") +
  ylab("Mean disorder score")


aiupred_idps_mean <- read_tsv(here::here("IDP decisions/idps_from_aiupred_mean_ensembl.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()

aiupred_idps_vote <- read_tsv(here::here("IDP decisions/idps_from_aiupred_vote_ensembl.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()

alphafold_idps_global <- read_tsv(here::here("IDP decisions/idps_from_global_plddt_scores.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()

alphafold_idps_local <- read_tsv(here::here("IDP decisions/idps_from_local_plddt_scores.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()


# AIUPred means vs AlphaFold global
ggvenn(list(aiupred = aiupred_idps_mean, alphafold = alphafold_idps_global))
commons <- dplyr::intersect(aiupred_idps_mean, alphafold_idps_global)
gost(commons,
     organism = "scerevisiae",
     correction_method = "fdr")$result$term_name

# AIUPred Votes vs AlphaFold
ggvenn(list(aiupred = aiupred_idps_vote, alphafold = alphafold_idps_global))
commons <- dplyr::intersect(aiupred_idps_vote, alphafold_idps_global)
gost(commons,
     organism = "scerevisiae",
     correction_method = "fdr")$result[c("term_name", "p_value")]
