library(tidyverse)
library(ggvenn)
library(gprofiler2)

alphafold_all_proteins_local <- readRDS("IDP decisions/local_plddt_means.Rds")
aiupred_all_proteins_mean <- read_tsv("IDP decisions/aiupred_all_proteins_mean.tsv") %>%
  mutate(mean_disorder = mean_disorder*100)

inner_join(alphafold_all_proteins_local, aiupred_all_proteins_mean) %>%
  ggplot(aes(x = mean_plddt, y = mean_disorder)) +
  geom_point()


aiupred_idps <- read_tsv(here::here("IDP decisions/idps_from_aiupred_ensembl_uniprot.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()

aiupred_idps_voted <- read_tsv(here::here("IDP decisions/idps_from_aiupred_votes.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()

alphafold_idps <- read_tsv(here::here("IDP decisions/idps_from_alphafold_ensembl_uniprot.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()

alphafold_idps_withmeans <- read_tsv(here::here("IDP decisions/idps_from_alphafold_ensembl_uniprot.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()


# for corr plots
alphafold_prots <- read_tsv(here::here("IDP decisions/all_proteins_from_alphafold_ensembl_uniprot.tsv")) %>%
  dplyr::select(uniprotswissprot, lddt) %>%
  dplyr::rename("uniprot" = "uniprotswissprot")
aiupred_all_mode <- read_tsv(here::here("aiupred/aiupred_modes_of_proteins.tsv"))
aiupred_all_vote <- read_tsv(here::here("aiupred/aiupred_majority_vote.tsv")) %>%
  dplyr::rename("uniprot" = "name_uniprot")
aiupred_all_mean <- read_tsv(here::here("aiupred/aiupred_means_of_proteins.tsv"))

af_mode <- alphafold_prots %>% inner_join(aiupred_all_mode)
af_vote <- alphafold_prots %>% inner_join(aiupred_all_vote)
af_mean <- alphafold_prots %>% inner_join(aiupred_all_mean)
plot(af_mode$lddt, af_mode$mode)
plot(af_vote$lddt, af_vote$ratio)
plot(af_mean$lddt, af_mean$mean)


# AIUPred modes vs AlphaFold
ggvenn(list(aiupred = aiupred_idps, alphafold = alphafold_idps))
commons <- dplyr::intersect(aiupred_idps, alphafold_idps)
gost(commons,
     organism = "scerevisiae",
     correction_method = "fdr")$result$term_name

# AIUPred Votes vs AlphaFold
ggvenn(list(aiupred = aiupred_idps_voted, alphafold = alphafold_idps))
commons <- dplyr::intersect(aiupred_idps_voted, alphafold_idps)
gost(commons,
     organism = "scerevisiae",
     correction_method = "fdr")$result$term_name
