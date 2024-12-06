library(tidyverse)
library(biomaRt)

# for recreating biomart
ensembl_up_mapping <- useEnsembl(biomart = "genes", dataset = "scerevisiae_gene_ensembl") %>%
  getBM(attributes = c("ensembl_gene_id", "uniprotswissprot"),
        mart = .)

# use prepared mapping from biomart
ensembl_up_mapping <- readRDS("IDP decisions/ensembl_up_mapping.Rds")



aiupred_all_proteins_mode <- read_tsv(here::here("aiupred/aiupred_modes_of_proteins.tsv"))
aiupred_all_proteins_mean <- read_tsv(here::here("aiupred/aiupred_means_of_proteins.tsv"))
aiupred_all_proteins_vote <- read_tsv(here::here("aiupred/aiupred_majority_vote.tsv"))

aiupred_all_proteins_mean %>%
  inner_join(ensembl_up_mapping, join_by(uniprot == uniprotswissprot)) %>%
  select(ensembl_gene_id, mean) %>%
  dplyr::rename("mean_disorder" = "mean") %>%
  write_tsv("IDP decisions/aiupred_all_proteins_mean.tsv")

aiupred_idps <- aiupred_all_proteins_mean %>%
  filter(mean > 0.8)
aiupred_votes <- read_tsv(here::here("aiupred/aiupred_majority_vote.tsv")) %>%
  filter(ratio > 0.5)



# save with mapped accessions
ensembl_up_mapping %>%
  right_join(aiupred_idps, join_by(uniprotswissprot == name_uniprot)) %>%
  drop_na() %>%
  write_tsv("IDP decisions/idps_from_aiupred_ensembl_uniprot.tsv")

# repeat for votes
ensembl_up_mapping %>%
  right_join(aiupred_votes, join_by(uniprotswissprot == name_uniprot)) %>%
  drop_na() %>%
  write_tsv("IDP decisions/idps_from_aiupred_votes.tsv")


