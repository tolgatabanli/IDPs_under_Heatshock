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

aiupred_all_proteins_mode %>%
  inner_join(ensembl_up_mapping, join_by(uniprot == uniprotswissprot)) %>%
  select(ensembl_gene_id, mode) %>%
  dplyr::rename("mode_disorder" = "mode") %>%
  write_tsv("IDP decisions/all_proteins_from_aiupred_with_mode.tsv")
aiupred_all_proteins_mean %>%
  inner_join(ensembl_up_mapping, join_by(uniprot == uniprotswissprot)) %>%
  select(ensembl_gene_id, mean) %>%
  dplyr::rename("mean_disorder" = "mean") %>%
  write_tsv("IDP decisions/all_proteins_from_aiupred_with_mean.tsv")
aiupred_all_proteins_vote %>%
  inner_join(ensembl_up_mapping, join_by(uniprot == uniprotswissprot)) %>%
  select(ensembl_gene_id, ratio) %>%
  dplyr::rename("disorder_ratio" = "ratio") %>%
  write_tsv("IDP decisions/all_proteins_from_aiupred_with_vote.tsv")


# IDPs
read_tsv("IDP decisions/all_proteins_from_aiupred_with_mean.tsv") %>%
  filter(mean_disorder > 0.8) %>%
  write_tsv("IDP decisions/idps_from_aiupred_mean_ensembl.tsv")

read_tsv("IDP decisions/all_proteins_from_aiupred_with_vote.tsv") %>%
  filter(disorder_ratio > 0.5) %>%
  write_tsv("IDP decisions/idps_from_aiupred_vote_ensembl.tsv")

# TODO: Mode also


