library(tidyverse)
library(biomaRt)

aiupred_idps <- read_tsv(here::here("aiupred/aiupred_modes_of_proteins.tsv")) %>%
  filter(score_mode > 0.8)

aiupred_votes <- read_tsv(here::here("aiupred/aiupred_majority_vote.tsv")) %>%
  filter(ratio > 0.5)

# for recreating biomart
ensembl_up_mapping <- useEnsembl(biomart = "genes", dataset = "scerevisiae_gene_ensembl") %>%
  getBM(attributes = c("ensembl_gene_id", "uniprotswissprot"),
        mart = .)

# use prepared mapping from biomart
ensembl_up_mapping <- readRDS("ensembl_up_mapping.Rds")

# save wtih mapped accessions
ensembl_up_mapping %>%
  right_join(aiupred_idps, join_by(uniprotswissprot == name_uniprot)) %>%
  drop_na() %>%
  write_tsv("idps_from_aiupred_ensembl_uniprot.tsv")

# repeat for votes
ensembl_up_mapping %>%
  right_join(aiupred_votes, join_by(uniprotswissprot == name_uniprot)) %>%
  drop_na() %>%
  write_tsv("idps_from_aiupred_votes.tsv")


