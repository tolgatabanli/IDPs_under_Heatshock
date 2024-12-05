library(tidyverse)
library(ggvenn)
library(gprofiler2)

aiupred_idps <- read_tsv(here::here("idps_from_aiupred_ensembl_uniprot.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()

aiupred_idps_voted <- read_tsv(here::here("idps_from_aiupred_votes.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()

alphafold_idps <- read_tsv(here::here("idps_from_alphafold_ensembl_uniprot.tsv")) %>%
  dplyr::select(ensembl_gene_id) %>%
  pull()

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
