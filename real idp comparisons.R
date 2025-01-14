library(readxl)
library(tidyverse)

# Prepare databases
biomart_mapping <- read_tsv("IDP decisions/biomart_mapping.tsv")
alphafold_mode <- read_tsv("IDP decisions/all_proteins_from_alphafold_with_mode.tsv")
aiupred_mode <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_mode.tsv")


d2p2_idps <- readxl::read_xlsx("d2p2_experimental_idps.xlsx",
                          range = "A1:A38",
                          col_names = F)
colnames(d2p2_idps) <- "ensembl_gene_id"
d2p2_idps <- d2p2_idps %>%
  left_join(alphafold_mode) %>%
  left_join(aiupred_mode)

disprot_idps <- read_tsv("DisProt/DisProt release_2024_12.tsv") %>%
  select(acc, disorder_content) %>%
  distinct() %>%
  inner_join(biomart_mapping, by = join_by(acc == uniprotswissprot)) %>%
  select(ensembl_gene_id, disorder_content) %>%
  arrange(disorder_content) %>%
  left_join(alphafold_mode) %>%
  left_join(aiupred_mode)

# Analyse parameters of real IDPs

# d2p2 AF Mode
d2p2_idps %>%
  ggplot(aes(x = mode_plddt)) +
  geom_histogram() +
  labs(title = "Mode of pLDDT of Experimental IDPs (d2p2)")

# d2p2 AIUPred Mode
d2p2_idps %>%
  ggplot(aes(x = mode_disorder)) +
  geom_histogram() +
  labs(title = "Mode of Disorder of Experimental IDPs (d2p2)")

# DisProt AF Mode
disprot_idps %>%
  ggplot(aes(x = mode_plddt)) +
  geom_histogram() +
  labs(title = "Mode of pLDDT of DisProt IDPs")

# DisProt AIUPred Mode
disprot_idps %>%
  right_join(disprot_idps) %>%
  ggplot(aes(x = mode_disorder)) +
  geom_histogram() +
  labs(title = "Mode of Disorder of DisProt IDPs")

# Scatter pLDDT and disorder scores of real IDPs
library(reticulate)
local_plddts <- readRDS("IDP decisions/local_plddts.Rds") %>%
  rename(local_plddts = scores)
source_python("aiupred/read_pickle.py") # special script to read pickles
pickle_data <- read_pickle("aiupred/aiupred_scores.p") # interface

local_disorders <- tibble(pickle_data) %>%
  mutate(uniprot = names(pickle_data)) %>%
  rename(local_disorder = pickle_data)
  
# Todo:
read_tsv("DisProt/DisProt release_2024_12.tsv") %>%
  select(acc, disorder_content) %>%
  rename(uniprot = acc) %>%
  distinct() %>%
  inner_join(local_plddts) %>%
  inner_join(local_disorders) %>%
  select(local_plddts, local_disorder) %>%
  map(reduce, c) %>%
  data.frame() %>%
  ggplot(aes(y = local_disorder, x = local_plddts)) +
  geom_point()
