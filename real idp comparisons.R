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

# Scatter plddt and disorder modes
disprot_idps %>%
  ggplot(aes(x = mode_plddt, y = mode_disorder)) +
  geom_point()

# Scatter pLDDT and disorder scores of real IDPs
library(reticulate)
local_plddts <- readRDS("IDP decisions/local_plddts.Rds") %>%
  rename(local_plddts = scores)
source_python("aiupred/read_pickle.py") # special script to read pickles
pickle_data <- read_pickle("aiupred/aiupred_scores.p") # interface

local_disorders <- tibble(pickle_data) %>%
  mutate(uniprot = names(pickle_data)) %>%
  rename(local_disorder = pickle_data)
  
# PREPARE Residue-wise
disprot_joined <- read_tsv("DisProt/DisProt release_2024_12.tsv") %>%
  select(uniprot = acc, disorder_content) %>%
  distinct() %>%
  inner_join(local_plddts) %>%
  inner_join(local_disorders) %>%
  select(uniprot, local_plddts, local_disorder)

# chaotic full scatter!!!
disprot_joined %>%
  unnest(cols = c(local_plddts, local_disorder)) %>%
  ggplot(aes(y = local_disorder, x = local_plddts)) +
  geom_point()

# violin plot
disprot_joined %>%
  unnest(cols = c(local_plddts, local_disorder)) %>%
  mutate(local_disorder = local_disorder * 100) %>%
  pivot_longer(names_to = "score_type",
               cols = c("local_plddts", "local_disorder")) %>%
  ggplot(aes(x = score_type, y = value)) +
  geom_violin() + 
  labs(title = "Violin Chart for Two Types of Scores")

# with sampling for scatter
disprot_joined %>%
  unnest(cols = c(local_plddts, local_disorder)) %>%
  group_by(uniprot) %>%
  sample_n(size = min(100, n()), replace = FALSE) %>%
  ungroup() %>%
  ggplot(aes(x = local_plddts, y = local_disorder)) +
  geom_point(alpha = 0.5) +
  labs(title = "Scatter Plot of all residue-wise scores")

disprot_joined %>%
  unnest(cols = c(local_plddts, local_disorder)) %>%
  ggplot(aes(x = local_plddts, y = local_disorder)) +
  geom_hex(bins = 10)

disprot_joined %>%
  unnest(cols = c(local_plddts, local_disorder)) %>%
  ggplot(aes(x = local_plddts, y = local_disorder)) +
  #stat_density_2d() +
  geom_density_2d()

# summarize per prot
disprot_joined %>%
  unnest(cols = c(local_plddts, local_disorder)) %>%
  group_by(uniprot) %>%
  summarize(plddt_mean = mean(local_plddts, na.rm = TRUE),
            disorder_mean = mean(local_disorder, na.rm = TRUE)) %>%
  ggplot(aes(x = plddt_mean, y = disorder_mean)) +
  geom_point(alpha = 0.7) +
  labs(title = "Mean Scores per Protein")
