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

# Change of the ratio of DisProt IDPs with disorder_mode > 60
# with disorder_content threshold
disorder_content_thresholds <- seq(0,100,5)
total_num <- nrow(disprot_idps)
ratio_alpha60 <- c()
i <- 1
for (threshold in disorder_content_thresholds) {
  df <- disprot_idps %>%
    filter(disorder_content > threshold)
  # disorder score ratios:
  ratio_alpha60[i] <- df %>%
    filter(mode_plddt > 60) %>%
    nrow() %>%
    `/`(nrow(df))
  i <- i + 1
}
plot(disorder_content_range, ratio_alpha60)
