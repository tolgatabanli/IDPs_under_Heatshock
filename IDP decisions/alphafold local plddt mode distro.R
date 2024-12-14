library(tidyverse)

local_plddts <- readRDS("IDP decisions/local_plddts.RDS")

find_mode <- function(scores, bins) {
  hist_data <- hist(scores, breaks = bins, plot = FALSE)
  
  index_of_mode <- which.max(hist_data$counts)
  mode <- (hist_data$breaks[index_of_mode] + hist_data$breaks[index_of_mode + 1]) / 2
  return(mode)
}

find_mode_python <- function(scores, bins) {
  bin_edges <- seq(min(scores), max(scores), length.out = bins + 1)
  bin_counts <- hist(scores, breaks = bin_edges, plot = FALSE)$counts
  bin_centers <- (bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2
  mode_value <- bin_centers[which.max(bin_counts)]
  return(mode_value)
}

with_modes <- local_plddts %>% mutate(mode_plddt =
                          sapply(scores, find_mode_python, bins = 20)) %>%
  select(uniprot, mode_plddt)

hist(with_modes$mode_plddt, breaks = 20)

# join with ensembl and save
ensembl_up_mapping <- readRDS("IDP decisions/ensembl_up_mapping.Rds")

# For all proteins
ensembl_up_mapping %>%
  inner_join(with_modes, join_by(uniprotswissprot == uniprot)) %>%
  select(ensembl_gene_id, mode_plddt) %>%
  write_tsv("IDP decisions/all_proteins_from_alphafold_with_mode.tsv")


# IDPs, threshold
read_tsv("IDP decisions/all_proteins_from_alphafold_with_mode.tsv") %>%
  filter(mode_plddt < 60) %>%
  write_tsv("IDP decisions/idps_from_alphafold_mode.tsv")

