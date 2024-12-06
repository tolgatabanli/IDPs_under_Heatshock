library(tidyverse)

dir <- here::here("alphafold cifs")
files <- list.files(dir)

local_plddts <- tibble(uniprot = character(), scores = list())
start <- "_ma_qa_metric_local.ordinal_id"
end <- "#"
score_pattern <- "\\b\\d+\\.\\d+\\b"  # regex for floating point numbers

i <- 1
for (file in files) {
  print(i) # for progress

  accession <- str_split_1(file, "-")[2] # get accession from the file name
  file_connection <- file(paste0(dir, "/", file), open = "r")
  found_start <- FALSE
  scores <- c()
  while(TRUE) {
    line <- readLines(file_connection, n = 1)
    if(!found_start) {
      if (line == start) found_start <- TRUE
    } else {
      if (line == end) {
        break
      }
      # extract the decimal number and concatenate to the scores
      scores <- c(scores, as.numeric(regmatches(line, regexpr(score_pattern, line))))
    }
  }
  local_plddts <- local_plddts %>% add_row(uniprot = accession, scores = list(scores))
  i = i + 1
}
saveRDS(local_plddts, "IDP decisions/local_plddts.Rds")



local_plddts <- readRDS("IDP decisions/local_plddts.Rds")

ensembl_up_mapping <- readRDS("IDP decisions/ensembl_up_mapping.Rds")

local_plddts %>%
  group_by(uniprot) %>%
  mutate(mean_plddt = mean(unlist(scores))) %>%
  ungroup() %>%
  inner_join(ensembl_up_mapping, join_by(uniprot == uniprotswissprot)) %>%
  select(ensembl_gene_id, mean_plddt) %>%
  write_tsv("IDP decisions/all_proteins_from_alphafold_local_ensembl.tsv")

# IDPs
means <- read_tsv("IDP decisions/all_proteins_from_alphafold_local_ensembl.tsv")
means %>%
  filter(mean_plddt < 50) %>%
  write_tsv("IDP decisions/idps_from_local_plddt_scores.tsv")
