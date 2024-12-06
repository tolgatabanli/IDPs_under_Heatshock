library(janitor)
library(biomaRt)
library(tidyverse)

dir <- here::here("alphafold cifs")
files <- list.files(dir)

# Read the already processed data
global_scores <- readRDS("IDP decisions/plddt_global_scores.Rds")

# Read the data new from the cif files, unnecessary if readRDS worked
global_scores <- c()
i <- 1
for (file in files) {
  print(i)
  lines <- readLines(paste0(dir, "/", file))
  plddt_line <- grep("_ma_qa_metric_global.metric_value", lines, value = TRUE)
  accession <- grep("_ma_target_ref_db_details.db_accession", lines, value = TRUE)
  accession <- strsplit(accession, "                 ")[[1]][2]

  val <- strsplit(plddt_line, " ")[[1]][2]
  global_scores[accession] <- as.numeric(val)
  i = i + 1
}
global_scores %>% saveRDS("IDP decisions/plddt_global_scores.Rds")

# biomart for accession mapping
ensembl_up_mapping <- useEnsembl(biomart = "genes", dataset = "scerevisiae_gene_ensembl") %>%
  getBM(attributes = c("ensembl_gene_id", "uniprotswissprot"),
      mart = .)

ensembl_up_mapping <- readRDS("IDP decisions/ensembl_up_mapping.Rds")

global_scores <- global_scores %>% enframe()


# For all proteins
ensembl_up_mapping %>%
  inner_join(global_scores, join_by(uniprotswissprot == name)) %>%
  dplyr::rename("lddt" = "value") %>%
  select(ensembl_gene_id, lddt) %>%
  write_tsv("IDP decisions/all_proteins_from_alphafold_global_ensembl.tsv")


# IDPs
read_tsv("IDP decisions/all_proteins_from_alphafold_global_ensembl.tsv") %>%
  filter(lddt < 50) %>%
  write_tsv("IDP decisions/idps_from_global_plddt_scores.tsv")

