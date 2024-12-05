library(janitor)
library(biomaRt)
library(tidyverse)
library(UniProt.ws)

dir <- here::here("alphafold cifs")
files <- list.files(dir)

# Read the already processed data
plddt_scores_cif <- readRDS("plddt_scores_cif.Rds")

# Read the data new from the cif files, unnecessary if readRDS worked
plddt_scores_cif <- c()
i <- 1
for (file in files) {
  print(i)
  lines <- readLines(paste0(dir, "/", file))
  plddt_line <- grep("_ma_qa_metric_global.metric_value", lines, value = TRUE)
  accession <- grep("_ma_target_ref_db_details.db_accession", lines, value = TRUE)
  accession <- strsplit(accession, "                 ")[[1]][2]

  val <- strsplit(plddt_line, " ")[[1]][2]
  plddt_scores_cif[accession] <- as.numeric(val)
  i = i + 1
}

idps <- plddt_scores_cif[which(plddt_scores_cif < 50)]


# Accession mapping
idp_df <- idps %>% enframe()

## from Evi's file, produces 5 NAs
id_mapping <- read_tsv("yeastMapping.tsv", col_names = TRUE) %>%
  janitor::clean_names() %>%
  dplyr::select(ensembl, uniprot_swissprot) %>%
  separate_longer_delim(cols = uniprot_swissprot,
                       delim = "|")

joined_mapping <- id_mapping %>%
  left_join(idp_df, join_by(uniprot_swissprot == name))
  

## using biomaRt, produces 3 NAs
ensembl_up_mapping <- useEnsembl(biomart = "genes", dataset = "scerevisiae_gene_ensembl") %>%
  getBM(attributes = c("ensembl_gene_id", "uniprotswissprot"),
      mart = .)

joined_biomart <- ensembl_up_mapping %>%
  right_join(idp_df, join_by(uniprotswissprot == name))                                

# save joined table from biomaRt
joined_biomart %>% drop_na() %>%
  rename("lddt" = "value") %>%
  write_tsv("idps_from_alphafold_ensembl_uniprot.tsv")
