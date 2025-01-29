library(tidyverse)
library(plotly)

aiupred_modes <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_mode.tsv") %>%
  mutate(mode_disorder = mode_disorder * 100)
alphafold_modes <- read_tsv("IDP decisions/all_proteins_from_alphafold_with_mode.tsv")

res <- data.frame(disorder_threshold = numeric(),
                     plddt_threshold = numeric(),
                     number_of_commons = numeric())

for (disorder_thr in seq(20,90, 5)) {
  for (plddt_thr in seq(90, 20, -5)) {
    aiupred <- aiupred_modes %>%
      filter(mode_disorder > disorder_thr) %>%
      pull(ensembl_gene_id)
    alphafold <- alphafold_modes %>%
      filter(mode_plddt < plddt_thr) %>%
      pull(ensembl_gene_id)
    n_commons <- dplyr::intersect(aiupred, alphafold) %>%
      length()
    res <- res %>%
      add_row(disorder_threshold = disorder_thr,
                  plddt_threshold = plddt_thr,
                  number_of_commons = n_commons)
  }
}

with(res, plot_ly(x = disorder_threshold,
                  y = plddt_threshold,
                  z = number_of_commons),
     type = "scatter3d",
     mode="markers")
