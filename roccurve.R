library(tidyverse)

# order the proteins by the disorder score 
# and calculate for each disorder treshold how many
# disordered proteins are also significant in the DEseq data

prot_disorder <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_mode.tsv") %>%
  mutate(mode_disorder = mode_disorder * 100) %>%
  arrange(desc(mode_disorder))

ggplot(prot_disorder, aes(x = mode_disorder)) +
  geom_histogram()

deseq_sigs <- readRDS("deseq_results.Rds") %>%
  map(data.frame) %>%
  map(rownames_to_column) %>%
  map(select, rowname, padj) %>%
  map(filter, padj < 0.05) %>%
  map(pull, rowname) %>%
  purrr::reduce(base::union)



# start with the empty set of IDPs
# (cutoff=highest disorder score -> no proteins more disordered than that)
# then iteratively expand the set of IDPs
# by lowering the cutoff to the next disorder score.

thresholds <- prot_disorder %>%
  pull(mode_disorder)
ratios <- c()

for (thr in thresholds) {
  idps <- prot_disorder %>%
    filter(mode_disorder > thr) %>%
    pull(ensembl_gene_id)
  n_total <- length(idps)
  
  n_sig <- deseq_sigs %>%
    intersect(idps) %>%
    length()
  
  ratios <- c(ratios, n_sig / n_total)
}

plot(thresholds, ratios)
