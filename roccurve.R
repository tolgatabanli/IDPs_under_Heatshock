library(tidyverse)

# order the proteins by the disorder score 
# and calculate for each disorder treshold how many
# disordered proteins are also significant in the DEseq data

prot_disorder <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_mode.tsv") %>%
  mutate(mode_disorder = mode_disorder * 100)

prot_plddt <- read_tsv("IDP decisions/all_proteins_from_alphafold_with_mode.tsv") %>%
  # arrange(mode_plddt) %>%
  filter(mode_plddt < 60)

prots <- inner_join(prot_disorder, prot_plddt)

ggplot(prot_disorder, aes(x = mode_disorder)) +
  geom_histogram()

deseq_sigs <- readRDS("deseq_results.Rds") %>%
  map(data.frame) %>%
  map(rownames_to_column, var = "ensembl_gene_id") %>%
  map(select, ensembl_gene_id, padj) %>%
  map(filter, padj < 0.05) %>%
  map(inner_join, prot_disorder) %>%
  map(arrange, desc(mode_disorder))


# start with the empty set of IDPs
# (cutoff=highest disorder score -> no proteins more disordered than that)
# then iteratively expand the set of IDPs
# by lowering the cutoff to the next disorder score.

thresholds <- prot_disorder %>%
  pull(mode_disorder) %>%
  unique()
thr_list <- list()
for (thr in thresholds) {
  idps <- prot_disorder %>%
    filter(mode_disorder > thr) %>%
    pull(ensembl_gene_id)
  thr_list[[toString(thr)]] <- list(idps)
}

# iwie zu langsam
filter_while <- function(sorted_df, condition) {
  cond <- substitute(condition)
  
  stopifnot(is.data.frame(sorted_df))
  result <- data.frame()
  
  for (i in seq_len(nrow(sorted_df))) {
    row <- sorted_df[i, , drop = F]
    if (!eval(cond, envir = row)) {
      break
    }
    result <- rbind(result, row)
  }
  
  return(result)
}

ratios <- data.frame(thr = numeric(),
                     ratio = numeric(),
                     cond = character())

for (cond in names(deseq_sigs)) {
  for (thr in thresholds) {
    idps <- thr_list[toString(thr)] %>% unlist() %>% unname()
    n_total <- length(idps)

    n_sig <- deseq_sigs[[cond]] %>%
      filter(mode_disorder > thr) %>%
      nrow()

    ratios <- ratios %>%
      add_row(thr = thr,
              ratio = n_sig / n_total,
              cond = cond)
  }
  print(cond)
}

ratios %>%
  ggplot(aes(x = thr, y = ratio)) +
  geom_line() +
  facet_wrap(~ cond)



# OR go through significant lfc s 

