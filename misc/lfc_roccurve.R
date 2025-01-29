library(tidyverse)

prot_disorder <- read_tsv("IDP decisions/all_proteins_from_aiupred_with_mode.tsv") %>%
  mutate(mode_disorder = mode_disorder * 100) %>%
  arrange(desc(mode_disorder))

deseq_all <- readRDS("deseq_results.Rds") %>%
  map(data.frame) %>%
  map(rownames_to_column, var = "ensembl_gene_id") %>%
  map(select, ensembl_gene_id, log2FoldChange, padj) %>%
  map(drop_na) %>%
  map(inner_join, prot_disorder) %>%
  map(arrange, desc(abs(log2FoldChange)))

idps <- readRDS("IDP decisions/commons_modes.Rds")

# ROC
ratios <- data.frame(thr = numeric(),
                     ratio = numeric(),
                     cond = character())

for (cond in names(deseq_all)) {
  thresholds <- deseq_all[[cond]] %>%
    pull(log2FoldChange) %>%
    abs() %>%
    unique()
  for (thr in thresholds) {
    lfc_subset <- deseq_all[[cond]] %>%
      filter(abs(log2FoldChange) > thr)
    
    # Find IDPs in the sub-experiment
    n_total <- deseq_all[[cond]] %>%
      pull(ensembl_gene_id) %>%
      intersect(idps) %>%
      length()
    
    # Find significant IDPs in the sub-experiment
    n_sig <- lfc_subset %>%
      pull(ensembl_gene_id) %>%
      intersect(idps) %>%
      length()
    
    ratios <- ratios %>%
      add_row(thr = thr,
              cond = cond,
              ratio = n_sig / n_total,
      )
  }
  print(cond)
}

ratios %>%
  ggplot(aes(x = thr, y = ratio)) +
  geom_line() +
  facet_wrap(~ cond) +
  xlab("abs(LFC) Cut-off") +
  ylab("Ratio (Significant IDP / All IDPs in the Condition)")

# Distro of LFC
deseq_all %>%
  map(drop_na) %>%
  map(filter, padj < 0.05) %>%
  bind_rows(.id = "Condition") %>%
  mutate(is_idp = ensembl_gene_id %in% idps) %>%
  ggplot(aes(x = log2FoldChange, fill = is_idp)) +
  geom_histogram(bins = 50) +
  scale_fill_manual(values=c("#feff7f", "#807fff")) +
  facet_wrap(~ Condition)

ggsave("significant plots/lfc_histograms_with_idps.png", width = 10, height = 6)