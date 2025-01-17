library(tidyverse)
library(UpSetR)
library(gprofiler2)
library(ggpubr)

# PREPARE

convert_to_counts <- function(df) {
  n_total <- df %>%
    group_by(direction) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(group = paste0("Total ", direction))
  
  n_idp <- df %>%
    filter(idp == TRUE) %>%
    group_by(direction) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(group = paste0("IDP ", direction))
  
  bind_rows(n_total, n_idp)
}

idps <- readRDS("IDP decisions/commons_modes.Rds")
deseq_results <- readRDS("deseq_results.Rds") %>% map(as.data.frame)
expr_direction <- deseq_results %>%
  map(drop_na) %>%
  map(arrange, padj) %>%
  map(~ dplyr::select(.x, log2FoldChange, padj)) %>%
  map(~ filter(.x, padj < 0.05)) %>%
  map(~ mutate(.x, direction = case_when(
    log2FoldChange < 0 ~ "Down",
    log2FoldChange > 0 ~ "Up"
  ))) %>%
  map(rownames_to_column) %>%
  map(~ mutate(.x, idp = case_when(
    rowname %in% idps ~ TRUE,
    !(rowname %in% idps) ~ FALSE
  ))) %>%
  map_dfr(convert_to_counts, .id = "dataset") %>%
  mutate(group = factor(group, levels = c("Total Up", "Total Down",
                                          "IDP Up", "IDP Down")))

ggplot(expr_direction, aes(x = group, y = count, fill = group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~dataset, scales = "free") +
  labs(
    title = "Comparison of Total and IDP Regulation",
    x = "Direction of Regulation in Total vs IDP",
    y = "Count",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c(
    "IDP Up" = "blue",
    "Total Up" = "lightblue",
    "IDP Down" = "red",
    "Total Down" = "pink"
  ))


ggsave("significant plots/expr_all_vs_idps.png", device = "png", height = 15, width = 20)


deseq_sig_idps <- deseq_results %>%
  map(drop_na) %>%
  map(arrange, padj) %>%
  map(~ rownames(.x[.x$padj < 0.05 & rownames(.x) %in% idps, ]))

all_sig <- purrr::reduce(deseq_sig_idps, union)
common_sig <- purrr::reduce(deseq_sig_idps, intersect)

idp_with_p <- deseq_results %>%
  map(drop_na) %>%
  map(select, "padj") %>%
  imap(~ rename(.x, !!paste0("padj_", .y) := padj)) %>%
  map(rownames_to_column) %>%
  purrr::reduce(inner_join, by = "rowname") %>%
  column_to_rownames()
  

dev.off()
png("significant plots/upset(ordered_degree)_from_deseq_on_commons_modes.png", width = 1920, height = 1080, res = 150)
upset(fromList(deseq_sig_idps), sets = names(deseq_sig_idps),
      keep.order = T,
      nsets = length(deseq_sig_idps),
      point.size = 3, line.size = 1,
      mainbar.y.label = "Significant IDP Intersection",
      sets.x.label = "Number of IDPs ",
      order.by = c("degree"),
      text.scale = c(1.7,1,1.3,1,1,1.3))
dev.off()

gost_list <- list()

for (sigs in deseq_sig_idps) {
  gost(query = sigs,
       organism = "scerevisiae",
       correction_method = "bonferroni",
       ordered_query = F)
  break
}

purrr::reduce(deseq_sig_idps, dplyr::intersect) %>%
  gost(query = ., organism = "scerevisiae",
       correction_method = "bonferroni", ordered_query = F) %>%
  gostplot()

gost_multi_query <- gost(query = deseq_sig_idps, multi_query = T,
     organism = "scerevisiae",
     correction_method = "bonferroni")

gost(all_sig,
     ordered_query = T,
     organism = "scerevisiae",
     correction_method = "bonferroni") %>%
  gostplot()


