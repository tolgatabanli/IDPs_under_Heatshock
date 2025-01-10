library(tidyverse)
library(UpSetR)
library(gprofiler2)
library(ggpubr)

idps <- readRDS("IDP decisions/commons_modes.Rds")
deseq_results <- readRDS("deseq_results.Rds") %>% map(as.data.frame)
expr_direction <- deseq_results %>%
  map(drop_na) %>%
  map(arrange, padj) %>%
  map(~ dplyr::select(.x, log2FoldChange, padj)) %>%
  map(~ filter(.x, padj < 0.05)) %>%
  map(~ mutate(.x, direction = case_when(
    log2FoldChange < 0 ~ "down",
    log2FoldChange > 0 ~ "up"
  )))

plots <- list()
for (i in seq_along(expr_direction)) {
  plots[[i]] <- expr_direction[[i]] %>%
    ggplot(aes(x = direction, fill = direction)) +
    geom_bar() +
    scale_fill_manual(values = c("down" = "#00bfc4",
                                 "up" = "#f9766e"))
}
ggarrange(plotlist = plots, ncol = 2, nrow = 8, common.legend = T)

  

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


non_wt <- deseq_sig_idps %>% `[`(!startsWith(names(deseq_sig_idps), "Wildtype"))
wt <- deseq_sig_idps %>% `[`(startsWith(names(deseq_sig_idps), "Wildtype"))

reduce(wt, intersect) %>% setdiff(reduce(non_wt, intersect))

reduce(deseq_sig_idps, intersect) %>%
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


