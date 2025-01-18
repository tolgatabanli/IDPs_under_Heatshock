library(tidyverse)
library(UpSetR)
library(gprofiler2)
library(ggpubr)
library(janitor)

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



  

dev.off()
png("significant plots/upset(ordered_freq)_from_deseq_on_commons_modes.png", width = 1920, height = 1080, res = 150)
upset(fromList(deseq_sig_idps), sets = names(deseq_sig_idps),
      keep.order = T,
      nsets = length(deseq_sig_idps),
      point.size = 3, line.size = 1,
      mainbar.y.label = "Significant IDP Intersection",
      sets.x.label = "Number of IDPs ",
      order.by = c("freq"),
      text.scale = c(1.7,1,1.3,1,1,1.3))
dev.off()



set_base_url("https://biit.cs.ut.ee/gprofiler_beta") # for when the server's down

# Against whole genome
gost_res <- deseq_sig_idps %>%
  janitor::clean_names() %>%
  map(~ gost(query = .x,
             organism = "scerevisiae",
             correction_method = "bonferroni"))

# GostTable
gost_res %>%
  map(`$`(result) %>%
        as.data.frame() %>%
        arrange(p_value) %>%
        select(term_id, term_name, p_value) %>%
        dplyr::slice(1:10)) %>%
  iwalk(~ publish_gosttable(gostres = .x,
                           show_columns = c("term_name"),
                           filename = paste0("gosttables/", .y, ".png")))

# GostPlots
terms <- gost_res %>%
  map(~ `$`(.x, "result")) %>%
  map(dplyr::pull, term_id) %>%
  purrr::reduce(base::union)

gost_res %>%
  map(~ {
    .x$result <- arrange(.x$result, p_value) %>%
      dplyr::slice(1:10)
    print(.x)
    .x
  }) %>%
  map(~ gostplot(.x, interactive = F)) %>%
  iwalk(~ publish_gostplot(p = .x,
                           filename = paste0("gostplots/", .y, ".png"),
                           highlight_terms = terms))

# Against IDPs: no result
gost_res <- deseq_sig_idps %>%
  janitor::clean_names() %>%
  map(~ gost(query = .x,
             organism = "scerevisiae",
             custom_bg = idps,
             correction_method = "bonferroni"))

gost_res %>%
  compact() %>%
  map(~ gostplot(.x, interactive = F)) %>%
  iwalk(~ publish_gostplot(p = .x,
                           filename = paste0("gostplots_idp_bg/", .y, ".png")))

# Wildtype 42_30 minus Double_42_30
wt_heatshock <- base::setdiff(deseq_sig_idps[["Wildtype_42_10"]],
                              deseq_sig_idps[["Wildtype_37_10"]])
exp_heatshock <- base::setdiff(deseq_sig_idps[["MSN24.KO_42_10"]],
                               deseq_sig_idps[["MSN24.KO_37_10"]])
diff <- base::setdiff(wt_heatshock, exp_heatshock)

diff %>%
  gost(query = .,
       organism = "scerevisiae",
       correction_method = "bonferroni") %>%
  `$`(result) %>%
  arrange(p_value) %>%
  dplyr::slice(1:10) %>%
  publish_gosttable(show_columns = c("term_name"))

