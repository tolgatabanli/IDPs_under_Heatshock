library(tidyverse)
library(UpSetR)
library(gprofiler2)
library(ggpubr)
library(janitor)

# PREPARE

idps <- readRDS("IDP decisions/commons_modes.Rds")

## Get DESeq Results
### Reference25
deseq_results <- readRDS("deseq_results.Rds") %>% map(as.data.frame)
ref <- 25
### Reference37
deseq_results <- readRDS("deseq_reference37.Rds") %>%
  map(as.data.frame)
ref <- 37

deseq_sig <- deseq_results %>%
  map(drop_na) %>%
  map(arrange, padj) %>%
  map(~ dplyr::select(.x, log2FoldChange, padj)) %>%
  map(~ filter(.x, padj < 0.05)) %>%
  map(~ filter(.x, abs(log2FoldChange) > 0.6)) # for lfc filtering

expr_direction <- deseq_sig %>%
  map(~ mutate(.x, direction = case_when(
    log2FoldChange < 0 ~ "Down",
    log2FoldChange > 0 ~ "Up"
  ))) %>%
  map(rownames_to_column) %>%
  map(~ mutate(.x, idp = case_when(
    rowname %in% idps ~ TRUE,
    !(rowname %in% idps) ~ FALSE
  )))

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

# Expression Barplot
expr_direction %>%
  map_dfr(convert_to_counts, .id = "dataset") %>%
  # filter(str_starts(dataset, "Wildtype")) %>%
  mutate(group = factor(group, levels = c("Total Up", "Total Down",
                                          "IDP Up", "IDP Down"))) %>%
  ggplot(aes(x = group, y = count, fill = group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ dataset, ncol = 2, scales = "free",
             labeller = labeller(dataset = function(x) {
               gsub("(\\w+?)_(\\d+?)_(\\d+?)", "\\1 \\2 °C, t = \\3", x)
             })) +
  labs(
    title = "Comparison of Total and IDP Regulation",
    x = "Direction of Regulation of Total vs IDP under heat shock",
    caption = paste("Reference is", ref, "°C of corresponding time and genotype."),
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
ggsave(paste0("significant plots/expr_reference", ref, ".png"),
       device = "png", height = 12, width = 7)
ggsave(paste0("significant plots/expr_reference", ref, "_wildtype.png"),
       device = "png", height = 3, width = 7)

# Wildtype early vs late IDPs
updown_wt <- expr_direction %>%
  keep(str_starts(names(.), "Wildtype")) %>%
  map(filter, idp) %>%
  map(~ mutate(.x, direction = factor(.x$direction, levels = c("Up", "Down")))) %>%
  map(~ split(.x$rowname, .x$direction)) %>%
  flatten() %>%
  set_names(c("Wildtype_42_10 Up", "Wildtype_42_10 Down",
              "Wildtype_42_30 Up", "Wildtype_42_30 Down"))

dev.off()
png(paste0("significant plots/upset(freq)_ref",ref,"_wildtype.png"),
    width = 1920, height = 1080, res = 200)
upset(fromList(updown_wt), sets = names(updown_wt),
      nsets = length(updown_wt),
      keep.order = T,
      point.size = 3, line.size = 1,
      mainbar.y.label = "Significant IDP Intersection",
      sets.x.label = "Number of IDPs",
      order.by = c("freq"),
      text.scale = c(1.7, 1, 1.3, 1, 1.3, 1.3))
dev.off()

deseq_sig_all <- deseq_sig %>%
  map(row.names)
deseq_sig_idps <- deseq_sig_all %>%
  map(intersect, idps)

# UPSET plot
dev.off()
png(paste0("significant plots/upset(freq)_ref",ref,".png"),
    width = 1920, height = 1080, res = 150)
upset(fromList(deseq_sig_idps), sets = names(deseq_sig_idps),
      keep.order = T,
      nsets = length(deseq_sig_idps),
      point.size = 3, line.size = 1,
      mainbar.y.label = "Significant IDP Intersection",
      sets.x.label = "Number of IDPs",
      order.by = c("freq"),
      text.scale = c(1.7,1,1.3,1,1,1.3))
dev.off()



# Gene Enrichment
dotplot_updown <- function(sub_deseq, reg) {
  xdirection <- sub_deseq %>%
    filter(direction == reg & idp) %>%
    pull(rowname)
  multi_gp = gost(query = xdirection,
                  organism = "scerevisiae",
                  multi_query = FALSE, evcodes = TRUE)
  gp_mod = multi_gp$result[,c("query", "source", "term_id",
                              "term_name", "p_value", "query_size", 
                              "intersection_size", "term_size", 
                              "effective_domain_size", "intersection")]
  gp_mod$GeneRatio = gp_mod$intersection_size / gp_mod$query_size
  gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
  names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                    "query_size", "Count", "term_size", "effective_domain_size", 
                    "geneID", "GeneRatio", "BgRatio")
  gp_mod$geneID = gsub(",", "/", gp_mod$geneID)
  return(gp_mod)
}

dot_data <- list()
expr_names <- names(expr_direction)
for (cond in seq_along(expr_direction)*2) {
  label <- expr_names[[cond/2]]
  dot_data[[cond - 1]] <- dotplot_updown(expr_direction[[cond/2]], "Up") %>%
    mutate(Cluster = paste0(label, "_", "Up"))
  dot_data[[cond]] <- dotplot_updown(expr_direction[[cond/2]], "Down") %>%
    mutate(Cluster = paste0(label, "_", "Down"))
  print(cond)
}

dot_data %>%
  map(filter, Category %in% c("GO:MF", "GO:BP", "GO:CC")) %>%
  map(arrange, term_size, desc(Count)) %>%
  map(slice_head, n = 10) %>%
  bind_rows() %>%
  filter(str_starts(Cluster, "Wildtype")) %>%
  ggplot(aes(x = Cluster, y = Description, size = Count, color = p.adjust)) +
  geom_point() +
  scale_size_continuous(range = c(3, 12)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)) +
  labs(x = "Cluster", y = "Enrichment Term", title = "Gene Enrichment Analysis") +
  facet_wrap(Category ~ ., scales = "free_y", ncol = 1)
#ggsave("significant plots/dotplot_reference37.jpeg", height = 20, width=15)
ggsave("significant plots/dotplot_reference37_wildtypes.jpeg", height = 12, width=12)





# Reference 37
# Gosts of Upregulated and Downregulated !!! CHANGE: Up / Down
deseq_sig_idps <- expr_direction %>%
  map(filter, rowname %in% idps) %>%
  map(pull, rowname)

xregulated_idps <- expr_direction %>%
  map(filter, direction == "Down" & idp)

gost_res <- xregulated_idps %>%
  janitor::clean_names() %>%
  map(pull, rowname) %>%
  map(~ gost(query = .x,
             organism = "scerevisiae",
             correction_method = "fdr"))

conditions <- gost_res %>%
  map(~ `$`(.x, "result") %>%
        as.data.frame() %>%
        arrange(p_value) %>%
        select(term_id, term_name, term_size, intersection_size, p_value) %>%
        filter(term_size < 50) %>%
        arrange(desc(intersection_size), term_size, p_value) %>%
        dplyr::slice(1:10)
  )

query_sizes <- xregulated_idps %>%
  map(nrow) %>% unlist() %>% unname()

write(html_content, file = "gosttables/downregulated_idps_reference37.html")










# miscellaneuous

# Against whole genome
gost_res <- deseq_sig_idps %>%
  janitor::clean_names() %>%
  map(~ gost(query = .x,
             organism = "scerevisiae",
             correction_method = "fdr"))

# GostTable
library(kableExtra)

generate_collapsible_panels <- function(conditions, query_sizes) {
  html_panels <- lapply(seq_along(conditions), function(i) {
    condition_name <- names(conditions)[i]
    query_size <- query_sizes[i]
    condition_table <- conditions[[i]] %>%
      select(term_id, term_name, p_value, term_size, intersection_size) %>%
      kbl(
        format = "html", 
        col.names = c("Term ID", "Term Name", "p-value", "Term Size", "Intersection Size")
      ) %>%
      kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed"), font_size = 12)
    
    # Create collapsible panel HTML
    paste0(
      "<div class='panel panel-default'>",
      "<div class='panel-heading'>",
      "<h4 class='panel-title'>",
      "<a data-toggle='collapse' href='#collapse", i, "'>",
      condition_name, " with ", query_size, " queries",
      "</a>",
      "</h4>",
      "</div>",
      "<div id='collapse", i, "' class='panel-collapse collapse'>",
      "<div class='panel-body'>", condition_table, "</div>",
      "</div>",
      "</div>"
    )
  })
  paste0(html_panels, collapse = "\n")
}

html_content <- paste0(
  "<!DOCTYPE html>
  <html lang='en'>
  <head>
    <meta charset='UTF-8'>
    <title>Condition Comparison</title>
    <link rel='stylesheet' href='https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css'>
    <script src='https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js'></script>
    <script src='https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js'></script>
  </head>
  <body>
  <div class='container'>
  <h1>Comparison of Conditions</h1>",
  generate_collapsible_panels(conditions, query_sizes),
  "</div>
  </body>
  </html>"
)
