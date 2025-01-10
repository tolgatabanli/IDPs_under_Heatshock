library(tidyverse)

local_plddts <- readRDS("IDP decisions/local_plddts.RDS")

# save in csv
all_scores <- local_plddts %>% select(scores) %>% pull() %>%
  reduce(c)
write.csv(all_scores, "alphafold/local_plddts.csv")
all_scores <- read.csv("alphafold/local_plddts.csv") %>% `colnames<-`(c("row_names", "locals"))
saveRDS(all_scores, "IDP decisions/all_local_plddt_scores.Rds")


# Hist of all local plddt scores
ggplot(all_scores, aes(locals)) +
  geom_histogram(bins = 100, fill = "#1f77b4") +
  xlab("Local pLDDT") +
  ylab("") +
  labs(title = "Distribution of all local pLDDTs from AlphaFold") +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14),
        title = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = margin(t = 20),
        aspect.ratio = 0.75) +
  scale_x_continuous(breaks = c(0,20,40,60,80,100)) +
  scale_y_continuous(expand = c(0,0,0.1,0)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("alphafold/hist (bins_100) of all plddt scores.png", width = 6.4, height = 4.8)

find_mode <- function(scores, bins) {
  hist_data <- hist(scores, breaks = bins, plot = FALSE)
  
  index_of_mode <- which.max(hist_data$counts)
  mode <- (hist_data$breaks[index_of_mode] + hist_data$breaks[index_of_mode + 1]) / 2
  return(mode)
}

find_mode_python <- function(scores, bins) {
  bin_edges <- seq(min(scores), max(scores), length.out = bins + 1)
  bin_counts <- hist(scores, breaks = bin_edges, plot = FALSE)$counts
  bin_centers <- (bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2
  mode_value <- bin_centers[which.max(bin_counts)]
  return(mode_value)
}

with_modes <- local_plddts %>% mutate(mode_plddt =
                          sapply(scores, find_mode_python, bins = 20)) %>%
  dplyr::select(uniprot, mode_plddt)

# Hist of modes
ggplot(with_modes, aes(x = mode_plddt)) +
  geom_histogram(bins = 50, fill = "#1f77b4") +
  xlab("Modes of proteins (pLDDT)") +
  ylab("") +
  labs(title = "Distribution of mode pLDDTs") +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14),
        title = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = margin(t = 30),
        aspect.ratio = 0.75) +
  scale_x_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,100)) +
  scale_y_continuous(expand = c(0,0,0.1,0)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("alphafold/hist (bins_50) of binned(20) mode plddts.png", width = 6.4, height = 4.8)

# Hist of means
means <- read_tsv("IDP decisions/all_proteins_from_alphafold_with_mean.tsv")

ggplot(means, aes(mean_plddt)) +
  geom_histogram(bins = 50, fill = "#1f77b4") +
  xlab("Mean pLDDT") +
  ylab("") +
  labs(title = "Distribution of mean pLDDTs") +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14),
        title = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = margin(t = 20),
        aspect.ratio = 0.75) +
  scale_x_continuous(breaks = c(0,20,40,60,80,100)) +
  scale_y_continuous(expand = c(0,0,0.1,0)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("alphafold/hist (bins_50) of mean plddt scores.png", width = 6.4, height = 4.8)

# Hist of perc60
perc <- read_tsv("IDP decisions/all_proteins_from_alphafold_with_perc_alpha60.tsv")
ggplot(perc, aes(alphafold_perc_60)) +
  geom_histogram(bins = 50, fill = "#1f77b4") +
  xlab("Percent disorder") +
  ylab("") +
  labs(title = expression(paste("Distribution of ", percent[alpha==60], " disorders with AlphaFold"))) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14),
        title = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = margin(t = 20),
        aspect.ratio = 0.75) +
  scale_y_continuous(expand = c(0,0,0.1,0)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("alphafold/hist (bins_50) of percent_alpha60 disorder alphafold.png", width = 6.4, height = 4.8)

# join with ensembl and save
ensembl_up_mapping <- readRDS("IDP decisions/ensembl_up_mapping.Rds")

# For all proteins
ensembl_up_mapping %>%
  inner_join(with_modes, join_by(uniprotswissprot == uniprot)) %>%
  dplyr::select(ensembl_gene_id, mode_plddt) %>%
  write_tsv("IDP decisions/all_proteins_from_alphafold_with_mode.tsv")


# IDPs, threshold
read_tsv("IDP decisions/all_proteins_from_alphafold_with_mode.tsv") %>%
  filter(mode_plddt < 60) %>%
  write_tsv("IDP decisions/idps_from_alphafold_mode.tsv")

