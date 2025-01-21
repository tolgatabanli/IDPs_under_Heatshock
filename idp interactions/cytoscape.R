library(tidyverse)
library(RCy3)


data <- read_tsv("idp interactions/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.4.241.tab3.txt")
ppi <- data %>%
  filter(`Organism Name Interactor A` == "Saccharomyces cerevisiae (S288c)" &
           `Organism Name Interactor B` == "Saccharomyces cerevisiae (S288c)") %>%
  select(`Systematic Name Interactor A`, `Systematic Name Interactor B`,
         `Experimental System Type`, Modification) %>%
  distinct()

common_sig_idps <- readRDS("common_sig_idps.Rds")

cytoscapePing()
ind_sub <- 

