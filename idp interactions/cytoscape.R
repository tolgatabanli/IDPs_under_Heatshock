library(tidyverse)
library(RCy3)
library(igraph)


data <- read_tsv("idp interactions/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.4.241.tab3.txt")
data <- data %>%
  filter(`Organism Name Interactor A` == "Saccharomyces cerevisiae (S288c)" &
           `Organism Name Interactor B` == "Saccharomyces cerevisiae (S288c)") %>%
  select(`Systematic Name Interactor A`, `Systematic Name Interactor B`,
         `Experimental System Type`, Modification) %>%
  distinct()

ppi <- data %>%
  filter(`Experimental System Type` == "physical") %>%
  select(starts_with("Systematic"))


common_sig_idps <- readRDS("common_sig_idps.Rds")
ppi_graph <- graph_from_data_frame(ppi)
neighs <- neighbors(ppi_graph, common_sig_idps)
induced_from_idps <- induced_subgraph(ppi_graph, common_sig_idps)

cytoscapePing()
createNetworkFromIgraph(induced_from_idps)
setNodeColorMapping(table.column="name",
                    table.column.values = common_sig_idps,
                    colors = "#FF0000",
                    mapping.type = "d")
