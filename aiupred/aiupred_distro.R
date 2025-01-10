library(reticulate)
library(tidyverse)

source_python("aiupred/read_pickle.py")
pickle_data <- read_pickle("aiupred/aiupred_scores.p")

all_disorder_scores <- pickle_data %>%
  reduce(c)

all_disorder_scores %>%
  as.tibble() %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 101, fill = "#1f77b4") +
  xlab("Local Disorder Score") +
  ylab("") +
  labs(title = "Distribution of local disorder scores") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5))
