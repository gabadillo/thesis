rm(list = ls())

library(tidyverse)

# summarize cluster results
phenotype.df = read_csv('../data/clean/simulated_metrics.csv')

path = '../data/output/cluster'
csv_files <- list.files(path = path, pattern = "\\.csv$", full.names = TRUE)
output.df <- map_df(csv_files, read_csv) %>% arrange(Y.approach, Reps)

h2s = c(1,0.8, 0.5, 0.2)
phenotypes = unique(phenotype.df$Phenotype)

Scenarios = list( # which observations are used to compute the metrics?
  all = 1:6,
  firsts = 1:3,
  mids = 3:5,
  evens = c(2,4,6),
  odds = c(1,3,5),
  limits = c(1,2,6),
  pair1 = c(2,5),
  pair2 = c(3,6),
  first = 1,
  fourth = 4,
  last = 6
)

output.df %>%
  ggplot(aes(x = Y.approach,
             y = abs(b.Value)))+
  geom_boxplot()+
  facet_grid(h2 ~ factor(Available,levels=names(Scenarios)))

#write_csv(output.df, '../data/output/scenario3.csv')

