#libraries -----------------------
library(dplyr)
library(tidyr)
library(ggplot2)

#load data -------------------------
load("../degeneracy/written/results_grid.RData")

#find appropriate radii and correlate with variance constants
res %>%
  filter(H == 4 & V == 4) %>%
  group_by(H, V, r1, r2) %>%
  do(data.frame(.$outside)) %>%
  group_by(H, V, r1, r2) %>%
  summarise(percent_degen = sum(near_hull)/n()) %>%
  filter(percent_degen <= .05 & r1 > r2) %>%
  mutate(C = (r1/3)^2, C_prime = (r2/3)^2) %>%
  ungroup() %>%
  select(percent_degen, C, C_prime) %>%
  data.frame -> variance_params

save(variance_params, file = "written/variance_params.Rdata")
  