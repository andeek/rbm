#libraries -----------------------
library(dplyr)
library(tidyr)
library(ggplot2)

#load data -------------------------
load("../degeneracy/written/results_grid.RData")

#sample --------------------------------
#samples the parameters to get one set from the "degenerate" and one from "non-degenerate"

set.seed(1001) #reproducible 

res %>%
  filter(H == 4 & V == 4) %>%
  group_by(H, V, r1, r2) %>%
  do(data.frame(.$outside)) %>%
  left_join(res %>%
              filter(H == 4 & V == 4) %>%
              group_by(H, V, r1, r2) %>%
              do(data.frame(.$g_theta))) %>%
  ungroup() %>%
  group_by(near_hull) %>%
  sample_n(1) %>%
  select(H, V, r1, r2, near_hull, starts_with("V")) -> sample.params

names(sample.params)[-(1:5)] <- colnames(res$stat[[nrow(res)]])

save(sample.params, file = "written/sample.params.Rdata")

