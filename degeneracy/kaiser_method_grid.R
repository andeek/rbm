#libraries -----------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("functions.R")

#load sample params ---------------------------------------------
n <- 100
load("apps/results_grid.RData") 

indep_params <- function(samp, H, V) {
  samp[, (H + V + 1):ncol(samp)] <- 0
  samp
}

res %>%
  group_by(H, V, N, r1, r2) %>%
  do(indep_exp = t(expected_value(t(indep_params(.$samp[[1]], .$H, .$V)), .$stat[[1]]))[, -((.$H + .$V + 1):(.$H + .$V + .$H*.$V))],
     marg_exp = .$g_theta[[1]][, (.$H + .$V + .$H*.$V + 1):(ncol(.$g_theta[[1]])-.$H*.$V)]) -> exp_vals

exp_vals %>%
  group_by(H, V, N, r1, r2) %>%
  do(data.frame(mean_abs_diff = apply(abs(.$indep_exp[[1]] - .$marg_exp[[1]]), 1, mean))) %>%
  group_by(H, V, N, r1, r2) %>%
  summarise(mean_abs_diff = mean(mean_abs_diff)) %>%
  mutate(Hiddens = paste0("Hiddens: ", H), Visibles = paste0("Visibles: ", V)) %>%
  ggplot() +
  geom_tile(aes(x = r1, y = r2, fill = mean_abs_diff)) +
  geom_contour(aes(x = r1, y = r2, z = mean_abs_diff), colour = "black", bins = 8) +
  #geom_contour(aes(x = r1, y = r2, z = mean_abs_diff), colour = "black", breaks = .1, size = 1.5) +
  #geom_abline(aes(intercept = 0, slope = 1), alpha = .5, lty = 2) +
  scale_fill_gradient("Mean absolute difference in expectation", low = "yellow", high = "red") +
  facet_grid(Visibles~Hiddens) +
  xlab(expression(group("||", theta[main], "||"))) +
  ylab(expression(group("||", theta[interaction], "||"))) +
  theme(aspect.ratio = 1)

exp_vals %>% 
  mutate(dist = sqrt(r1^2 + r2^2)) %>%
  group_by(dist, H, V) %>%
  do(data.frame(mean_abs_diff = mean(abs(.$indep_exp[[1]] - .$marg_exp[[1]])))) %>%
  mutate(Hiddens = as.character(H), Visibles = paste0("Visibles: ", V)) %>%
  ggplot() +
  geom_point(aes(dist, mean_abs_diff, colour = Hiddens)) + ylim(c(0, .0025))
  #geom_smooth(aes(dist, mean_abs_diff, colour = Hiddens, group = Hiddens)) +
  facet_grid(~Visibles) +
  xlab(expression(paste(group("||", theta, "||"), ", Distance from the origin"))) +
  ylab("Mean absolute difference in expectation")