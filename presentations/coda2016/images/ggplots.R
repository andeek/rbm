library(knitr)
library(ggplot2)
library(dplyr)
library(tidyr)

opts_chunk$set(echo=FALSE, message=FALSE, warnings=FALSE)
theme_set(theme_bw(base_family="serif"))

source("../../degeneracy/functions.R")

##variable encoding ----------------------
test_cases <- expand.grid(data.frame(1:4)[,rep(1,2)]) %>% 
  rename(H = X1.4, V = X1.4.1) %>%
  filter(H <= V) %>% #remove those cases with less visibles than hiddens. Is this necessary?
  mutate(n_param = H*V + H + V) %>%
  mutate(max_facets = (2^(H+V))^(floor(n_param/2))) %>%
  filter(n_param <= 11) #calc_hull can't handle any higher dimensions currently

bin_res <- plyr::dlply(test_cases, plyr::.(n_param), function(x) calc_hull(x$V, x$H, "binary"))
neg_res <- plyr::dlply(test_cases, plyr::.(n_param), function(x) calc_hull(x$V, x$H, "negative"))

plyr::ldply(bin_res, function(x) x$c_hull$vol) %>% 
  mutate(frac_vol = V1/(1^n_param)) %>%
  inner_join(plyr::ldply(neg_res, function(x) x$c_hull$vol) %>% 
               mutate(frac_vol = V1/(2^n_param)),
             by="n_param") %>%
  rename(vol.bin = V1.x, vol.neg = V1.y, frac_vol.bin = frac_vol.x, frac_vol.neg = frac_vol.y) %>%
  gather(vars, value, -n_param) %>%
  separate(vars, c("type", "encoding"), "\\.") %>%
  spread(type, value) %>%
  ggplot() +
  geom_point(aes(x=n_param, y=frac_vol, colour=encoding)) +
  geom_line(aes(x=n_param, y=frac_vol, colour=encoding, group=encoding)) +
  ylab("Fraction of unrestricted volume") +
  xlab("Number of parameters") +
  scale_colour_discrete("Encoding", labels=c("Binary (1,0)", "Negative (1,-1)")) +           
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)) -> p

ggsave("images/frac_volume.pdf",
       plot = p,
       bg = "transparent",
       width = 7,
       height = 3.9,
       units = "in")

#manageable examples ---------
#reshape data functions
plot_data <- function(res, grid = FALSE) {
  
  plot.data <- data.frame()
  
  if(!grid) {
    for(i in 1:nrow(res)) {
      tmp <- res$g_theta[[i]] %>% data.frame()
      H <- res[i,]$H
      V <- res[i,]$V
      N <- res[i,]$N
      r1 <- res[i,]$r1
      r2 <- res[i,]$r2
      C <- res[i,]$C
      epsilon <- res[i,]$epsilon
      
      
      tmp %>% 
        rowwise() %>% 
        mutate_(ss_interaction = paste0("(", paste(paste0(names(.)[(H + V + 1):(H + V + H*V)], "^2"), collapse = " + "), ")"),
                ss_main = paste0("sum(", paste(paste0(names(.)[1:(H+V)],"^2"), collapse = " + "), ")")) %>%    
        ungroup() -> ratio
      
      
      inner_join(ratio, res$outside[[i]] %>% ungroup()) %>%
        select(ss_interaction, ss_main, near_hull) %>%
        mutate(H = H, V = V, n_param = H + V + H*V, N = H + V, N = N, r1 = r1, r2 = r2, C = C, epsilon = epsilon) %>%
        rbind(plot.data) -> plot.data
    }
  } else {
    for(i in 1:nrow(res)) {
      tmp <- res$g_theta[[i]] %>% data.frame()
      H <- res[i,]$H
      V <- res[i,]$V
      N <- res[i,]$N
      r1 <- res[i,]$r1
      r2 <- res[i,]$r2
      
      tmp %>% 
        rowwise() %>% 
        mutate_(ss_interaction = paste0("(", paste(paste0(names(.)[(H + V + 1):(H + V + H*V)], "^2"), collapse = " + "), ")"),
                ss_main = paste0("sum(", paste(paste0(names(.)[1:(H+V)],"^2"), collapse = " + "), ")")) %>%  
        ungroup() -> ratio
      
      
      inner_join(ratio, res$outside[[i]] %>% ungroup()) %>%
        select(ss_interaction, ss_main, near_hull) %>%
        mutate(H = H, V = V, n_param = H + V + H*V, N = H + V, N = N, r1 = r1, r2 = r2) %>%
        rbind(plot.data) -> plot.data
    }
  }
  
  return(plot.data)
}
indep_params <- function(samp, H, V) {
  samp[, (H + V + 1):ncol(samp)] <- 0
  samp
}
max_Q <- function(theta, stats) {
  apply(crossprod(t(stats), theta), 2, max)
}

#grid data
load("../../writing/data/results_grid.RData")

#near-degeneracy
plot_dat_grid <- res %>% plot_data(grid = TRUE)

plot_dat_grid %>%
  group_by(r1, r2, H, V) %>%
  summarise(frac_degen = sum(near_hull)/n(), count = n()) %>% 
  ungroup() %>% 
  mutate(Hiddens = paste0("Hiddens: ", H), Visibles = paste0("Visibles: ", V)) -> convex_hull_summary

#uninterpretability
res %>%
  group_by(H, V, N, r1, r2) %>%
  do(indep_exp = t(expected_value(t(indep_params(.$samp[[1]], .$H, .$V)), .$stat[[1]]))[, -((.$H + .$V + 1):(.$H + .$V + .$H*.$V))],
     marg_exp = .$g_theta[[1]][, (.$H + .$V + .$H*.$V + 1):(ncol(.$g_theta[[1]])-.$H*.$V)]) -> exp_vals

exp_vals %>%
  group_by(H, V, N, r1, r2) %>%
  do(data.frame(mean_abs_diff = apply(abs(.$indep_exp[[1]] - .$marg_exp[[1]]), 1, mean))) %>%
  group_by(H, V, N, r1, r2) %>%
  summarise(mean_abs_diff = mean(mean_abs_diff)) %>%
  mutate(Hiddens = paste0("Hiddens: ", H), Visibles = paste0("Visibles: ", V)) -> exp_vals_summary

#instability
res %>%
  group_by(H, V, N, r1, r2) %>%
  do(data.frame(LHS = max_Q(t(.$samp[[1]]), .$stat[[1]])/(.$H + .$V))) %>%
  summarise(mean_max_q = mean(LHS)) -> max_q_summary


convex_hull_summary %>%
  left_join(max_q_summary) %>%
  left_join(exp_vals_summary) -> three_ways

three_ways %>%
  ggplot() +
  geom_tile(aes(x = r1, y = r2, fill = frac_degen)) +
  geom_contour(aes(x = r1, y = r2, z = frac_degen), colour = "black", bins = 8) +
  geom_contour(aes(x = r1, y = r2, z = frac_degen), colour = "black", breaks = .05, size = 1.5) +
  geom_abline(aes(intercept = 0, slope = 1), alpha = .5, lty = 2) +
  scale_fill_gradient("Fraction near-degenerate", low = "yellow", high = "red") +
  facet_grid(Visibles~Hiddens) +
  xlab(expression(group("||", theta[main], "||"))) +
  ylab(expression(group("||", theta[interaction], "||"))) +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  ggtitle("Near-degeneracy") +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)) -> p.degen

three_ways %>%
  ggplot() +
  geom_tile(aes(x = r1, y = r2, fill = mean_abs_diff)) +
  geom_contour(aes(x = r1, y = r2, z = mean_abs_diff), colour = "black", bins = 8) +
  #geom_contour(aes(x = r1, y = r2, z = mean_abs_diff), colour = "black", breaks = .05, size = 1.5) +
  geom_abline(aes(intercept = 0, slope = 1), alpha = .5, lty = 2) +
  scale_fill_gradient(expression(group("|", E(bold(X), "|", bold(theta)) - E(bold(X), "|", plain(independence)), "|")), low = "yellow", high = "red") +
  facet_grid(Visibles~Hiddens) +
  xlab(expression(group("||", theta[main], "||"))) +
  ylab(expression(group("||", theta[interaction], "||"))) +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  ggtitle("Uninterpretability") +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)) -> p.exp_diff

three_ways %>%
  ggplot() +
  geom_tile(aes(x = r1, y = r2, fill = mean_max_q)) +
  geom_contour(aes(x = r1, y = r2, z = mean_max_q), colour = "black", bins = 8) +
  #geom_contour(aes(x = r1, y = r2, z = mean_max_q), colour = "black", breaks = .05, size = 1.5) +
  geom_abline(aes(intercept = 0, slope = 1), alpha = .5, lty = 2) +
  scale_fill_gradient(expression(frac(max(Q), H + V)), low = "yellow", high = "red") +
  facet_grid(Visibles~Hiddens) +
  xlab(expression(group("||", theta[main], "||"))) +
  ylab(expression(group("||", theta[interaction], "||"))) +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  ggtitle("Instability") +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)) -> p.max_q

ggsave("images/degeneracy.pdf",
       plot = p.degen,
       bg = "transparent",
       width = 4.5,
       height = 5.25,
       units = "in")

ggsave("images/uninterpretability.pdf",
       plot =p.exp_diff,
       bg = "transparent",
       width = 4.5,
       height = 5.25,
       units = "in")

ggsave("images/instability.pdf",
       plot =p.max_q,
       bg = "transparent",
       width = 4.5,
       height = 5.25,
       units = "in")
