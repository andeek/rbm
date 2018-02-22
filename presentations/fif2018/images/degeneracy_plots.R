library(knitr)
library(ggplot2)
library(dplyr)
library(tidyr)

opts_chunk$set(echo=FALSE, message=FALSE, warnings=FALSE)
theme_set(theme_bw(base_family="serif"))

source("../../../writing/resources/code/functions.R")

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
min_Q <- function(theta, stats) {
  apply(crossprod(t(stats), theta), 2, min)
}
min_max_Q <- function(theta, stats) {
  require(dplyr)
  tcrossprod(stats, theta) %>% data.frame() %>%
    cbind(stats) %>%
    group_by_(.dots = colnames(stats)[grepl("v", colnames(stats)) & !grepl("theta", colnames(stats))]) %>%
    summarise_each(funs(max), contains("X")) %>%
    ungroup() %>%
    summarise_each(funs(min), contains("X")) %>%
    select(contains("X")) %>%
    data.matrix() %>%
    t()
}
elpr <- function(theta, stats) {
  require(dplyr)
  exp(tcrossprod(stats, theta)) %>% data.frame() %>%
    cbind(stats) %>%
    group_by_(.dots = colnames(stats)[grepl("v", colnames(stats)) & !grepl("theta", colnames(stats))]) %>%
    summarise_each(funs(sum), contains("X")) %>%
    ungroup() -> marginalized
  
  marginalized %>%
    summarise_each(funs(max), contains("X")) %>%
    select(contains("X")) %>%
    data.matrix() -> max_marg
  
  marginalized %>%
    summarise_each(funs(min), contains("X")) %>%
    select(contains("X")) %>%
    data.matrix() -> min_marg
  
  t(log(max_marg/min_marg))
}


#grid data
load("../../../writing/resources/data/results_grid.RData")

#near-degeneracy
plot_dat_grid <- res %>% plot_data(grid = TRUE)

plot_dat_grid %>%
  group_by(r1, r2, H, V) %>%
  summarise(frac_degen = sum(near_hull)/n(), count = n()) %>% 
  ungroup() %>% 
  mutate(Hiddens = paste0("paste(n[h], '= ", H, "')"), Visibles = paste0("paste(n[v], '= ", V, "')")) -> convex_hull_summary

#uninterpretability
res %>%
  group_by(H, V, N, r1, r2) %>%
  do(indep_exp = t(expected_value(t(indep_params(.$samp[[1]], .$H, .$V)), .$stat[[1]]))[, -((.$H + .$V + 1):(.$H + .$V + .$H*.$V))],
     marg_exp = .$g_theta[[1]][, (.$H + .$V + .$H*.$V + 1):(ncol(.$g_theta[[1]])-.$H*.$V)]) -> exp_vals

exp_vals %>%
  group_by(H, V, N, r1, r2) %>%
  do(data.frame(max_abs_diff = apply(abs(.$indep_exp[[1]] - .$marg_exp[[1]]), 1, max))) %>%
  group_by(H, V, N, r1, r2) %>%
  summarise(max_abs_diff = mean(max_abs_diff)) %>%
  mutate(Hiddens = paste0("paste(n[h], '= ", H, "')"), Visibles = paste0("paste(n[v], '= ", V, "')")) -> exp_vals_summary

#instability
res %>%
  group_by(H, V, N, r1, r2) %>%
  do(data.frame(max = max_Q(t(.$samp[[1]]), .$stat[[1]]),
                min = min_Q(t(.$samp[[1]]), .$stat[[1]]),
                min_max = min_max_Q(.$samp[[1]], .$stat[[1]]))) %>%
  ungroup() %>% 
  group_by(H, V, N, r1, r2) %>%
  mutate(LHS1 = (max - min)/V,
         LHS2 = (max - min_max - H*log(2))/V) %>%
  summarise_each(funs(mean), LHS1, LHS2) -> max_q_summary

res %>%
  group_by(H, V, N, r1, r2) %>%
  do(data.frame(elpr = elpr(.$samp[[1]], .$stat[[1]]))) %>%
  ungroup() %>%
  group_by(H, V, N, r1, r2) %>% 
  summarise(mean_elpr = mean(elpr)) %>%
  mutate(scaled_mean_elpr = mean_elpr/V) -> elpr_summary


convex_hull_summary %>%
  left_join(elpr_summary) %>%
  left_join(exp_vals_summary) -> three_ways

three_ways$Visibles <- factor(three_ways$Visibles, levels=rev(unique(three_ways$Visibles)))

three_ways %>%
  mutate(r1 = round(r1/(H + V), 8), r2 = round(r2/(H * V), 8)) %>% 
  ggplot() +
  geom_tile(aes(x = r1, y = r2, fill = frac_degen)) +
  geom_contour(aes(x = r1, y = r2, z = frac_degen), colour = "black", bins = 8) +
  geom_contour(aes(x = r1, y = r2, z = frac_degen), colour = "black", breaks = .05, size = 1.5) +
  geom_abline(aes(intercept = 0, slope = 1), alpha = .5, lty = 2) +
  scale_fill_gradient("Fraction near-degenerate", low = "yellow", high = "red") +
  facet_grid(Visibles~Hiddens, labeller = label_parsed) +
  xlab(expression(paste(paste("||", theta[main], "||", "/(", n[h]," + ", n[v], ")")))) +
  ylab(expression(paste(paste("||", theta[interaction], "||", "/(", n[h], "*", n[v], ")")))) +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  ggtitle("Near-degeneracy") +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    plot.title=element_text(size=22)) -> p.degen

three_ways %>%
  mutate(r1 = round(r1/(H + V), 8), r2 = round(r2/(H * V), 8)) %>% 
  ggplot() +
  geom_tile(aes(x = r1, y = r2, fill = max_abs_diff)) +
  geom_contour(aes(x = r1, y = r2, z = max_abs_diff), colour = "black", bins = 8) +
  geom_abline(aes(intercept = 0, slope = 1), alpha = .5, lty = 2) +
  scale_fill_gradient("Mean max absolute difference", low = "yellow", high = "red", limits = c(0,2)) +
  facet_grid(Visibles~Hiddens, labeller = label_parsed) +
  xlab(expression(paste(paste("||", theta[main], "||", "/(", n[h]," + ", n[v], ")")))) +
  ylab(expression(paste(paste("||", theta[interaction], "||", "/(", n[h], "*", n[v], ")")))) +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  ggtitle("Uninterpretability") +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    plot.title=element_text(size=22)) -> p.exp_diff

three_ways %>%
  mutate(r1 = round(r1/(H + V), 8), r2 = round(r2/(H * V), 8)) %>% 
  ggplot() +
  geom_tile(aes(x = r1, y = r2, fill = scaled_mean_elpr)) +
  geom_contour(aes(x = r1, y = r2, z = scaled_mean_elpr), colour = "black", bins = 8) +
  geom_abline(aes(intercept = 0, slope = 1), alpha = .5, lty = 2) +
  scale_fill_gradient(expression(paste("Mean ", frac(LREP(theta), V))), low = "yellow", high = "red") +
  facet_grid(Visibles~Hiddens, labeller = label_parsed) +
  xlab(expression(paste(paste("||", theta[main], "||", "/(", n[h]," + ", n[v], ")")))) +
  ylab(expression(paste(paste("||", theta[interaction], "||", "/(", n[h], "*", n[v], ")")))) +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  ggtitle("Instability") +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    plot.title=element_text(size=22)) -> p.elpr


ggsave("degeneracy.pdf",
       plot = p.degen,
       bg = "transparent",
       width = 4.5,
       height = 5.5,
       units = "in")

ggsave("uninterpretability.pdf",
       plot =p.exp_diff,
       bg = "transparent",
       width = 4.5,
       height = 5.5,
       units = "in")

ggsave("instability.pdf",
       plot =p.elpr,
       bg = "transparent",
       width = 4.5,
       height = 5.5,
       units = "in")
