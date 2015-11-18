#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
source("ml_functs.R")


#data and params ------------------------
H <- 4
V <- 4
load("written/fitted_ml_margin.Rdata")
load("written/sample_images.Rdata")
load("written/sample.params.Rdata")

parms <- list(visibles = flat_images_good$visibles, 
              hiddens = flat_images_good$hiddens,
              H = H, V = V)

theta_good <- sample.params %>% filter(!near_hull) %>% ungroup() %>% select(starts_with("v"), starts_with("h"), starts_with("theta"), -H, -V) %>% data.matrix
theta_ml <- cbind(theta, likelihood) %>% data.frame()
names(theta_ml) <- c(colnames(theta_good), "likelihood")

#plots--------------------------
theta_ml %>%
  mutate(iter = 1:nrow(theta)) %>%
  gather(variable, value, -iter) %>%
  ggplot() +
  geom_point(aes(x = iter, y = value)) +
  facet_wrap(~variable, scales = "free_y") 

plot_two <- function(idx1, idx2, theta_true, theta_ml, parms) {

  expand.grid(seq(min(theta_ml[, idx1]), max(theta_ml[, idx1]), length.out = 30), 
              seq(min(theta_ml[, idx2]), max(theta_ml[, idx2]), length.out = 30)) %>% 
    data.frame() -> grid
  
  names(grid) <- colnames(theta_true)[c(idx1, idx2)]
  
  matrix(theta_true, nrow = nrow(grid), ncol = length(theta_true), byrow = TRUE) %>% 
    data.frame() -> theta_true_expand
  names(theta_true_expand) <- colnames(theta_true)
  
  theta_true_expand[, colnames(theta_true)[c(idx1, idx2)]] <- grid
  theta_true_expand$likelihood <- apply(theta_true_expand, 1, function(x) loglik(theta = x, parms = parms))
  
  theta_true_expand %>%
    ggplot(aes_string(x = colnames(theta_true)[idx1], y = colnames(theta_true)[idx2], z = "likelihood")) +
    stat_contour(aes(colour = ..level..)) +
    geom_path(data = theta_ml, aes(colour = likelihood)) +
    theme_bw() +
    theme(legend.position = "none")
}

plots <- list()
for(i in 1:(H + V + H*V)) {
  for(j in 1:(H + V + H*V)) {
    if(i < j) {
        plots[[i + (H + V + H*V)*(j - 1)]] <- plot_two(i, j, theta_good, theta_ml, parms)
        names(plots)[i + (H + V + H*V)*(j - 1)] <- paste(names(theta_ml)[i], names(theta_ml)[j], sep = "-")
    }
  }
}

save(plots, file = "written/ml_plots_margin.RData")

view_grid <- function(idx1, idx2) {
  expand.grid(names(theta_ml)[idx1], names(theta_ml)[idx2]) %>%
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
    filter(Var1 != Var2) %>%
    mutate(names = paste(Var1, Var2, sep = "-")) %>%
    select(names) %>%
    data.frame -> idx
  
  grid.arrange(grobs = plots[names(plots) %in% idx$names & !is.na(names(plots))])
}

## Visibles-Visibles
view_grid(1:V, 1:V)

## Visibles-Hiddens
view_grid(1:V, (V+1):(V+H))

## Visibles-Interaction
view_grid(1:V, (V+H+1):(V+H+H*V))

## Hiddens-Hiddens
view_grid((V+1):(V+H), (V+1):(V+H))

## Hiddens-Interactions
view_grid((V+1):(V+H), (V+H+1):(V+H+H*V))

## Interactions-Interactions
view_grid((V+H+1):(V+H+H*V), (V+H+1):(V+H+H*V))





