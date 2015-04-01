#libraries -----------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("functions.R")

## main function --------------------------------------------
find_prop <- function(g_theta, stat, r_exp) {
  g_theta %>%
    data.frame() %>%
    group_by_(.dots = c(colnames(g_theta))) %>%
    do(samp = sample_points_outside(as.numeric(.), stat, r_exp)$in_hull) %>%
    group_by_(.dots = c(colnames(g_theta))) %>%
    do(as.data.frame(sum(.$samp[[1]]) > 0)) -> is_outside # count of not in hull greater than zero => there exists a point outside hull for this model
  
  is_outside %>% 
    rename_(near_hull = as.name(names(is_outside)[ncol(is_outside)])) %>%
    mutate(r = r_exp) -> inside_outside
  
  return(inside_outside)
}

## additional useful functions -------------------------------------
sample_points_outside <- function(g_theta, stats, r, n = 100) {
  #g_theta is the expected value mapping of a point theta in R^H+V+H*V
  x0 <- data.frame(samp = 1:n) %>%
    mutate(vals = samp_exp(length(g_theta), g_theta, r)) %>%
    separate(col = vals, sep = ",", into = colnames(stats)) %>%
    select(-samp)
  
  x0$in_hull <- !apply(x0, 1, in_hull, stats)
  
  return(x0)
}

samp_exp <- function(space, center, r) {
  vec <- rnorm(space)
  paste0(r*vec/sqrt(sum(vec^2)) + center, collapse=",")
}

in_hull <- function(point, hull_points) {
  require(lpSolveAPI)
  P <- data.matrix(point)
  A <- t(data.matrix(hull_points))
  
  lp_obj <- make.lp(nrow = 0, ncol = ncol(A))
  sapply(1:nrow(A), function(x) add.constraint(lp_obj, A[x,], type = "=", rhs = P[x]))
  add.constraint(lp_obj, rep(1, ncol(A)), type = "=", rhs = 1)
  sapply(1:ncol(A), function(x) add.constraint(lp_obj, 1, type = ">=", rhs = 0, indices = x))
  
  return(solve(lp_obj) == 0)
}

#sample split radii ---------------------------------------------
n <- 100
C <- seq(0.2, 3, by = 0.2)
r_exp <- .05
epsilon <- seq(0, .5, by = .05)

expand.grid(H = 1:4, V = 1:4, C = C, epsilon = epsilon) %>%
  mutate(N = H + V) %>%
  mutate(n = n, r1 = C, r2 = C*(N/(H*V))^(1+epsilon)) %>%
  group_by(H, V, C, epsilon, r1, r2, N) %>%
  do(samp = cbind(matrix(rnorm(n*.$N, mean = 0, sd = .$r1/3), nrow = n), matrix(rnorm(n*.$H*.$V, mean = 0, sd = .$r2/3), nrow = n))) -> split_sample

expand.grid(H = 1:4, V = 1:4, C = C, epsilon = epsilon) %>%
  mutate(N = H + V) %>%
  mutate(n = n, r1 = C, r2 = C*(N/(H*V))^(1+epsilon)) %>%
  group_by(H, V, C, epsilon, r1, r2, N) %>%
  do(stat = stats(.$H, .$V, "negative")) -> split_sample_stat

split_sample <- inner_join(split_sample, split_sample_stat)

split_sample %>%
  group_by(H, V, N, C, epsilon, r1, r2) %>%
  do(g_theta = t(rbind(t(.$samp[[1]]), expected_value(theta = t(.$samp[[1]]), stats = data.matrix(.$stat[[1]]))))) -> split

split <- inner_join(split, split_sample)

split %>%
  group_by(H, V, N,C, epsilon, r1, r2) %>%
  do(outside = find_prop(.$g_theta[[1]] %>% data.frame() %>% select(starts_with("exp")), .$stat[[1]], r_exp)) -> tmp

res <- inner_join(tmp, split) 
save(res, file = "apps/results_split.RData")


# plots ------------------------
plot_data <- function(res) {
  
  plot.data <- data.frame()
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
      mutate_(ss_ratio = paste0("(", paste(paste0(names(.)[(H + V + 1):(H + V + H*V)], "^2"), collapse = " + "), ")/sum(", paste(paste0(names(.)[1:(H+V)],"^2"), collapse = " + "), ")")) %>%    
      ungroup() -> ratio
    
    
    inner_join(ratio, res$outside[[i]] %>% ungroup()) %>%
      select(ss_ratio, near_hull) %>%
      mutate(H = H, V = V, n_param = H + V + H*V, N = H + V, N = N, r1 = r1, r2 = r2, C = C, epsilon = epsilon) %>%
      rbind(plot.data) -> plot.data
  }
  return(plot.data)
}

plot_dat <- res %>% plot_data

plot_dat %>%
  group_by(H, V, n_param, N, r1, r2, C, epsilon) %>%
  summarise(percent_degen = sum(near_hull)/n()) %>%
  ggplot() +
  geom_point(aes(n_param, percent_degen)) +
  facet_grid(epsilon~C)
  
plot_dat %>%
  filter(epsilon == 0) %>%
  ggplot() +
  geom_boxplot(aes(factor(n_param), ss_ratio, colour = near_hull)) + 
  facet_wrap(~C, scales = "free_y")

