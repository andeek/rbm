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
epsilon <- seq(-1, 1, length.out = 11)

expand.grid(H = 1:4, V = 1:4, C = C, epsilon = epsilon) %>%
  mutate(N = H + V) %>%
  mutate(n = n, r1 = sqrt(C), r2 = sqrt(C*(N/(H*V))^(epsilon))) %>%
  group_by(H, V, C, epsilon, r1, r2, N) %>%
  do(samp = cbind(matrix(rnorm(n*.$N, mean = 0, sd = .$r1/3), nrow = n), matrix(rnorm(n*.$H*.$V, mean = 0, sd = .$r2/3), nrow = n))) -> split_sample

expand.grid(H = 1:4, V = 1:4, C = C, epsilon = epsilon) %>%
  mutate(N = H + V) %>%
  mutate(n = n, r1 = sqrt(C), r2 = sqrt(C*(N/(H*V))^(epsilon))) %>%
  group_by(H, V, C, epsilon, r1, r2, N) %>%
  do(stat = stats(.$H, .$V, "negative")) -> split_sample_stat

split_sample <- inner_join(split_sample, split_sample_stat)

split_sample %>%
  group_by(H, V, N, C, epsilon, r1, r2) %>%
  do(g_theta = t(rbind(t(.$samp[[1]]), expected_value(theta = t(.$samp[[1]]), stats = data.matrix(.$stat[[1]]))))) -> split

split <- inner_join(split, split_sample)

split %>%
  group_by(H, V, N, C, epsilon, r1, r2) %>%
  do(outside = find_prop(.$g_theta[[1]] %>% data.frame() %>% select(starts_with("exp")), .$stat[[1]], r_exp)) -> tmp

res <- inner_join(tmp, split)
save(res, file = "written/results_split.RData")
