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

sample_sphere_unif <- function(space, rep, center, r) {
  vec <- matrix(rnorm(space*rep), nrow = rep)
  r*vec/t(matrix(rep(sqrt(rowSums(vec^2)), each = space), nrow = space))
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
r_exp <- .05

expand.grid(H = 1:4, V = 1:4, r1 = seq(0.001, 3, length.out = 20), r2 = seq(0.001, 3, length.out = 20)) %>%
  mutate(r1 = r1*(H + V), r2 = r2*(H*V)) %>%
  mutate(N = H + V) %>%
  group_by(H, V, r1, r2, N) %>%
  do(samp = cbind(sample_sphere_unif(.$N, n, 0, .$r1), sample_sphere_unif(.$H*.$V, n, 0, .$r2))) -> grid_sample

expand.grid(H = 1:4, V = 1:4, r1 = seq(0.001, 3, length.out = 20), r2 = seq(0.001, 3, length.out = 20)) %>%
  mutate(r1 = r1*(H + V), r2 = r2*(H*V)) %>%
  mutate(N = H + V) %>%
  group_by(H, V, r1, r2, N) %>%
  do(stat = stats(.$H, .$V, "negative")) -> grid_sample_stat

grid_sample <- inner_join(grid_sample, grid_sample_stat)

grid_sample %>%
  group_by(H, V, N, r1, r2) %>%
  do(g_theta = t(rbind(t(.$samp[[1]]), expected_value(theta = t(.$samp[[1]]), stats = data.matrix(.$stat[[1]]))))) -> grid

grid <- inner_join(grid, grid_sample)

grid %>% 
  ungroup() %>%
  group_by(H, V, N, r1, r2) %>%
  mutate(epsilon = (1-(1-2*r_exp)^(3/(H + V + H*V)))/2) %>%
  do(outside = find_prop(.$g_theta[[1]] %>% data.frame() %>% select(starts_with("exp")), .$stat[[1]], .$epsilon)) -> tmp

res <- inner_join(tmp, grid)
save(res, file = "written/results_grid.RData")
