## load libraries --------------------------- 
source("functions.R")
library(dplyr)
library(tidyr)
library(ggplot2)

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

sample_sphere_main <- function(stat, n, r = 1) {
  require(dplyr)
  require(tidyr)

  samp <- function(space, r) {
    vec <- rnorm(space)
    paste0(r*vec/sqrt(sum(vec^2)), collapse=",")
  }

  res <- data.frame(samp = 1:n) %>%
    group_by(samp) %>%
    mutate(vals = samp(ncol(stat %>% data.frame() %>% select(-contains("theta"))), r)) %>%
    separate(col = vals, sep = ",", into = colnames(stat %>% data.frame() %>% select(-contains("theta")))) %>%
    select(-samp) %>%
    data.matrix() %>%
    t()

  return(res)
}

sample_sphere_interact <- function(stat, n, r = 1) {
  require(dplyr)
  require(tidyr)

  samp <- function(space, r) {
    vec <- rnorm(space)
    paste0(r*vec/sqrt(sum(vec^2)), collapse=",")
  }

  res <- data.frame(samp = 1:n) %>%
    group_by(samp) %>%
    mutate(vals = samp(ncol(stat %>% data.frame() %>% select(contains("theta"))), r)) %>%
    separate(col = vals, sep = ",", into = colnames(stat %>% data.frame() %>% select(contains("theta")))) %>%
    select(-samp) %>%
    data.matrix() %>%
    t()

  return(res)
}

## params ----------------------------------
n <- 100
r_center <- seq(0.2, 3, by = 0.2)
#r_center <- 1
r_exp <- .05
power <- seq(0, 1, by = .1)

## function for running simulation -------------------------------
get_res <- function(n, r_center, r_exp, power) {
  ## perform functions
  test_cases <- expand.grid(data.frame(1:4)[,rep(1,2)]) %>%
    rename(H = X1.4, V = X1.4.1) %>%
    filter(H <= V) %>% #remove those cases with less visibles than hiddens. Is this necessary?
    mutate(n_param = H*V + H + V) %>%
    mutate(max_facets = (2^(H+V))^(floor(n_param/2))) %>%
    mutate(n = n, r1 = r_center*(H + V), r2 = r_center*(H + V)^(power)) %>%
    group_by(H, V, n_param, n, r1, r2) %>%
    do(stat = stats(.$H, .$V, "negative")) %>%
    ungroup() %>%
    group_by(H, V, n_param, n, r1, r2) %>%
    do(samp = rbind(sample_sphere_main(.$stat[[1]], .$n, .$r1), sample_sphere_interact(.$stat[[1]], .$n, .$r2)))

  test_cases %>%
    left_join(test_cases %>%
                group_by(H, V, n_param, n, r1, r2) %>%
                do(true_r = apply(.$samp[[1]], 2, function(x) sqrt(sum(x^2))))) -> test_cases


  test_cases_stat <- expand.grid(data.frame(1:4)[,rep(1,2)]) %>%
    rename(H = X1.4, V = X1.4.1) %>%
    filter(H <= V) %>% #remove those cases with less visibles than hiddens. Is this necessary?
    mutate(n_param = H*V + H + V) %>%
    mutate(max_facets = (2^(H+V))^(floor(n_param/2))) %>%
    mutate(n = n, r1 = r_center*(H + V), r2 = r_center*(H + V)^(power)) %>%
    group_by(H, V, n_param, n, r1, r2) %>%
    do(stat = stats(.$H, .$V, "negative"))

  test_cases <- inner_join(test_cases, test_cases_stat)

  test_cases %>%
    group_by(H, V, n_param, n, r1, r2) %>%
    do(g_theta = t(rbind(.$samp[[1]], expected_value(theta = .$samp[[1]], stats = data.matrix(.$stat[[1]]))))) -> cases

  cases <- inner_join(cases, test_cases)

  cases %>%
    group_by(H, V, r1, r2) %>%
    do(outside = find_prop(.$g_theta[[1]] %>% data.frame() %>% select(starts_with("exp")), .$stat[[1]], r_exp)) -> tmp

  res <- inner_join(tmp, cases %>% group_by(H, V, r1, r2, n_param) %>% do(samp = .$g_theta[[1]]))
}

## run simulations -------------------------------------
for(i in power) {
  cat(paste("power:", i))
  res <- lapply(r_center, function(r) get_res(n, r, r_exp, i))
  save(res, file = paste0("written/results_split_", i, ".RData"))
}
