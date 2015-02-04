source("functions.R")
library(dplyr)
library(tidyr)

n <- 100
r_center <- 2

test_cases <- expand.grid(data.frame(1:4)[,rep(1,2)]) %>% 
  rename(H = X1.4, V = X1.4.1) %>%
  filter(H <= V) %>% #remove those cases with less visibles than hiddens. Is this necessary?
  mutate(n_param = H*V + H + V) %>%
  mutate(max_facets = (2^(H+V))^(floor(n_param/2))) %>%
  mutate(n = n, r = r_center) %>%
  group_by(H, V, n_param, n, r) %>%
  do(stat = stats(.$H, .$V, "negative")) %>%
  ungroup() %>%
  group_by(H, V, n_param, n, r) %>%
  do(g_theta = t(expected_value(theta = sample_sphere(.$stat[[1]], .$n, .$r), stats = data.matrix(.$stat[[1]]))))

test_cases_stat <- expand.grid(data.frame(1:4)[,rep(1,2)]) %>% 
  rename(H = X1.4, V = X1.4.1) %>%
  filter(H <= V) %>% #remove those cases with less visibles than hiddens. Is this necessary?
  mutate(n_param = H*V + H + V) %>%
  mutate(max_facets = (2^(H+V))^(floor(n_param/2))) %>%
  mutate(n = n, r = r_center) %>%
  group_by(H, V, n_param, n, r) %>%
  do(stat = stats(.$H, .$V, "negative"))

cases <- inner_join(test_cases, test_cases_stat)

cases %>%
  group_by(H, V) %>%
  do(min_dist = find_min_dist(.$g_theta[[1]], .$stat[[1]], .05, .005, .$H + .$V)) %>%
  mutate(min_dist = as.numeric(min_dist))-> dists

find_min_dist <- function(g_theta, stat, step, backstep, n_param) {
  r_exp_0 <- 0
  num_0 <- 0
  
  success <- FALSE
  while(!success) {
    r_exp_1 <- r_exp_0 + step
    
    for(i in 1:nrow(g_theta)) {
      g_theta[i,] %>%
        sample_points_outside(stat, r_exp_1) %>%
        filter(!in_hull) %>%
        dim() -> num_1
      
      if(num_1[1] > 0) break
    }
    

    if(num_1[1] == 0 & num_0 > 0) success <- TRUE
    
    r_exp_0 <- r_exp_1
    num_0 <- num_1[1]
    
    if(num_1[1] > 0) {
      step <- -backstep
    }
  } 
  
  return(data.frame(r = r_exp_0, num_outside = num_0))
}

sample_points_outside <- function(g_theta, stats, r, n = 100) {
  #g_theta is the expected value mapping of a point theta in R^H+V+H*V
  x0 <- data.frame(samp = 1:n) %>%
        mutate(vals = samp_exp(length(g_theta), g_theta, r)) %>%
        separate(col = vals, sep = ",", into = colnames(stats)) %>%
        select(-samp)
     
  x0$in_hull <- apply(x0, 1, in_hull, stats)

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