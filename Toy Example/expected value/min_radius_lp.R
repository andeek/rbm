source("functions.R")
library(dplyr)
library(tidyr)
library(lpSolveAPI)

## constants
n <- 500 #num points to sample on the sphere
cutoff <- .05
step <- .5
backstep <- .05
v <- seq(1, 3, by=.01)

test_cases <- expand.grid(data.frame(1:4)[,rep(1,2)]) %>% 
  rename(H = X1.4, V = X1.4.1) %>%
  filter(H <= V) %>% #remove those cases with less visibles than hiddens. Is this necessary?
  mutate(n_param = H*V + H + V) %>%
  #filter(n_param <= 11) %>%
  rowwise() %>%
  mutate(min_rad = min_radius_lp(H, V, n=n, cutoff=cutoff, backstep=backstep, step=step, v=v))

save(test_cases, file="min_radius_lp.Rdata")


min_radius_lp <- function(H, V, type = "negative", n, cutoff, backstep, step, v) {
  stat <- stats(H, V, type)
  r0 <- 0
  theta0 <- sample_sphere(stat, n, r0)
  exp_vals0 <- expected_value(stats = stat, theta = theta0, normalized = TRUE)
  d0 <- 1
  
  success <- FALSE
  while(!success){  
    r1 <- r0 + step
    theta1 <- sample_sphere(stat, n, r1)
    exp_vals1 <- expected_value(stats = stat, theta = theta1, normalized = TRUE)
    
#     which_v <- v %>%
#       data.frame() %>%
#       rowwise() %>%
#       do(exp_v = t(exp_vals1*.$v)) %>%    
#       do(num = sum(apply(.$exp_v, 1, function(y) !in_hull(stat, y)))) %>%
#       summarise(num=as.numeric(num))
#     
    outside <- FALSE
    i <- 1
    while(!outside) {
      num_outside <- sum(apply(exp_vals1*v[i], 2, function(y) !in_hull(stat, y)))
      i <- i + 1
      if(num_outside > 0) outside <- TRUE
    }
    
    exp_vals1_v <- exp_vals1*v[i]    
    d1 <- min(sapply(1:nrow(t(exp_vals1)), function(x) dist(matrix(c(exp_vals1[,x], exp_vals1_v[,x]), nrow=2, byrow = TRUE))))
    
    if(d0 < cutoff & d1 >= cutoff) success <- TRUE
    
    if(d1 < cutoff) {
      r0 <- r1
      d0 <- d1
      step <- -backstep
    } else {
      r0 <- r1
      d0 <- d1
    }
    cat(paste0("n_dim: ", H + V + H*V, " r: ", r0, " d:", d0, "\n"))   
  }
  return(r0)
}



in_hull <- function(hull_points, point) {
  require(lpSolveAPI)
  P <- data.matrix(point)
  A <- t(data.matrix(hull_points))
  
  lp_obj <- make.lp(nrow = 0, ncol = ncol(A))
  sapply(1:nrow(A), function(x) add.constraint(lp_obj, A[x,], type = "=", rhs = P[x]))
  add.constraint(lp_obj, rep(1, ncol(A)), type = "=", rhs = 1)
  sapply(1:ncol(A), function(x) add.constraint(lp_obj, 1, type = ">=", rhs = 0, indices = x))
  
  return(solve(lp_obj) == 0)
}


