source("functions.R")
library(dplyr)
library(tidyr)

H <- 1
V <- 2

stat <- stats(H, V, "negative")
n <- 100
r_center <- 2

theta <- sample_sphere(stat, n, r_center)
g_theta <- t(expected_value(theta, stat))

r_exp <- .015

apply(g_theta, 1, function(x){
  x %>%
    sample_points_outside(stat, r_exp) %>%
    filter(!in_hull) %>%
    nrow() -> got_out
  res <- ifelse(got_out > 0, r_exp, NA)   
})




sample_points_outside <- function(g_theta, stats, r, n = 100) {
  #g_theta is the expected value mapping of a point theta in R^H+V+H*V
  x0 <- data.frame(samp = 1:n) %>%
        mutate(vals = samp_exp(length(g_theta), g_theta, r)) %>%
        separate(col = vals, sep = ",", into = colnames(stats)) %>%
        select(-samp)
     
  x0$in_hull <- apply(x0, 1, in_hull, stat)

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