source("functions.R")
library(dplyr)
library(reshape2)
library(ggplot2)

## constants
n <- 1000 #num points to sample on the sphere

min_radius <- function(H, V, type='negative', n) {
  stat <- stats(H, V, type)
  lines <- generate_hull_from_stats(stat)
  
  ## inits
  r0 <- 0
  theta0 <- sample_sphere(space = ncol(stat), n, r0)
  exp_vals0 <- expected_value(stats = stat, theta = theta0, normalized = TRUE)
  d0 <- min(apply(exp_vals0, 2, function(x) min_distance_to_hull(x, lines)))
  step0 <- 
    
  search <- function(inits, lines, cutoff) {
    r <- inits$r0
    d <- inits$d0
    step <- inits$step0
    
    if(d < cutoff) return(r)
    
    r1 <- r + step
    theta1 <- sample_sphere(space = ncol(stat), n, r1)
    exp_vals1 <- expected_value(stats = stat, theta = theta1, normalized = TRUE)
    d1 <- min(apply(exp_vals1, 2, function(x) min_distance_to_hull(x, lines)))
    
    r2 <- r1 - step/2
    theta2 <- sample_sphere(space = ncol(stat), n, r2)
    exp_vals2 <- expected_value(stats = stat, theta = theta2, normalized = TRUE)
    d2 <- min(apply(exp_vals2, 2, function(x) min_distance_to_hull(x, lines)))
  
    if(d1 >= cutoff) {
      search(list(r0=r1, d0=d1, step0=step), lines, cutoff)
    } else {
      search(list(r0=r2, d0=d2, step0=step/2), lines, cutoff)
    } 
  }
  
  search(list(r0=r0, d0=d0, step0=step0), lines, .001)
}


