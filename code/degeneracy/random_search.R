source("functions.R")
library(scatterplot3d)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)


###############################################
## Generate data
###############################################

## Num nodes
H <- 1
V <- 1

## Num to sample
n <- 10000

## Radius to sample
r <- c(seq(0,1, by=0.05), seq(1.5, 10, by=0.5))

stat <- stats(H, V, "negative")
theta <- mlply(r, function(x) sample_sphere(space = ncol(stat), n, x))
exp_vals <- llply(theta, function(x) expected_value(stats = stat, theta = x, normalized = TRUE))
lines <- generate_hull_from_stats(stat)

dist <- ldply(exp_vals, function(y) apply(y, 2, function(x){ min_distance_to_hull(x, lines) }))
min_dist <- apply(dist[,-1], 1, min)
min_point <- apply(dist[,-1], 1, which.min)

## Data frame of points at each radius with minimum distance to a closest line
min_points <- ldply(1:length(theta), function(i) theta[[i]][, min_point[i]])
names(min_points) <- names(stat)
min_points$r <- r
min_points$min_dist <- min_dist

## Start with first radius that has a point within .01 of a line
cut <- .01
sigma <- .5
n0 <- 500

mu0 <- min_points[min_dist < cut,][1, 1:ncol(stat)] 
sphere = list(center = c(0,0,0), radius = min(min_points[min_dist < cut, "r"]))


local_search <- function(mu0, sigma, sphere, dim = 3, n = 500, cutoff = .01) {
  store_mus <- rbind(mu0, c(0,0,0))
  
  search <- function(mu0, sigma, sphere, dim, n, cutoff) {
    local_theta <- t(matrix(rnorm(n*dim, sd = sigma), nrow=dim) + mu0)
    proj_theta <- t(apply(local_theta, 1, function(point) project_to_sphere(sphere=sphere, point=point)))
    exp_vals <- expected_value(stats = stat, theta = t(proj_theta), normalized = TRUE)
    dist <- apply(exp_vals, 2, function(x)  min_distance_to_hull(x, lines))
    mu1 <- proj_theta[which.min(dist),]
    store_mus <<- rbind(mu1, store_mus)
    if(sqrt(sum((mu1 - mu0)^2)) >= cutoff) {
      search(proj_theta[which.min(dist),], sigma, sphere, dim, n, cutoff) 
    } else {
      return(store_mus)
    }
  }
  search(mu0, sigma, sphere, dim, n, cutoff)
  return(store_mus)
}


mus <- lapply(1:100, function(i){
  cat("iteration = ", i, "\n")
  local_search(mu0, sigma, sphere, cutoff=.005)
})


###############################################
## Plots
###############################################

s <- scatterplot3d(mus[[1]], pch=19)
segs <- NULL
for(i in 1:(nrow(mus[[1]])-1)) {
  start <- s$xyz.convert(mus[[1]][i,])
  end <- s$xyz.convert(mus[[1]][i + 1,])
  segs <- rbind(segs, c(start$x, start$y, end$x, end$y))
}
segments(segs[,1], segs[,2], segs[,3], segs[,4], lty=2)


scatterplot3d(ldply(mus, function(x) x[1,]), pch=19)

