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

###############################################
## Plots
###############################################

## Plot the progression of "worst" parameter
s <- scatterplot3d(min_points[, 1:3], pch=19)
segs <- NULL
for(i in 1:(nrow(min_points)-1)) {
  start <- s$xyz.convert(min_points[i, 1:3])
  end <- s$xyz.convert(min_points[i + 1, 1:3])
  segs <- rbind(segs, c(start$x, start$y, end$x, end$y))
}
segments(segs[,1], segs[,2], segs[,3], segs[,4], lty=2)
## Doesn't look like any sort of pattern

## Plot of min_dist by radius
qplot(r, min_dist, colour = min_dist < 1e-4, data=min_points) + geom_hline(aes(yintercept=1e-4))
