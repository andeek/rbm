stats <- function(H, V, type="negative") {
  if(type == "binary") t <- data.frame(c(0,1))[,rep(1, H + V)]
  else if(type == "negative") t <- data.frame(c(-1,1))[,rep(1, H + V)]
  else stop("Type must be 'binary', 'negative'")
  
  names(t) <- c(paste0("h", 1:H), paste0("v", 1:V))  
  
  #all possibilities of statistics in (H + V)-space
  t.grid <- expand.grid(t)
  
  for(i in 1:H){
    for(j in (H+1):(H+V)){
      t.grid <- cbind(t.grid, t.grid[,i]*t.grid[,j])
      names(t.grid)[ncol(t.grid)] <- paste0("theta", i, j - H)
    }
  }
  
  return(data.matrix(t.grid))
}

expected_value <- function(theta, stats, normalized = TRUE) {
  result <- crossprod(t(crossprod(stats, exp(crossprod(t(stats), theta)))), diag(1/apply(exp(crossprod(t(stats), theta)), 2, sum)))
  rownames(result) <- paste0("exp_", rownames(result))
  return(result)
}

sample_sphere <- function(stat, n, r = 1) {
  require(dplyr)
  require(tidyr)
  
  samp <- function(space, r) {
    vec <- rnorm(space)
    paste0(r*vec/sqrt(sum(vec^2)), collapse=",")
  }
  
  res <- data.frame(samp = 1:n) %>%
    group_by(samp) %>%
    mutate(vals = samp(ncol(stat), r)) %>%
    separate(col = vals, sep = ",", into = colnames(stat)) %>%
    select(-samp) %>%
    data.matrix() %>%
    t()
    
  
  return(res)  
}

calc_hull <- function(stats, options="FA") {
  #calculate convex hull
  require(geometry)
  C <- convhulln(stats, options=options)
  return(C)
}

generate_hull_from_stats <- function(stats) {
  ## stat is the df of possible statistics as generated from stats
  lines <- NULL
  for(i in 1:nrow(stats)){
    for(j in 1:nrow(stats)){
      if(i < j) {
        lines <- c(lines, list(list(x1=stats[i,], x2=stats[j,])))
      }
    }
  }
  
  return(lines)
}

min_distance_to_hull <- function(x0, lines) {
  ## x0 is a single point (numeric vector)
  ## lines is a list of lines connecting all the points in the statistic space as generated from generate_hull_from_stats

  x0 <- unlist(x0)
  data.frame(min_dist = min(unlist(lapply(lines, function(line){
    x1 <- as.numeric(line$x1)
    x2 <- as.numeric(line$x2)
    
    rej <- (x1 - x0) - (sum((x1 - x0)*(x2-x1))/sum((x2-x1)^2))*(x2 - x1)
    
    sqrt(sum(rej^2))    
  }))))
}

project_to_sphere <- function(sphere, point) {
  ## sphere is a named list with two arguments: center and radius
  r <- sphere$radius
  center <- sphere$center  
  return((point - center)*r/sqrt(sum((point-center)^2)) + center  )
}

distance_to_simplex <- function(simplex, P) {
  require(dplyr)
  require(tidyr)
  ## P is a point in H + V + H*V space
  ## simplex a list of three points that make a triangle 
  ## algorithm implemented from http://data-science-radio.com/pubs/2012-eccad-distance-to-simplex.pdf
  
  d <- length(simplex) - 1

  simplex_copy <- simplex
  P_copy <- P

  P <- as.numeric(strsplit(P, ",")[[1]])
  simplex <- simplex %>% unlist %>%
    data.frame() %>%
    rename(point=simplex.....unlist) %>%
    separate(point, paste0("coord", 1:length(P)), sep=",") %>%
    data.matrix()
  
  if(d == 0) return(dist(matrix(c(P, simplex), nrow=2, byrow=TRUE)))
  
  ## translate to that S0 is the origin
  P <- P - simplex[1,]
  simplex <- sweep(simplex, 2, simplex[1,], "-")[-1,]
  
  ## project onto the linear subspace
  rows <- expand.grid(1:nrow(simplex), 1:nrow(simplex))
  alpha <- solve(matrix(apply(rows, 1, function(x){ crossprod(simplex[x[1],], simplex[x[2],])}), nrow=nrow(simplex)), apply(simplex, 1, crossprod, P))
  P0 <- colSums(alpha*simplex)
  
  if(sum(alpha) <= 1 & sum(alpha < 0) == 0) {
    return(dist(matrix(c(P, P0), nrow=2, byrow=TRUE)))
  } else if(sum(alpha < 0) > 0) {
    simplex0 <- simplex_copy[c(1, 1+ which(alpha > 0))]
    distance_to_simplex(P_copy, simplex0)    
  } else {
    distance_to_simplex(P_copy, simplex_copy[-1]) 
  }
  
}

