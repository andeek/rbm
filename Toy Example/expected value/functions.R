
stats <- function(H, V, type="negative") {
  if(type == "binary") t <- data.frame(c(0,1))[,rep(1, H + V)]
  else if(type == "negative") t <- data.frame(c(-1,1))[,rep(1, H + V)]
  else if(type == "experiment") t <- data.frame(c(1,3))[,rep(1, H + V)]
  else stop("Type must be 'binary', 'negative', or 'experiment'")
  
  names(t) <- c(paste0("h", 1:H), paste0("v", 1:V))  
  
  #all possibilities of statistics in (H + V)-space
  t.grid <- expand.grid(t)
  
  for(i in 1:H){
    for(j in (H+1):(H+V)){
      t.grid <- cbind(t.grid, t.grid[,i]*t.grid[,j])
      names(t.grid)[ncol(t.grid)] <- paste0("theta", i, j - H)
    }
  }
  
  return(t.grid)
}

expected_value <- function(stats, theta, normalized = TRUE) {
  if(normalized) result <- t(as.matrix(stats)) %*% exp(as.matrix(stats) %*% as.matrix(theta)) %*% diag(1/apply(exp(as.matrix(stats) %*% as.matrix(theta)), 2, sum))
  else result <- t(as.matrix(stats)) %*% exp(as.matrix(stats) %*% as.matrix(theta))
  
  rownames(result) <- paste0("exp_", rownames(result))
  return(result)
}

sample_sphere <- function(space, n, r = 1) {
  return(data.frame(sapply(1:n, function(i){
    vec <- rnorm(space)
    r*vec/sqrt(sum(vec^2))
  })))  
}

calc_hull <- function(stats) {
  #calculate convex hull
  require(geometry)
  C <- convhulln(stats, options="FA")
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
  
  min(unlist(lapply(lines, function(line){
    x1 <- as.numeric(line$x1)
    x2 <- as.numeric(line$x2)
    
    rej <- (x1 - x0) - (sum((x1 - x0)*(x2-x1))/sum((x2-x1)^2))*(x2 - x1)
    
    sqrt(sum(rej^2))    
  })))
}

project_to_sphere <- function(sphere, point) {
  ## sphere is a named list with two arguments: center and radius
  r <- sphere$radius
  center <- sphere$center  
  return((point - center)*r/sqrt(sum((point-center)^2)) + center  )
}

