##function for creating all possible values of the t(x)
stats <- function(H, V, type="negative") {
  if(type == "binary") t <- data.frame(c(0,1))[,rep(1, H + V)]
  else if(type == "negative") t <- data.frame(c(-1,1))[,rep(1, H + V)]
  else stop("Type must be 'binary', 'negative'")
  
  names(t) <- c(paste0("h", 1:H), paste0("v", 1:V))  
  
  #all possibilities of statistics in (H + V)-space
  t.grid <- expand.grid(t)
  
  for(i in 1:V){
    for(j in (V+1):(H+V)){
      t.grid <- cbind(t.grid, t.grid[,i]*t.grid[,j])
      names(t.grid)[ncol(t.grid)] <- paste0("theta", i, j - H)
    }
  }
  
  return(data.matrix(t.grid))
}

##function for calculating the convex hull of t(x)
calc_hull <- function(H, V, type="binary") {
  #only possible values of t are 0 or 1
  if(type == "binary") t <- data.frame(c(0,1))[,rep(1, H + V)]
  else if(type == "negative") t <- data.frame(c(-1,1))[,rep(1, H + V)]
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
  
  #calculate convex hull
  library(geometry)
  C <- convhulln(t.grid, options="FA")
  return(list(possible_t=t.grid, c_hull=C))
}