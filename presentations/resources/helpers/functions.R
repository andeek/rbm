##function for creating all possible values of the t(x)
stats <- function(H, V, type="negative") {
  if(type == "binary") t <- data.frame(c(0,1))[,rep(1, H + V)]
  else if(type == "negative") t <- data.frame(c(-1,1))[,rep(1, H + V)]
  else stop("Type must be 'binary', 'negative'")
  
  names(t) <- c(paste0("v", 1:H), paste0("h", 1:V))  
  
  #all possibilities of statistics in (H + V)-space
  t.grid <- expand.grid(t)
  
  for(i in 1:V){
    for(j in (V+1):(H+V)){
      t.grid <- cbind(t.grid, t.grid[,i]*t.grid[,j])
      names(t.grid)[ncol(t.grid)] <- paste0("theta", i, j - V)
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

expected_value <- function(theta, stats, normalized = TRUE) {
  result <- crossprod(t(crossprod(stats, exp(crossprod(t(stats), theta)))), diag(1/apply(exp(crossprod(t(stats), theta)), 2, sum), nrow=length(apply(exp(crossprod(t(stats), theta)), 2, sum))))
  rownames(result) <- paste0("exp_", rownames(result))
  return(result)
}

visible_distn <- function(params) {
  #params is a list containing three named elements of initial parameter values: main_hidden, main_visible, and interaction
  H <- length(params$main_hidden)
  V <- length(params$main_visible)
  
  theta <- matrix(c(params$main_visible, params$main_hidden, as.numeric(t(params$interaction))))
  possibles <- stats(H, V) 
  
  e_to_the_junk <- exp(possibles %*% theta)
  
  data.frame(possibles, prob = e_to_the_junk/sum(e_to_the_junk)) %>% 
    group_by_(.dots = paste0("v", 1:V)) %>%
    summarise(prob = sum(prob)) %>%
    ungroup() %>%
    mutate(image_id = 1:(2^V)) %>%
    data.frame()
}
