#libraries --------------------
library(dplyr)
library(tidyr)

loglik_num <- function(theta, parms) {
  #theta is a H + V + H*V length vector containing main_visible, then main_hidden, then interactions
  
  #parms is a named list containing data and H and V
  visibles <- parms$visibles
  hiddens <- parms$hiddens
  H <- parms$H
  V <- parms$V
  
  #things for derivative, build data and possible values
  N <- visibles %>% nrow()
  data <- cbind(visibles, hiddens)
  for(i in 1:V) {
    for(j in 1:H) {
      data <- cbind(data, visibles[i,]*hiddens[j,])
    }
  }
  
  sum(theta*colSums(data))
}


loglik <- function(theta, parms) {
  #theta is a H + V + H*V length vector containing main_visible, then main_hidden, then interactions
  
  #parms is a named list containing data and H and V
  visibles <- parms$visibles %>% data.frame()
  H <- parms$H
  V <- parms$V
  names(visibles) <- paste0("v", 1:V)
  
  #things for derivative, build data and possible values
  N <- visibles %>% nrow()
  possibles <- stats(H, V) 
  W <- exp(possibles %*% theta)
  
  # marginalize out the hiddens
  data <- visibles %>%
    left_join(possibles %>% data.frame(), by = paste0("v", 1:V))

  data %>%
    group_by_(.dots = paste0("h", 1:H)) %>%
    do(parms = list(visibles = data.frame(.) %>% select(starts_with("v")) %>% data.matrix(),
                    hiddens = data.frame(.) %>% select(starts_with("h")) %>% data.matrix(),
                    H = H,
                    V = V)) %>%
    mutate(A = loglik_num(theta = theta, parms)) %>%
    mutate(num = exp(A)) %>%
    ungroup() %>%
    summarise(numerator = sum(num)) %>%
    as.numeric() -> numerator
  
    log(numerator) - N*log(sum(W))
}

#not working
loglik_derivs <- function(theta, parms) {
  #theta is a H + V + H*V length vector containing main_visible, then main_hidden, then interactions
  
  #parms is a named list containing data and H and V
  visibles <- parms$visibles %>% data.frame()
  H <- parms$H
  V <- parms$V
  names(visibles) <- paste0("v", 1:V)
  
  #things for derivative, build data and possible values
  N <- visibles %>% nrow()
  possibles <- stats(H, V) 
  
  # marginalize out the hiddens
  data <- visibles %>%
    left_join(possibles %>% data.frame(), by = paste0("v", 1:V))

  #vector of derivatives set to zero
  W <- exp(possibles %*% theta)

  data %>%
    group_by_(.dots = paste0("h", 1:H)) %>%
    summarise_each(funs(sum), everything()) - 
    colSums(apply(possibles, 2, function(x) { W*x }))/sum(W)*N -> tmp
  
  colSums(tmp)

}

loglik_hessian <- function(theta, parms) {
  #theta is a H + V + H*V length vector containing main_visible, then main_hidden, then interactions
  
  #parms is a named list containing data and H and V
  visibles <- parms$visibles
  hiddens <- parms$hiddens
  H <- parms$H
  V <- parms$V
  
  #things for derivative, build data and possible values
  N <- visibles %>% nrow()
  data <- cbind(visibles, hiddens)
  for(i in 1:V) {
    for(j in 1:H) {
      data <- cbind(data, visibles[i,]*hiddens[j,])
    }
  }
  possibles <- stats(H, V)
  
  W <- exp(possibles %*% theta)
  
  res <- array(NA, dim = c(nrow(possibles), ncol(possibles), ncol(possibles)))
  for(i in 1:ncol(possibles)) {
    for(j in 1:ncol(possibles)) {
      res[, i, j] <- possibles[, i]*possibles[, j]*W
    }
  }
  res_sum <- matrix(0, nrow = ncol(possibles), ncol = ncol(possibles))
  for(i in 1:nrow(possibles)) {
    res_sum <- res_sum + res[i, , ]
  }
  first_term <- res_sum/sum(W)
  
  second_term <- matrix(NA, nrow = ncol(possibles), ncol = ncol(possibles))
  for(i in 1:ncol(possibles)) {
    for(j in 1:ncol(possibles)) {
      second_term[i, j] <- sum(possibles[, i]*W)*sum(possibles[, j]*W)/(sum(W))^2
    }
  }
  
  second_term - first_term
}

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