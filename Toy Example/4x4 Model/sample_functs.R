sample_main_hidden <- function(hiddens, C) {
  #hiddens is an NxH matrix of hidden node values
  #C is the variance constant
  res <- sapply(colSums(hiddens), function(mean) rnorm(1, mean*C, sqrt(C)))
  names(res) <- paste0("H", 1:ncol(hiddens))
  return(res)
}

sample_main_visible <- function(visibles, C) {
  #visibles is an NxV matrix of visible node values
  #C is the variance constant
  res <- sapply(colSums(visibles), function(mean) rnorm(1, mean*C, sqrt(C)))
  names(res) <- paste0("V", 1:ncol(visibles))
  return(res)
}

sample_interaction <- function(hiddens, visibles, C) {
  #hiddens is an NxH matrix of hidden node values
  #visibles is an NxV matrix of visible node values
  #C is the variance constant
  H <- ncol(hiddens)
  V <- ncol(visibles)
  
  stopifnot(nrow(hiddens) == nrow(visibles)) #something wrong with the data structure if these don't match
  
  res <- matrix(NA, nrow = H, ncol = V)
  for(i in 1:V) {
    for(j in 1:H) {
      res[i, j] <- rnorm(1, sum(hiddens[, j]*visibles[, i])*C, sqrt(C))
    }
  }
  
  rownames(res) <- paste0("V", 1:V)
  colnames(res) <- paste0("H", 1:H)
  
  return(res)
}

sample_main_mh <- function(theta, hiddens, visibles, tau, C) {
  #theta is an H + V + H*V vector of parameter values
  #hiddens is an NxH matrix of hidden node values
  #visibles is an NxV matrix of visible node values
  #tau is a variance parameter 
  V <- ncol(visibles)
  H <- ncol(hiddens)
  
  theta_star_main <- MASS::mvrnorm(n = 1, mu = theta[1:(V + H)], Sigma = diag(tau, nrow = V + H, ncol = H + V))
  theta_star <- c(theta_star_main, theta[(V + H + 1):(V + H + V*H)])
  u <- runif(1)
  
  N <- visibles %>% nrow()
  data <- cbind(visibles, hiddens)
  for(i in 1:V) {
    for(j in 1:H) {
      data <- cbind(data, visibles[i,]*hiddens[j,])
    }
  }
  possibles <- stats(H, V)
  W <- exp(possibles %*% theta)
  W_star <- exp(possibles %*% theta_star)
  
  A <- theta[1:(V + H)] %*% colSums(data[, 1:(V + H)]) - theta[1:(V + H)]^2 %*% (1/(2*rep(C, H + V)))
  A_star <- theta_star_main %*% colSums(data[, 1:(V + H)]) - theta_star_main^2 %*% (1/(2*rep(C, H + V)))
  
  if(u <= min(exp(-N*log(sum(W_star)) + N*log(sum(W)) + A_star - A), 1)) {
    return(theta_star_main)
  } else {
    return(theta[1:(V + H)])
  }
  
}

sample_interaction_mh <- function(theta, hiddens, visibles, tau, C_prime) {
  #theta is an H + V + H*V vector of parameter values
  #hiddens is an NxH matrix of hidden node values
  #visibles is an NxV matrix of visible node values
  #covar is an (H + V + H*V)x(H + V + H*V) covariance matrix
  V <- ncol(visibles)
  H <- ncol(hiddens)
  
  theta_star_interaction <- MASS::mvrnorm(n = 1, mu = theta[(V + H + 1):(V + H + V*H)], Sigma = diag(tau, nrow = V*H, ncol = H*V))
  theta_star <- c(theta[1:(V + H)], theta_star_interaction)
  u <- runif(1)
  
  N <- visibles %>% nrow()
  data <- cbind(visibles, hiddens)
  for(i in 1:V) {
    for(j in 1:H) {
      data <- cbind(data, visibles[i,]*hiddens[j,])
    }
  }
  possibles <- stats(H, V)
  W <- exp(possibles %*% theta)
  W_star <- exp(possibles %*% theta_star)
  
  A <- theta[(V + H + 1):(V + H + V*H)] %*% colSums(data[, (V + H + 1):(V + H + V*H)]) - theta[(V + H + 1):(V + H + V*H)]^2 %*% (1/(2*rep(C_prime, H*V)))
  A_star <- theta_star_interaction %*% colSums(data[, (V + H + 1):(V + H + V*H)]) - theta_star_interaction^2 %*% (1/(2*rep(C_prime, H*V)))
  
  if(u <= min(exp(-N*log(sum(W_star)) + N*log(sum(W)) + A_star - A), 1)) {
    return(theta_star_interaction)
  } else {
    return(theta[(V + H + 1):(V + H + V*H)])
  }
  
}

sample_single_theta_adaptive_mh <- function(theta, hiddens, visibles, C, trunc_const, shrink, h, s1, s2, index) {
  #theta is an H + V + H*V vector of parameter values
  #hiddens is an NxH matrix of hidden node values
  #visibles is an NxV matrix of visible node values
  #covar is an (H + V + H*V)x(H + V + H*V) covariance matrix
  require(dplyr)
  
  V <- ncol(visibles)
  H <- ncol(hiddens)
  
  parms <- list(visibles = visibles, hiddens = hiddens, H = H, V = V)
  
  h_vec <- rep(0, length(theta))
  h_vec[index] <- h
  
  N <- nrow(visibles)
  data <- cbind(visibles, hiddens)
  for(i in 1:V) {
    for(j in 1:H) {
      data <- cbind(data, visibles[i,]*hiddens[j,])
    }
  }
  possibles <- stats(H, V)
  W <- log_conditional(possibles, theta, data, index, C)

  H_h <- (W - 2*log_conditional(possibles, theta - h_vec, data, index, C) + log_conditional(possibles, theta - 2*h_vec, data, index, C))/h^2
  sigma2_tilde <- -1/H_h
  
  c_opt <- 2.4
  
  theta_star_single <- rnorm(n = 1, mean = theta[index], sd = sqrt(c_opt^2 * min(max(sigma2_tilde, s1), s2)))
  theta_star <- theta
  theta_star[index] <- theta_star_single
  
  const <- ifelse(shrink, C*trunc_const, C)
  if(theta_star_single^2 > const) return(list(theta = theta[index], var = NA)) #truncation
  
  u <- runif(1)
  
  val <- exp(log_conditional(possibles, theta_star, data, index, C) - W)
  val <- ifelse(!is.finite(val), sign(val) * .Machine$double.xmax, val)
  
  if(u <= min(val, 1)) { #accept/reject ratio
    return(list(theta = theta_star_single, var = min(max(sigma2_tilde, s1), s2)))
  } else {
    return(list(theta = theta[index], var = min(max(sigma2_tilde, s1), s2)))
  }
  
}

sample_hidden <- function(main_hidden, interaction, visibles) {
  #main_hidden is a 1xH matrix of hidden main effects param values
  #interaction is a HxV matrix of interaction effects param values
  #visibles is an NxV matrix of visible node values
  a <- exp(2*(matrix(main_hidden, nrow = nrow(visibles), ncol = length(main_hidden), byrow = TRUE) + crossprod(t(visibles), interaction)))
  a[!is.finite(a)] <- sign(a[!is.finite(a)]) * .Machine$double.xmax 
  res <- apply(a, 1:2, function(b) rbinom(1, 1, max(min(b/(1 + b), 1),0))*2 - 1)
  
  #a <- exp(2*(matrix(main_hidden, nrow = nrow(visibles), ncol = length(main_hidden), byrow = TRUE) + matrix(rowSums(interaction), nrow = nrow(visibles), ncol = nrow(interaction), byrow = TRUE)*visibles))
  #res <- apply(a, 1:2, function(b) rbinom(1, 1, b/(b+1))*2 - 1)
  colnames(res) <- paste0("H", 1:H)
  return(res)
}

sample_visible <- function(main_visible, interaction, hiddens) {
  #main_visible is a 1xH matrix of visible main effects param values
  #interaction is a HxV matrix of interaction effects param values
  #hiddens is an NxH matrix of hidden node values
  #a <- exp(2*(matrix(main_visible, nrow = nrow(hiddens), ncol = length(main_visible), byrow = TRUE) + matrix(colSums(interaction), nrow = nrow(hiddens), ncol = nrow(interaction), byrow = TRUE)*hiddens))
  #res <- apply(a, 1:2, function(b) rbinom(1, 1, b/(b+1))*2 - 1)
  
  a <- exp(2*(matrix(main_visible, nrow = nrow(hiddens), ncol = length(main_visible), byrow = TRUE) + tcrossprod(hiddens, interaction)))
  res <- apply(a, 1:2, function(b) rbinom(1, 1, b/(1 + b))*2 - 1)
  
  colnames(res) <- paste0("V", 1:H)
  return(res)
}

sample_images <- function(params, hidden0, visible0, mc.iter = 1e4) {
  #params is a list containing three named elements: main_hidden, main_visible, and interaction
  #hidden0 is an initialized NxH hidden node values
  #visible0 is an initialized NxV visible node values
  #mc.iter is the number of iterations to run
  
  # initialize data frame to save chains
  visible_save <- matrix(NA, ncol = length(visible0), nrow = mc.iter + 1) #flatten for storage
  hidden_save <- matrix(NA, ncol = length(hidden0), nrow = mc.iter + 1)
  
  # store initial values
  visible_save[1, ] <- visible0 #flattened so first row, then second row, etc.
  hidden_save[1, ] <- hidden0 #flattened so first row, then second row, etc.
  
  for(i in 2:(mc.iter + 1)) {
    hidden_save[i, ] <- sample_hidden(main_hidden = params$main_hidden, interaction = params$interaction, visibles = matrix(visible_save[i - 1, ], ncol = ncol(hidden0), byrow = TRUE))
    visible_save[i, ] <- sample_visible(main_visible = params$main_visible, interaction = params$interaction, hiddens = matrix(hidden_save[i, ], ncol = ncol(visible0), byrow = TRUE))
  }
  return(list(visibles = visible_save, hiddens = hidden_save))
}

sample_gibbs <- function(visibles, params0, C, C_prime, mc.iter = 1e4) {
  #visibles is an NxV matrix of images
  #params0 is a list containing three named elements of initial parameter values: main_hidden, main_visible, and interaction
  #C is the variance constant for the main effect
  #C_prime is the variance constant for the interaction
  #mc.iter is the number of iterations to run

  # store initial values
  params <- list(main_hidden = matrix(NA, nrow = mc.iter + 1, ncol = length(params0$main_hidden)),
                 main_visible = matrix(NA, nrow = mc.iter + 1, ncol = length(params0$main_visible)),
                 interaction = array(NA, dim=c(mc.iter + 1, length(params0$main_hidden), length(params0$main_visible))))
  hiddens <- array(NA, dim = c(mc.iter, nrow(visibles), length(params0$main_hidden)))
  #post_pred <- array(NA, dim = c(mc.iter, nrow(visibles), length(params0$main_visible)))
  distn <- array(NA, dim = c(mc.iter, 2^ncol(visibles), ncol(visibles) + 2))
  
  params$main_hidden[1, ] <- params0$main_hidden
  params$main_visible[1, ] <- params0$main_visible
  params$interaction[1, , ] <- params0$interaction
  
  for(i in 1:mc.iter) {
    cat(paste0(i, "\r"))
    #sample each type of variable in turn
    hiddens[i, , ] <- sample_hidden(main_hidden = params$main_hidden[i, ], interaction = params$interaction[i, , ], visibles = visibles)
    params$main_hidden[i + 1, ] <- sample_main_hidden(hiddens = hiddens[i, , ], C = C)
    params$main_visible[i + 1, ] <- sample_main_visible(visibles = visibles, C = C)
    params$interaction[i + 1, , ] <- sample_interaction(hiddens = hiddens[i, , ], visibles = visibles, C = C_prime)
    #post_pred[i, , ] <- sample_visible(main_visible = params$main_visible[i + 1, ], interaction = params$interaction[i + 1, , ], hiddens = hiddens[i, , ])
    distn[i, , ] <- visible_distn(list(main_hidden = params$main_hidden[i + 1, ],
                                       main_visible = params$main_visible[i + 1, ],
                                       interaction = params$interaction[i + 1, , ])) %>%
      data.matrix()
  }
  
  #remove initial value
  params$main_hidden <- params$main_hidden[-1, ]
  params$main_visible <- params$main_visible[-1, ]
  params$interaction <- params$interaction[-1, , ]
  
  return(list(params = params, hiddens = hiddens, distn = distn))
}

sample_mh_within_gibbs <- function(visibles, params0, C, C_prime, tau_main, tau_interaction, mc.iter = 1e4) {
  #visibles is an NxV matrix of images
  #params0 is a list containing three named elements of initial parameter values: main_hidden, main_visible, and interaction
  #C is the variance constant for the main effect
  #C_prime is the variance constant for the interaction
  #mc.iter is the number of iterations to run
  V <- length(params0$main_visible)
  H <- length(params0$main_hidden)
  
  # store initial values
  theta <- matrix(NA, nrow = mc.iter + 1, ncol = H + V + H*V)
  hiddens <- array(NA, dim = c(mc.iter, nrow(visibles), length(params0$main_hidden)))
  #post_pred <- array(NA, dim = c(mc.iter, nrow(visibles), length(params0$main_visible)))
  distn <- array(NA, dim = c(mc.iter, 2^ncol(visibles), ncol(visibles) + 2))
  
  theta[1, ] <- c(params0$main_visible, params0$main_hidden, t(params0$interaction))
  
  for(i in 1:mc.iter) {
    cat(paste0(i, "\r"))
    #sample each type of variable in turn
    hiddens[i, , ] <- sample_hidden(main_hidden = theta[i, (V+1):(V+H)], interaction = matrix(theta[i, (V+H+1):(V+H+H*V)], nrow = H, byrow = TRUE), visibles = visibles)
    theta[i + 1, ] <- c(sample_main_mh(theta = theta[i, ], hiddens = hiddens[i, ,], visibles = visibles, C = C, tau = tau_main), sample_interaction_mh(theta = theta[i, ], hiddens = hiddens[i, ,], visibles = visibles, C_prime = C_prime, tau = tau_interaction))
    #post_pred[i, , ] <- sample_visible(main_visible = params$main_visible[i + 1, ], interaction = params$interaction[i + 1, , ], hiddens = hiddens[i, , ])
    distn[i, , ] <- visible_distn(list(main_hidden = theta[i + 1, (V+1):(V+H)],
                                       main_visible = theta[i + 1, 1:V],
                                       interaction = matrix(theta[i + 1, (V+H+1):(V+H+H*V)], nrow = H, byrow = TRUE))) %>%
      data.matrix()
  }
  
  #remove initial value
  theta <- theta[-1, ]
  
  
  return(list(theta = theta, hiddens = hiddens, distn = distn))
}

sample_single_adaptive_mh_within_gibbs <- function(visibles, params0, C, C_prime, trunc_const, h, s1, s2, mc.iter = 1e4) {
  #visibles is an NxV matrix of images
  #params0 is a list containing three named elements of initial parameter values: main_hidden, main_visible, and interaction
  #C is the variance constant for the main effect
  #C_prime is the variance constant for the interaction
  #mc.iter is the number of iterations to run
  V <- length(params0$main_visible)
  H <- length(params0$main_hidden)
  
  # store initial values
  theta <- matrix(NA, nrow = mc.iter + 1, ncol = H + V + H*V)
  hiddens <- array(NA, dim = c(mc.iter, nrow(visibles), length(params0$main_hidden)))
  distn <- array(NA, dim = c(mc.iter, 2^ncol(visibles), ncol(visibles) + 2))
  var <- matrix(NA, nrow = mc.iter, ncol = H + V + H*V)
  
  theta[1, ] <- c(params0$main_visible, params0$main_hidden, t(params0$interaction))
  
  for(i in 1:mc.iter) {
    cat(paste0(i, "\r"))
    #sample each type of variable in turn
    hiddens[i, , ] <- sample_hidden(main_hidden = theta[i, (V+1):(V+H)], interaction = matrix(theta[i, (V+H+1):(V+H+H*V)], nrow = H, byrow = TRUE), visibles = visibles)
    for(j in 1:(H + V)) {
      samp <- sample_single_theta_adaptive_mh(theta = theta[i, ], hiddens = hiddens[i, ,], visibles = visibles, C = C, trunc_const = trunc_const, shrink = FALSE, h = h, s1 = s1, s2 = s2, index = j)
      theta[i + 1, j] <- samp$theta
      var[i, j] <- samp$var
    }
    for(j in (H + V + 1):(H + V + H*V)) {
      samp <- sample_single_theta_adaptive_mh(theta = theta[i, ], hiddens = hiddens[i, ,], visibles = visibles, C = C_prime, trunc_const = trunc_const, shrink = TRUE, h = h, s1 = s1, s2 = s2, index = j)
      theta[i + 1, j] <- samp$theta
      var[i, j] <- samp$var
    }
    distn[i, , ] <- visible_distn(list(main_hidden = theta[i + 1, (V+1):(V+H)],
                                       main_visible = theta[i + 1, 1:V],
                                       interaction = matrix(theta[i + 1, (V+H+1):(V+H+H*V)], nrow = H, byrow = TRUE))) %>%
      data.matrix()
  }
  
  #remove initial value
  theta <- theta[-1, ]
  
  
  return(list(theta = theta, hiddens = hiddens, distn = distn, var = var))
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

log_conditional <- function(possibles, theta, data, index, C) {
  N <- nrow(data)
  -N*log(sum(exp(possibles %*% theta))) + theta[index] * sum(data[, index]) - theta[index]^2 * (1/(2*C/9))
}

