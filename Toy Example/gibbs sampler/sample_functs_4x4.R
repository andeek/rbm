sample_main_hidden <- function(hiddens, C) {
  #hiddens is an NxH matrix of hidden node values
  #C is the variance constant
  res <- sapply(colSums(hiddens), function(mean) rnorm(1, mean/(2*C), 1/(2*C^2)))
  names(res) <- paste0("H", 1:ncol(hiddens))
  return(res)
}

sample_main_visible <- function(visibles, C) {
  #visibles is an NxV matrix of visible node values
  #C is the variance constant
  res <- sapply(colSums(visibles), function(mean) rnorm(1, mean/(2*C), 1/(2*C^2)))
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
  for(i in 1:H) {
    for(j in 1:V) {
      res[i, j] <- rnorm(1, sum(hiddens[, i]*hiddens[, j])/(2*C), 1/(2*C^2))
    }
  }
  
  colnames(res) <- paste0("V", 1:V)
  rownames(res) <- paste0("H", 1:H)
  
  return(res)
}

sample_hidden <- function(main_hidden, interaction, visibles) {
  #main_hidden is a 1xH matrix of hidden main effects param values
  #interaction is a HxV matrix of interaction effects param values
  #visibles is an NxV matrix of visible node values
  
  a <- exp(2*(sum(main_hidden) + rowSums(t(colSums(interaction) * t(visibles)))))
  res <- t(sapply(a, function(b) { rbinom(length(main_hidden), 1, b/(b + 1))*2 - 1 })) #make 0,1 -> -1,1
  colnames(res) <- paste0("H", 1:H)
  return(res)
}

sample_visible <- function(main_visible, interaction, hiddens) {
  #main_visible is a 1xH matrix of visible main effects param values
  #interaction is a HxV matrix of interaction effects param values
  #hiddens is an NxH matrix of hidden node values
  
  a <- exp(2*(sum(main_visible) + rowSums(t(colSums(interaction) * t(hiddens)))))
  res <- t(sapply(a, function(b) { rbinom(length(main_visible), 1, b/(b + 1))*2 - 1 })) #make 0,1 -> -1,1
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

  
  # store initial values
  visible_save[1, ] <- matrix(t(visible0), nrow = 1) #flattened so first row, then second row, etc.
  hiddens <- hidden0
  
  for(i in 2:(mc.iter + 1)) {
    hiddens <- sample_hidden(main_hidden = params$main_hidden, interaction = params$interaction, visibles = matrix(visible_save[i - 1, ], ncol = ncol(visible0), byrow = TRUE))
    visible_save[i, ] <- matrix(t(sample_visible(main_visible = params$main_visible, interaction = params$interaction, hiddens = hiddens)), nrow = 1) #flattened so first row, then second row, etc.
  }
  return(visible_save)
}

sample_gibbs <- function(visibles, hiddens0, params0, C, C_prime, mc.iter = 1e4) {
  #visibles is an NxV matrix of images
  #hidden0 is an initialized NxH hidden node values
  #params0 is a list containing three named elements of initial parameter values: main_hidden, main_visible, and interaction
  #C is the variance constant for the main effect
  #C_prime is the variance constant for the interaction
  #mc.iter is the number of iterations to run
  
  
  #TODO finish this
  #relationship of C, C_prime to variances in sampler?
  #which Cs to use? -- see old analysis?
  
  
  # initialize data frame to save chains
  visible_save <- matrix(NA, ncol = length(visible0), nrow = mc.iter + 1) #flatten for storage
  
  
  # store initial values
  visible_save[1, ] <- matrix(t(visible0), nrow = 1) #flattened so first row, then second row, etc.
  hiddens <- hidden0
  
  for(i in 2:(mc.iter + 1)) {
    hiddens <- sample_hidden(main_hidden = params$main_hidden, interaction = params$interaction, visibles = matrix(visible_save[i - 1, ], ncol = ncol(visible0), byrow = TRUE))
    visible_save[i, ] <- matrix(t(sample_visible(main_visible = params$main_visible, interaction = params$interaction, hiddens = hiddens)), nrow = 1) #flattened so first row, then second row, etc.
  }
  return(visible_save)
}


