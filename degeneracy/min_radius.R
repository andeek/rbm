source("functions.R")
library(dplyr)
library(reshape2)
library(ggplot2)

## constants
n <- 500 #num points to sample on the sphere
cutoff <- .01
step <- 1
backstep <- .01

test_cases <- expand.grid(data.frame(1:4)[,rep(1,2)]) %>% 
  rename(H = X1.4, V = X1.4.1) %>%
  filter(H <= V) %>% #remove those cases with less visibles than hiddens. Is this necessary?
  mutate(n_param = H*V + H + V) %>%
  mutate(max_facets = (2^(H+V))^(floor(n_param/2))) %>%
  filter(n_param <= 11) %>% #can't handle any higher dimensions currently
  rowwise() %>%
  mutate(min_radius(H, V, n=n, cutoff=cutoff, backstep=backstep, step=step))

min_radius <- function(H, V, type='negative', n, cutoff, backstep, step) {
  stat <- stats(H, V, type)
  
  h <- calc_hull(stat, options="Qt")
  triangles <- h %>%
    data.frame() %>%
    rowwise() %>%
    do(data.frame(stat[unlist(.),])) %>%
    ungroup() %>%
    transmute_(point=paste0("paste(",paste(colnames(stat), collapse=","),", sep=',')")) %>%
    mutate(triangle=rep(1:(n()/3), each = 3)) %>%
    mutate(point_num=paste0("point", rep(1:3, (n()/3)))) %>%
    dcast(triangle ~ point_num, value.var = "point")
  
  
  ## inits
  r0 <- 0
  theta0 <- sample_sphere(stat, n, r0)
  exp_vals0 <- expected_value(stats = stat, theta = theta0, normalized = TRUE)
  d0 <- exp_vals0 %>% t() %>%
    data.frame() %>%
    rowwise() %>% 
    transmute_(P=paste0("paste(",paste(rownames(exp_vals0), collapse=","),", sep=',')")) %>%
    merge(triangles) %>%    
    select(-triangle) %>%
    ungroup() %>%
    apply(1, function(x) distance_to_simplex(P=x["P"], simplex=as.list(x[c("point1", "point2", "point3")]))) %>%
    min()
 


  success <- FALSE
  while(!success){

    r1 <- r0 + step
    theta1 <- sample_sphere(stat, n, r1)
    exp_vals1 <- expected_value(stats = stat, theta = theta1, normalized = TRUE)
    d1 <- exp_vals1 %>% t() %>%
      data.frame() %>%
      rowwise() %>% 
      transmute_(P=paste0("paste(",paste(rownames(exp_vals0), collapse=","),", sep=',')")) %>%
      merge(triangles) %>%    
      select(-triangle) %>%
      ungroup() %>%
      apply(1, function(x) distance_to_simplex(P=x["P"], simplex=as.list(x[c("point1", "point2", "point3")]))) %>%
      min()

    if(d0 < cutoff & d1 >= cutoff) success <- TRUE
    
    if(d1 < cutoff) {
      r0 <- r1
      d0 <- d1
      step <- -backstep
    } else {
      r0 <- r1
      d0 <- d1
    }
    cat(d0)

  }
  return(r0)
}


