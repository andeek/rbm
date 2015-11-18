#libraries -----------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("functions.R")

## main function --------------------------------------------
find_prop <- function(g_theta, stat, r_exp) {
  g_theta %>%
    data.frame() %>%
    group_by_(.dots = c(colnames(g_theta))) %>%
    do(samp = sample_points_outside(as.numeric(.), stat, r_exp)$in_hull) %>%
    group_by_(.dots = c(colnames(g_theta))) %>%
    do(as.data.frame(sum(.$samp[[1]]) > 0)) -> is_outside # count of not in hull greater than zero => there exists a point outside hull for this model
  
  is_outside %>% 
    rename_(near_hull = as.name(names(is_outside)[ncol(is_outside)])) %>%
    mutate(r = r_exp) -> inside_outside
  
  return(inside_outside)
}

## additional useful functions -------------------------------------
sample_points_outside <- function(g_theta, stats, r, n = 100) {
  #g_theta is the expected value mapping of a point theta in R^H+V+H*V
  x0 <- data.frame(samp = 1:n) %>%
    mutate(vals = samp_exp(length(g_theta), g_theta, r)) %>%
    separate(col = vals, sep = ",", into = colnames(stats)) %>%
    select(-samp)
  
  x0$in_hull <- !apply(x0, 1, in_hull, stats)
  
  return(x0)
}

samp_exp <- function(space, center, r) {
  vec <- rnorm(space)
  paste0(r*vec/sqrt(sum(vec^2)) + center, collapse=",")
}

in_hull <- function(point, hull_points) {
  require(lpSolveAPI)
  P <- data.matrix(point)
  A <- t(data.matrix(hull_points))
  
  lp_obj <- make.lp(nrow = 0, ncol = ncol(A))
  sapply(1:nrow(A), function(x) add.constraint(lp_obj, A[x,], type = "=", rhs = P[x]))
  add.constraint(lp_obj, rep(1, ncol(A)), type = "=", rhs = 1)
  sapply(1:ncol(A), function(x) add.constraint(lp_obj, 1, type = ">=", rhs = 0, indices = x))
  
  return(solve(lp_obj) == 0)
}

#sample split radii ---------------------------------------------
n <- 100
C <- seq(0.2, 3, by = 0.2)
r_exp <- .05
epsilon <- seq(0, .5, by = .05)

expand.grid(H = 1:4, V = 1:4, C = C, epsilon = epsilon) %>%
  mutate(N = H + V) %>%
  mutate(n = n, r1 = sqrt(C), r2 = sqrt(C*(N/(H*V))^(1+epsilon))) %>%
  group_by(H, V, C, epsilon, r1, r2, N) %>%
  do(samp = cbind(matrix(rnorm(n*.$N, mean = 0, sd = .$r1/3), nrow = n), matrix(rnorm(n*.$H*.$V, mean = 0, sd = .$r2/3), nrow = n))) -> split_sample

expand.grid(H = 1:4, V = 1:4, C = C, epsilon = epsilon) %>%
  mutate(N = H + V) %>%
  mutate(n = n, r1 = C, r2 = C*(N/(H*V))^(1+epsilon)) %>%
  group_by(H, V, C, epsilon, r1, r2, N) %>%
  do(stat = stats(.$H, .$V, "negative")) -> split_sample_stat

split_sample <- inner_join(split_sample, split_sample_stat)

split_sample %>%
  group_by(H, V, N, C, epsilon, r1, r2) %>%
  do(g_theta = t(rbind(t(.$samp[[1]]), expected_value(theta = t(.$samp[[1]]), stats = data.matrix(.$stat[[1]]))))) -> split

split <- inner_join(split, split_sample)

split %>%
  group_by(H, V, N, C, epsilon, r1, r2) %>%
  do(outside = find_prop(.$g_theta[[1]] %>% data.frame() %>% select(starts_with("exp")), .$stat[[1]], r_exp)) -> tmp

res <- inner_join(tmp, split) 
save(res, file = "apps/results_split.RData")


# plots ------------------------
load("apps/results_split.RData")

plot_data <- function(res) {
  
  plot.data <- data.frame()
  for(i in 1:nrow(res)) {
    tmp <- res$g_theta[[i]] %>% data.frame()
    H <- res[i,]$H
    V <- res[i,]$V
    N <- res[i,]$N
    r1 <- res[i,]$r1
    r2 <- res[i,]$r2
    C <- res[i,]$C
    epsilon <- res[i,]$epsilon
    
    
    tmp %>% 
      rowwise() %>% 
      mutate_(ss_ratio = paste0("(", paste(paste0(names(.)[(H + V + 1):(H + V + H*V)], "^2"), collapse = " + "), ")/sum(", paste(paste0(names(.)[1:(H+V)],"^2"), collapse = " + "), ")")) %>%    
      ungroup() -> ratio
    
    
    inner_join(ratio, res$outside[[i]] %>% ungroup()) %>%
      select(ss_ratio, near_hull) %>%
      mutate(H = H, V = V, n_param = H + V + H*V, N = H + V, N = N, r1 = r1, r2 = r2, C = C, epsilon = epsilon) %>%
      rbind(plot.data) -> plot.data
  }
  return(plot.data)
}

plot_dat <- res %>% plot_data

plot_dat %>%
  group_by(H, V, n_param, N, r1, r2, C, epsilon) %>%
  summarise(percent_degen = sum(near_hull)/n()) %>%
  ggplot() +
  geom_point(aes(C, percent_degen, colour = factor(n_param))) +
  geom_line(aes(C, percent_degen, colour = factor(n_param))) +
  facet_wrap(~epsilon) +
  ylim(c(0, 1))

plot_dat %>%
  group_by(H, V, n_param, N, r1, r2, C, epsilon) %>%
  summarise(percent_degen = sum(near_hull)/n()) %>%
  ggplot() +
  geom_point(aes(n_param, percent_degen)) +
  facet_grid(epsilon~C) +
  ylim(c(0, 1))
  
plot_dat %>%
  filter(epsilon == 0) %>%
  ggplot() +
  geom_boxplot(aes(factor(n_param), ss_ratio, colour = near_hull)) + 
  facet_wrap(~C, scales = "free_y")

# models ---------------------------------------------
plot_dat %>%
  mutate(C = r1^2, C_prime = r2^2) %>%
  select(near_hull, H, V, C, C_prime) %>%
  group_by(H, V, C, C_prime) %>%
  summarise(percent_degen = sum(near_hull)/n(),
            weights = n()) -> degen.dat

degen.m0 <- glm(formula = percent_degen ~ 1, family = "binomial", data = degen.dat, weights = weights)
degen.m1 <- glm(formula = percent_degen ~ H + V + C + C_prime, family = "binomial", data = degen.dat, weights = weights)
degen.m2 <- glm(formula = percent_degen ~ (H + V + C + C_prime)^2, family = "binomial", data = degen.dat, weights = weights)

#goodness of fit
-2*(logLik(degen.m0) - logLik(degen.m1)) >= qchisq(.05, df = degen.m0$df.residual - degen.m1$df.residual)
-2*(logLik(degen.m1) - logLik(degen.m2)) >= qchisq(.05, df = degen.m1$df.residual - degen.m2$df.residual) #use full model

find_C_C_prime <- function(p, H, V, Cs, C_primes, model) {
  LHS <- log(p/(1-p))
  grid <- expand.grid(C = Cs, C_prime = C_primes)
  new.x <- data.frame(H = rep(H, nrow(grid)), V = rep(V, nrow(grid)), grid)
  RHS <- data.frame(grid, pred = model.matrix(update(model$formula, NULL ~ .), new.x) %*% coef(model))
  return(list(LHS = LHS, RHS = RHS))
}


p <- .05
data.frame(expand.grid(H = seq(1, 11, by = 2), V = seq(1, 11, by = 2))) %>%
  group_by(H, V) %>%
  do(finding = find_C_C_prime(p, .$H, .$V, Cs = seq(0, 3, by = .1), C_primes =  seq(0, 3, by = .1), degen.m2)) %>%
  do(plotting = data.frame(H = .$H,
                           V = .$V, C = .$finding$RHS$C, 
                           C_prime = .$finding$RHS$C_prime, 
                           pred = .$finding$RHS$pred, 
                           pred_trans = exp(.$finding$RHS$pred)/(1 + exp(.$finding$RHS$pred)),
                           p = p)) -> explore


do.call(rbind, explore$plotting) %>% 
  ggplot(aes(x = C, y = C_prime, z = pred_trans)) +
  stat_contour(aes(colour = ..level..), size = 0.5) +
  stat_contour(breaks = c(p), size = 2, colour = "black") +
  facet_grid(H ~ V) +
  theme_bw()

data.frame(expand.grid(H = 4, V = 4)) %>%
  group_by(H, V) %>%
  do(finding = find_C_C_prime(p, .$H, .$V, Cs = seq(0, 3, by = .1), C_primes =  seq(0, 3, by = .1), degen.m2)) %>%
  do(plotting = data.frame(H = .$H,
                           V = .$V, C = .$finding$RHS$C, 
                           C_prime = .$finding$RHS$C_prime, 
                           pred = .$finding$RHS$pred, 
                           pred_trans = exp(.$finding$RHS$pred)/(1 + exp(.$finding$RHS$pred)),
                           p = p)) -> explore_44

do.call(rbind, explore_44$plotting) %>% 
  ggplot(aes(x = C, y = C_prime, z = pred_trans)) +
  stat_contour(aes(colour = ..level..), binwidth = .01, size = 0.5) +
  stat_contour(breaks = c(p), size = 2, colour = "black") +
  facet_grid(H ~ V) +
  theme_bw()

explore_44[[1, "plotting"]] %>% mutate(close_p = abs(p - pred_trans)) %>% arrange(close_p) %>% head
