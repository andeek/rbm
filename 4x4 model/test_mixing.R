library(dplyr)
library(tidyr)

load("written/sample_images.Rdata")

#data and params ------------------------
H <- 4
V <- 4

load("written/fitted_models_trunc_marginal_full.Rdata")
marginal <- models_good

load("written/fitted_models_trunc_full.Rdata")
trunc <- models_good

rm(models_good)
rm(models_bad)

source("functs_sample.R")
load("written/params_theta.Rdata")

params_good <- list(main_hidden = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                    main_visible = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                    interaction = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))
distn_good <- visible_distn(params = params_good)


reshape_sample_distn <- function(model) {
  sample_distn <- model$distn
  dim(sample_distn) <- c(dim(sample_distn)[1]*dim(sample_distn)[2], dim(sample_distn)[3])
  sample_distn %>% data.frame() -> sample_distn
  names(sample_distn) <- names(distn_good)
  
  sample_distn %>%
    group_by(image_id) %>%
    mutate(iter = 1:n()) -> sample_distn
  
  return(sample_distn)
}

marginal_sample <- reshape_sample_distn(marginal)
trunc_sample <- reshape_sample_distn(trunc)


## acf -------------------------------
marginal_sample %>% 
  group_by(v1, v2, v3, v4, image_id) %>% 
  do(data.frame(acf = acf(.$prob, plot=FALSE)$acf,
     lag = acf(.$prob, plot=FALSE)$lag)) -> marginal_acfs

trunc_sample %>% 
  group_by(v1, v2, v3, v4, image_id) %>% 
  do(data.frame(acf = acf(.$prob, plot=FALSE)$acf,
                lag = acf(.$prob, plot=FALSE)$lag)) -> trunc_acfs

marginal_acfs %>% ungroup() %>% select(-image_id) %>% rename(marginal = acf) %>%
  left_join(trunc_acfs %>% rename(truncated = acf)) %>%
  ggplot() +
  geom_segment(aes(x=lag, xend=lag, y=0, yend=truncated), colour = "black") +
  geom_segment(aes(x=lag, xend=lag, y=0, yend=marginal), colour = "red") +
  # geom_hline(aes(yintercept = qnorm((1 + .95)/2)/sqrt(2000)), colour="red", lty=2) +
  facet_wrap(~image_id) +
  xlab("Lag") +
  ylab("ACF")

## block bootstrap effective sample size ----------------
###
### MBB Statistic
###
mbb_mean <- function(data_star, data, b) {
  m <- length(data_star)
  n <- length(data)
  N <- n - b + 1
  
  blockmeans <- rep(0, N)
  for(i in 1:N) {
    blockmeans[i] <- mean(data[i:(b + i - 1)])
  }
  
  sqrt(m)*(mean(data_star) - mean(blockmeans))
}

###
### MBB Sample
###
mbb_sample <- function(data, b) {
  n <- length(data)
  k <- floor(n/b)
  
  start <- sample.int(n - b + 1, k, replace = TRUE)
  data_star <- NULL
  for(i in 1:k) {
    data_star <- c(data_star, data[start[i]:(start[i] + b - 1)])
  }
  
  return(data_star)
}

###
### Block Bootstrap function
###
bb <- function(x, b, B, boot_stat, boot_sample) {
  #ii. For blocksize b, generate B MBB versions of statistic
  bb.stat <- array(NA, B)
  for(i in seq_len(B)) {
    bb.stat[i] <- boot_stat(boot_sample(x, b), x, b)
  }
  return(bb.stat)
}

marginal_sample %>% 
  group_by(v1, v2, v3, v4, image_id) %>%
  do(mbb = bb(.$prob, 2000^(1/3), 2000, mbb_mean, mbb_sample)) -> marginal_mbb

trunc_sample %>% 
  group_by(v1, v2, v3, v4, image_id) %>%
  do(mbb = bb(.$prob, 5000^(1/3), 2000, mbb_mean, mbb_sample)) -> trunc_mbb

epsilon <- 1/1000

marginal_mbb %>%
  group_by(v1, v2, v3, v4, image_id) %>%
  do(data.frame(C = var(.$mbb[[1]]))) %>%
  left_join(marginal_sample %>% group_by(v1, v2, v3, v4) %>% summarise(sigma2 = var(prob))) %>%
  mutate(T_star = 1/epsilon * C/sigma2) %>%
  ungroup() %>%
  summarise(T_star = max(T_star)) %>% as.numeric() %>% ceiling() -> marginal_T

trunc_mbb %>%
  group_by(v1, v2, v3, v4, image_id) %>%
  do(data.frame(C = var(.$mbb[[1]]))) %>%
  left_join(trunc_sample %>% group_by(v1, v2, v3, v4) %>% summarise(sigma2 = var(prob))) %>%
  mutate(T_star = 1/epsilon * C/sigma2) %>%
  ungroup() %>%
  summarise(T_star = max(T_star)) %>% as.numeric() %>% ceiling() -> trunc_T

## timing ---------------------------
H <- 4
V <- 4
mc.iter <- 100
params <- list(main_hidden = rep(0, H),
               main_visible = rep(0, V),
               interaction = rep(0, H*V))

h <- .01
s1 <- 1e-4
s2 <- 1e4
trunc_const <- 1

#rule of thumb for variance params
C <- 1/(H + V)
C_prime = 1/(H*V)

time <- Sys.time()
sample_single_adaptive_mh_within_gibbs(visibles = flat_images_good$visibles, params0 = params, C = C, C_prime = C_prime, trunc_const = trunc_const, h = h, s1 = s1, s2 = s2, mc.iter = mc.iter)
trunc_time <- difftime(Sys.time(), time, units = "secs")

time <- Sys.time()
sample_single_adaptive_mh_within_gibbs(visibles = flat_images_good$visibles, params0 = params, C = C, C_prime =  C_prime, trunc_const = trunc_const, h = h, s1 = s1, s2 = s2, mc.iter = mc.iter, conditional_function = log_conditional_marginal)
marginal_time <- difftime(Sys.time(), time, units = "secs")


trunc_time/mc.iter*trunc_T
marginal_time/mc.iter*marginal_T

## note: even still, marginal is way slower.

