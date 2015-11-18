library(ggplot2)
library(dplyr)
library(tidyr)
source("functions.R")

model_data <- function(res) {
  
  model.data <- list()
  for(i in 1:nrow(res)) {
    tmp <- res$samp[[i]]
    H <- res[i,]$H
    V <- res[i,]$V
    
    tmp %>% 
      data.frame() %>% 
      rowwise() %>% 
      mutate_(ss_ratio = paste0("(", paste(paste0(names(.)[(H + V + 1):(H + V + H*V)], "^2"), collapse = " + "), ")/sum(", paste(paste0(names(.)[1:(H+V)],"^2"), collapse = " + "), ")")) %>%    
      ungroup() -> ratio
    
    inner_join(ratio, res$outside[[i]] %>% ungroup()) %>%
      mutate(H = H, V = V, n_param = H + V + H*V) -> model.data[[i]]
  }
  names(model.data) <- res$n_param
  return(model.data)
}

n <- 100
r_center <- seq(0.2, 3, by = 0.2)
r_exp <- .05

set.seed(10-22-1985) #best day ever!

model.dat <- list()
for(i in seq(0, 1, by = 0.1)) {
  load(paste0("apps/results_", i, ".RData")) 
  model_dat <- lapply(res, model_data)
  names(model_dat) <- r_center
  
  model.dat[[i*10 + 1]] <- model_dat
}
names(model.dat) <- seq(0, 1, by = 0.1)

expand.grid(power = as.numeric(names(model.dat)), multiplier = as.numeric(names(model.dat[[1]])), n_param = as.numeric(names(model.dat[[1]][[1]]))) %>%
  cbind(expand.grid(i = 1:length(model.dat), j = 1:length(model.dat[[1]]), k = 1:length(model.dat[[1]][[1]]))) %>%
  group_by(power, multiplier, n_param) %>%
  do(models = model.dat[[.$i]][[.$j]][[.$k]]) -> frame

frame %>%
  group_by(power, multiplier, n_param) %>%
  mutate(frac_degen = sum(models[[1]][, "near_hull"])/nrow(models[[1]]), 
         radius = (models[[1]][1, "H"] + models[[1]][1, "V"])^power,
         happy = frac_degen < .05) %>%
  group_by(n_param, happy) %>% 
  sample_n(20) -> models.samp


models.samp %>%
  group_by(n_param, happy) %>%
  do(models = do.call(rbind, .$models)) -> models.comb

no_int <- models.comb
for(i in 1:nrow(no_int)) {
  no_int$models[[i]][, grep("theta", names(no_int$models[[i]]))[!(grep("theta", names(no_int$models[[i]])) %in% grep("exp", names(no_int$models[[i]])))]] <- 0
  no_int$models[[i]]$happy <- no_int[i, ]$happy
}

no_int %>%
  group_by(n_param) %>%
  do(models = do.call(rbind, .$models)) %>%
  group_by(n_param) %>%
  do(compare = .$models[[1]] %>% 
       select(h1:exp_h1) %>%
       select(-starts_with("exp")) %>%
       data.matrix() %>%
       t()) %>%
  left_join(no_int %>%
              group_by(n_param) %>%
              do(stats = stats(.$models[[1]][1, c("H")], .$models[[1]][1, c("V")]))) %>%
  group_by(n_param) %>%
  do(no_int = expected_value(theta = .$compare[[1]], stats = .$stats[[1]]) %>% t()) %>%
  left_join(no_int %>% 
              group_by(n_param) %>% 
              do(int = do.call(rbind, .$models) %>% select(starts_with("exp")) %>% data.matrix())) %>%
  group_by(n_param) %>%
  do(diff = .$no_int[[1]] - .$int[[1]]) %>%
  group_by(n_param) %>%
  do(mean_abs_diff = .$diff[[1]] %>% apply(2, abs) %>% data.frame() %>% select(-contains("theta")) %>% apply(1, mean) %>% data.frame() %>% mutate(happy = c(rep(FALSE, 2000), rep(TRUE, 2000)))) %>%
  group_by(n_param) %>%
  do(mean_abs_diff = data.frame(.$mean_abs_diff[[1]], n_param = .$n_param)) -> abs

res <- data.frame()
for(i in 1:nrow(abs)) {
  res <- rbind(res, abs$mean_abs_diff[[i]])
}
names(res)[1] <- "value"

res %>%
  ggplot() +
  geom_histogram(aes(value)) +
  facet_grid(happy ~ n_param) +
  theme_bw(base_family = "serif")

res %>%
  group_by(n_param, happy) %>%
  summarise(mean_abs_diff = mean(value)) %>%
  spread(n_param, mean_abs_diff)
