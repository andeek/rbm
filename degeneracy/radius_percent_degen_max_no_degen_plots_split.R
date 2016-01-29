library(ggplot2)
library(dplyr)
library(tidyr)

plot_data <- function(res) {
  
  plot.data <- data.frame()
  for(i in 1:nrow(res)) {
    tmp <- res$samp[[i]]
    H <- res[i,]$H
    V <- res[i,]$V
    
    tmp %>% 
      data.frame() %>% 
      rowwise() %>%
      mutate_(ss_ratio = paste0("(", paste(paste0(names(.)[(H + V + 1):(H + V + H*V)], "^2"), collapse = " + "), ")/sum(", paste(paste0(names(.)[1:(H+V)],"^2"), collapse = " + "), ")")) %>%    
      ungroup() -> ratio
    
    tmp %>%
      data.frame() %>%
      select(-starts_with("exp")) %>%
      rowwise() %>%
      mutate_(true_r = paste0("sqrt(", paste(paste0(names(.), "^2"), collapse = " + "), ")")) -> radius
      
    inner_join(radius, ratio) %>%
    inner_join(res$outside[[i]] %>% ungroup()) %>%
      select(ss_ratio, near_hull, true_r) %>%
      mutate(H = H, V = V, n_param = H + V + H*V) %>%
      rbind(plot.data) -> plot.data
  }
  return(plot.data)
}

n <- 100
r_center <- seq(0.2, 3, by = 0.2)
r_exp <- .05


plot.dat <- list()
for(i in seq(0, 1, by = 0.1)) {
  load(paste0("apps/results_split_", i, ".RData")) 
  plot_dat <- lapply(res, plot_data)
  names(plot_dat) <- gsub("[.]", "_", r_center)
  
  do.call(rbind, plot_dat) %>%
    ungroup() %>%
    mutate(name = rownames(do.call(rbind, plot_dat))) %>%
    separate(name, into = c("multiplier", "num"), "[.]") %>% 
    mutate(multiplier = as.numeric(gsub("_", ".", multiplier))) %>%
    mutate(n_param_f = factor(n_param)) -> plot_dat
  
  plot.dat[[i*10 + 1]] <- plot_dat
}

names(plot.dat) <- gsub("[.]", "_", seq(0, 1, by = 0.1))
do.call(rbind, plot.dat) %>%
  mutate(name2 = rownames(do.call(rbind, plot.dat))) %>% 
  separate(name2, into = c("power", "num2"), "[.]") %>% 
  mutate(power = as.numeric(gsub("_", ".", power))) -> plot.dat


plot.dat %>%
  mutate(split_r = multiplier*(H + V)^power) %>%
  group_by(true_r, split_r, n_param) %>%
  summarise(frac_degen = sum(near_hull)/n()) %>%
  filter(frac_degen < .05) %>% 
  group_by(n_param) %>%
  summarise(best_radius = max(true_r),
            split_radius = max(split_r)) -> best_rad_nodegen

best_rad_nodegen %>%
  gather(radii, value, best_radius, split_radius) %>%
  ggplot() +
  geom_point(aes(n_param, value, colour = radii)) +
  geom_smooth(aes(n_param, value, colour = radii), method = "lm") +
  theme_bw(base_family = "serif") +
  theme(legend.position = "bottom")

plot.dat %>%
  mutate(split_r = multiplier*(H + V)^power) %>%
  group_by(true_r, split_r, n_param_f, multiplier, power) %>%
  summarise(frac_degen = sum(near_hull)/n()) %>%
  ggplot() +
  geom_jitter(aes(true_r, frac_degen, colour = factor(power)), alpha = .5) +
  facet_grid(power~n_param_f) +
  theme_bw(base_family = "serif") +
  theme(legend.position = "bottom")

plot.dat %>%
  mutate(split_r = multiplier*(H + V)^power) %>%
  group_by(true_r, split_r, n_param, multiplier, power) %>%
  summarise(frac_degen = sum(near_hull)/n()) %>%
  mutate(happy = frac_degen < .05) %>%
  group_by(n_param, happy) %>%
  summarise(count = n())