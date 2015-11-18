#libraries -----------------------
library(dplyr)
library(tidyr)
library(ggplot2)

#reshape data ---------------------
load("../expected value/apps/results_split.RData")
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

#choose C, C_prime --------------
plot_dat %>%
  filter(H == 4 & V == 4) %>%
  group_by(H, V, C, epsilon) %>%
  summarise(percent_degen = sum(near_hull)/n()) %>%
  mutate(C_prime = C*((H + V)/(H*V))^(1+epsilon)) %>%
  mutate(percent_degen = round(percent_degen, 2)) %>%
  group_by(percent_degen) %>%
  mutate(C_C_prime = C + C_prime) %>%
  top_n(1, C_C_prime) %>%
  arrange(percent_degen, C_prime, C) %>%
  ungroup() %>%
  top_n(1, C_C_prime) %>%
  select(percent_degen, C, C_prime) %>%
  data.frame -> variance_params

save(variance_params, file = "written/variance_params.Rdata")
