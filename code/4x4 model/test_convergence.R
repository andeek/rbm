library(dplyr)
library(ggplot2)
library(tidyr)

load("written/fitted_models_trunc_marginal_full.Rdata")
load("written/params_theta.Rdata")
theta_good <- sample.params %>% 
  ungroup() %>% 
  filter(!near_hull) %>% 
  select(starts_with("v"), starts_with("h"), starts_with("theta"), -H, -V) %>% 
  data.matrix()


if("params" %in% names(models_good)) {
  theta_est <- data.frame(models_good$params$main_visible, models_good$params$main_hidden, models_good$params$interaction)
} else if("theta" %in% names(models_good)) {
  theta_est <- models_good$theta %>% data.frame()
}
names(theta_est) <- colnames(theta_good)


theta_est %>%
  mutate(iter = 1:n()) %>% 
  gather(variable, value, -iter) %>% 
  ggplot() + 
  geom_abline(aes(intercept = true_value, slope = 0), data = theta_good %>% data.frame %>% gather(variable, true_value), colour = "red") +
  geom_line(aes(iter, value)) + 
  facet_wrap(~variable, scales = "free")

theta_est %>%
  gather(variable, value) %>% 
  ggplot() + 
  geom_histogram(aes(value)) + 
  facet_wrap(~variable, scales = "free")

var <- models_good$var %>% data.frame
names(var) <- colnames(theta_good)

var %>%
  mutate(iter = 1:n()) %>% 
  gather(variable, value, -iter) %>% 
  ggplot() + 
  #geom_abline(aes(intercept = true_value, slope = 0), data = theta_good %>% data.frame %>% gather(variable, true_value), colour = "red") +
  geom_line(aes(iter, value)) + 
  facet_wrap(~variable, scales = "free")


theta_est %>%
  mutate(iter = 1:n()) %>% 
  mutate_(ss_main = paste(paste0(names(theta_est)[1:8], "^2"), collapse=" + "),
          ss_inter = paste(paste0(names(theta_est)[9:24], "^2"), collapse=" + ")) %>%
  gather(type, ss, ss_main, ss_inter) %>%
  ggplot() +
  geom_line(aes(iter, ss, colour = type))


