#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("sample_functs.R")

#data and params ------------------------
H <- 4
V <- 4
set.seed(102285) #reproducible seed
mc.iter <- 200

load("written/sample_images.Rdata")
load("written/variance_params.Rdata")
load("written/sample.params.Rdata")
load("written/tune_models_mh.RData")

params <- list(main_hidden = rnorm(H, mean = 0, sd = sqrt(variance_params$C)),
               main_visible = rnorm(V, mean = 0, sd = sqrt(variance_params$C)),
               interaction = matrix(rnorm(H*V, mean = 0, sd = sqrt(variance_params$C_prime)), nrow = H))

# data.frame(expand.grid(tau_main = seq(.001, .02, by = .002), tau_interaction = seq(.001, .02, by = .002))) %>%
#   group_by(tau_main, tau_interaction) %>%
#   do(model = sample_mh_within_gibbs(visibles = flat_images_good$visibles, params0 = params, C = variance_params$C, C_prime = variance_params$C_prime, tau_main = .$tau_main, tau_interaction = .$tau_interaction, mc.iter = mc.iter)) -> models
# 
# save(models, file = "written/tune_models_mh.RData")

apply(models, 1, function(x) {
  data.frame(tau_main = x$tau_main,
            tau_interaction = x$tau_interaction,
            acceptance_rate_main = sum(diff(x$model$theta[, 1]) != 0)/mc.iter,
            acceptance_rate_interaction = sum(diff(x$model$theta[, 9]) != 0)/mc.iter,
            acceptance_rate_combined = sum(diff(x$model$theta[, 1]) != 0 | diff(x$model$theta[, 9]) != 0)/mc.iter)
}) -> models_tune

models_tune <- do.call(rbind, models_tune)
  

models_tune %>%
  group_by(tau_main) %>%
  summarise(acceptance_rate_main = mean(acceptance_rate_main)) %>%
  ggplot() +
  geom_line(aes(tau_main, acceptance_rate_main)) +
  geom_abline(aes(intercept = .2, slope = 0), colour = "red") +
  ylim(c(0,1))

models_tune %>%
  group_by(tau_interaction) %>%
  summarise(acceptance_rate_interaction = mean(acceptance_rate_interaction)) %>%
  ggplot() +
  geom_line(aes(tau_interaction, acceptance_rate_interaction)) +
  geom_abline(aes(intercept = .2, slope = 0), colour = "red") +
  ylim(c(0,1))

models_tune %>%
  ggplot() +
  #geom_tile(aes(x = tau_main, y = tau_interaction, fill = acceptance_rate_combined, alpha = acceptance_rate_combined)) +
  stat_contour(aes(x = tau_main, y = tau_interaction, z = acceptance_rate_combined, fill = ..level.., alpha = ..level..), binwidth = .005, geom = "polygon") +
  stat_contour(aes(x = tau_main, y = tau_interaction, z = acceptance_rate_combined), breaks = .2) +
  scale_fill_gradient(low = "yellow", high = "red")