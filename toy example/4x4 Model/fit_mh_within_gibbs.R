#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("sample_functs.R")

#data and params ------------------------
H <- 4
V <- 4
mc.iter <- 10000
set.seed(102285) #reproducible seed

load("written/sample_images.Rdata")
load("written/variance_params.Rdata")
load("written/sample.params.Rdata")

params <- list(main_hidden = rnorm(H, mean = 0, sd = sqrt(variance_params$C)),
               main_visible = rnorm(V, mean = 0, sd = sqrt(variance_params$C)),
               interaction = matrix(rnorm(H*V, mean = 0, sd = sqrt(variance_params$C_prime)), nrow = H))
tau1 <- 1
tau2 <- .25 #tuning parameter

models_good <- sample_mh_within_gibbs(visibles = flat_images_good$visibles, params0 = params, C = variance_params$C, C_prime = variance_params$C_prime, tau_main = tau1, tau_interaction = tau2, mc.iter = mc.iter)
models_bad <- sample_mh_within_gibbs(visibles = flat_images_good$visibles, params0 = params, C = variance_params$C, C_prime = variance_params$C, tau_main = tau1, tau_interaction = tau2, mc.iter = mc.iter)
models_degen <- sample_mh_within_gibbs(visibles = flat_images_degen$visibles, params0 = params, C = variance_params$C, C_prime = variance_params$C_prime, tau_main = tau1, tau_interaction = tau2, mc.iter = mc.iter)

#burnin + thinning
burnin_thinning <- function(models) {
  idx <- which(1:mc.iter %% 5 == 0 & 1:mc.iter > mc.iter/2)
  models$theta <- models$theta[idx, ]
  models$hiddens <- models$hiddens[idx, , ]
  #models$post_pred <- models$post_pred[idx, , ]
  models$distn <- models$distn[idx, , ]
  
  return(models)
}

models_good <- burnin_thinning(models_good)
models_bad <- burnin_thinning(models_bad)
models_degen <- burnin_thinning(models_degen)

save(models_good, models_bad, models_degen, file = "written/fitted_models_mh_distn.Rdata")
