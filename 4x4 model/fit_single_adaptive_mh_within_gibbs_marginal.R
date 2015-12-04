#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("sample_functs.R")

#data and params ------------------------
H <- 4
V <- 4
mc.iter <- 2000
set.seed(102285) #reproducible seed

load("written/sample_images.Rdata")
load("written/variance_params.Rdata")
load("written/sample.params.Rdata")

params <- list(main_hidden = rep(0, H),
               main_visible = rep(0, V),
               interaction = rep(0, H*V))

h <- .01
s1 <- 1e-4
s2 <- 1e4
trunc_const <- 1

## variance params to match Jing's prior results
#C <- 1/8 
#C_prime <- (.4)^2/16

models_good <- sample_single_adaptive_mh_within_gibbs(visibles = flat_images_good$visibles, params0 = params, C = variance_params$C, C_prime =  variance_params$C_prime, trunc_const = trunc_const, h = h, s1 = s1, s2 = s2, mc.iter = mc.iter, conditional_function = log_conditional_marginal)
models_bad <- sample_single_adaptive_mh_within_gibbs(visibles = flat_images_good$visibles, params0 = params, C = variance_params$C, C_prime = variance_params$C, trunc_const = trunc_const, h = h, s1 = s1, s2 = s2, mc.iter = mc.iter, conditional_function = log_conditional_marginal)
#models_degen <- sample_single_adaptive_mh_within_gibbs(visibles = flat_images_degen$visibles, params0 = params, C = variance_params$C, C_prime = variance_params$C_prime, trunc_const = trunc_const, h = h, s1 = s1, s2 = s2, mc.iter = mc.iter, conditional_function = log_conditional_marginal)

save(models_good, models_bad, file = "written/fitted_models_adaptive_mh_full_trunc_distn_marginal.Rdata")

#burnin + thinning
burnin_thinning <- function(models) {
  idx <- which(1:mc.iter %% 5 == 0 & 1:mc.iter > mc.iter/2)
  models$theta <- models$theta[idx, ]
  models$hiddens <- models$hiddens[idx, , ]
  #models$post_pred <- models$post_pred[idx, , ]
  models$distn <- models$distn[idx, , ]
  models$var <- models$var[idx, ]
  
  return(models)
}

models_good <- burnin_thinning(models_good)
models_bad <- burnin_thinning(models_bad)
#models_degen <- burnin_thinning(models_degen)

save(models_good, models_bad, file = "written/fitted_models_adaptive_mh_trunc_distn_marginal.Rdata")