#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("sample_functs.R")

#data and params ------------------------
H <- 4
V <- 4
mc.iter <- 5000
set.seed(102285) #reproducible seed

load("written/sample_images.Rdata")
load("written/params_theta.Rdata")
load("written/params_variance.Rdata")

params <- list(main_hidden = rnorm(H, mean = 0, sd = sqrt(variance_params$C)),
               main_visible = rnorm(V, mean = 0, sd = sqrt(variance_params$C)),
               interaction = matrix(rnorm(H*V, mean = 0, sd = sqrt(variance_params$C_prime)), nrow = H))

N <- nrow(flat_images_good$visibles)
consts <- 6:8 #tune constant 

for(const in consts) {
  models_good <- sample_gibbs(visibles = flat_images_good$visibles, params0 = params, C = const*variance_params$C/N, C_prime = const*variance_params$C_prime/N, mc.iter = mc.iter)
  models_bad <- sample_gibbs(visibles = flat_images_good$visibles, params0 = params, C = const*variance_params$C/N, C_prime = const*variance_params$C/N, mc.iter = mc.iter)
  #models_degen <- sample_gibbs(visibles = flat_images_degen$visibles, params0 = params, C = const*variance_params$C/N, C_prime = const*variance_params$C_prime/N, mc.iter = mc.iter)
  
  #burnin + thinning
  burnin_thinning <- function(models) {
    idx <- which(1:mc.iter %% 5 == 0 & 1:mc.iter > mc.iter/2)
    models$params$main_hidden <- models$params$main_hidden[idx, ]
    models$params$main_visible <- models$params$main_visible[idx, ]
    models$params$interaction <- models$params$interaction[idx, , ]
    models$hiddens <- models$hiddens[idx, , ]
    #models$post_pred <- models$post_pred[idx, , ]
    models$distn <- models$distn[idx, , ]
    
    return(models)
  }
  
  models_good <- burnin_thinning(models_good)
  models_bad <- burnin_thinning(models_bad)
  #models_degen <- burnin_thinning(models_degen)
  
  save(models_good, models_bad, file = paste0("written/fitted_models_jing_", const, ".Rdata"))
  
  
}
