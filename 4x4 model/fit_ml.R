#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("sample_functs.R")
source("ml_functs.R")

#data and params ------------------------
H <- 4
V <- 4
load("written/sample_images.Rdata")
load("written/sample.params.Rdata")
set.seed(102285)


#system of nonlinear eqns to solve -----------------
parms <- list(visibles = flat_images_good$visibles, 
              hiddens = flat_images_good$hiddens,
              H = H, V = V)

max.iter <- 10
theta <- matrix(NA, nrow = max.iter + 1, ncol = H + V + H*V)
likelihood <- matrix(NA, nrow = max.iter + 1, ncol = 1)

#inits
theta[1, ] <- sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("v"), starts_with("h"), starts_with("theta"), -H, -V) %>% data.matrix()
likelihood[1, ] <- loglik(theta = theta[1, ], parms = parms)

for(i in 1:max.iter) {
  cat(paste0(i, "\r"))
  optim(par = theta[i, ],
        fn = loglik,
        #gr = loglik_derivs,
        parms = parms,
        method = "BFGS",
        control = list(fnscale = -1, maxit = 1)) -> mle
  
  if(mle$convergence == 1) {
    theta[i + 1, ] <- mle$par
    likelihood[i + 1, ] <- mle$value
  } else {
    return
  }
}

save(theta, likelihood, file = "written/fitted_ml_margin.Rdata")


cbind(theta, likelihood) %>%
  data.frame() %>%
  mutate(iter = 1:(max.iter + 1)) %>%
  rename(likelihood = X25) %>%
  gather(variable, value, -iter) %>%
  ggplot() +
  geom_line(aes(x = iter, y = value)) +
  facet_wrap(~variable, scales = "free_y") 








