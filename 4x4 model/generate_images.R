#libraries --------------------
library(dplyr)
library(tidyr)
source("sample_functs.R")

#data and params ------------------------
H <- 4
V <- 4
N <- 1
mc.iter <- 10e4
set.seed(102285)

load("written/sample.params.Rdata")
params_degen <- list(main_hidden = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                     main_visible = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                     interaction = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))
params_good <- list(main_hidden = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                    main_visible = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                    interaction = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))

hidden0 <- matrix(sample(c(-1,1), H*N, replace = TRUE), nrow = N)
visible0 <- matrix(sample(c(-1,1), V*N, replace = TRUE), nrow = N)

flat_images_degen <- sample_images(params = params_degen, hidden0 = hidden0, visible0 = visible0, mc.iter = mc.iter)
flat_images_good <- sample_images(params = params_good, hidden0 = hidden0, visible0 = visible0, mc.iter = mc.iter)

for(i in 1:2) {
  idx <- 1:(mc.iter/2) %% 10 == 0
  flat_images_degen[[i]] <- flat_images_degen[[i]][-(1:(mc.iter/2 + 1)), ][idx, ] #burnin + thinning
  flat_images_good[[i]] <- flat_images_good[[i]][-(1:(mc.iter/2 + 1)), ][idx, ] #burnin + thinning
}

save(flat_images_good, flat_images_degen, file = "written/sample_images.Rdata")
