#libraries --------------------
library(dplyr)
library(tidyr)
source("sample_functs_4x4.R")

#data and params ------------------------
H <- 4
V <- 4
N <- 1

load("written/sample.params.Rdata")
params_degen <- list(main_hidden = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                     main_visible = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                     interaction = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4, byrow = TRUE))
params_good <- list(main_hidden = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                    main_visible = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                    interaction = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4, byrow = TRUE))

hidden0 <- matrix(sample(c(-1,1), H*N, replace = TRUE), nrow = N)
visible0 <- matrix(sample(c(-1,1), V*N, replace = TRUE), nrow = N)

flat_images_degen <- sample_images(params = params_degen, hidden0 = hidden0, visible0 = visible0, mc.iter = 5e4)
flat_images_good <- sample_images(params = params_good, hidden0 = hidden0, visible0 = visible0, mc.iter = 5e4)

for(i in 1:2) {
  idx <- 1:2.5e4 %% 5 == 0
  flat_images_degen[[i]] <- flat_images_degen[[i]][-(1:(2.5e4 + 1)), ][idx, ] #burnin + thinning
  flat_images_good[[i]] <- flat_images_good[[i]][-(1:(2.5e4 + 1)), ][idx, ] #burnin + thinning
}

save(flat_images_good, flat_images_degen, file = "written/sample_images.Rdata")
