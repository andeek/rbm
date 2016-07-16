#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("sample_functs.R")

#params ------------------------
H <- 4
V <- 4

load("written/sample_images.Rdata")
load("written/params_theta.Rdata")
params_degen <- list(main_hidden = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                     main_visible = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                     interaction = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))
params_good <- list(main_hidden = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                    main_visible = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                    interaction = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))

#computing actual distributions --------------
distn_good <- visible_distn(params = params_good)
distn_degen <- visible_distn(params = params_degen)

#data ----------------------
consts <- c(0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5, 5.5, 5.6, 5.7, 5.8, 5.9, 6, 6.1, 6.2, 7, 8)
reshape_sample_distn <- function(model) {
  sample_distn <- model$distn
  dim(sample_distn) <- c(dim(sample_distn)[1]*dim(sample_distn)[2], dim(sample_distn)[3])
  sample_distn %>% data.frame() -> sample_distn
  names(sample_distn) <- names(distn_good)
  
  sample_distn %>%
    group_by(image_id) %>%
    mutate(iter = 1:n()) -> sample_distn
  
  return(sample_distn)
}

res <- data.frame()
for(const in consts) {
  load(paste0("written/fitted_models_jing_", const, ".Rdata"))
  sample_good <- reshape_sample_distn(models_good)
  sample_bad <- reshape_sample_distn(models_bad)
  
  sample_good %>% rename(good = prob) %>%
    left_join(sample_bad %>% rename(bad = prob) %>% ungroup() %>% select(-image_id)) %>%
    left_join(distn_good %>% rename(true_prob = prob) %>% select(-image_id)) %>%
    mutate(const = const) -> samp
  
  res <<- bind_rows(res, samp)
}

#error rates --------------
res %>%
  gather(which, est_prob, good, bad) %>%
  mutate(error = (true_prob - est_prob)^2) %>%
  group_by(which, const) %>%
  summarise(mse = mean(error)) -> error

error %>%
  ungroup() %>%
  mutate(min_mse = min(mse)) %>%
  filter(mse == min_mse)

error %>%
  ggplot() +
  geom_line(aes(const, mse, colour = which)) +
  geom_point(aes(const, mse, colour = which))
