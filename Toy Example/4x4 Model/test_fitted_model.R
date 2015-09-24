#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("sample_functs.R")

#data and params ------------------------
H <- 4
V <- 4
N <- 1

load("written/sample_images.Rdata")
load("written/fitted_models_distn_1.Rdata")
load("written/sample.params.Rdata")
params_degen <- list(main_hidden = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                     main_visible = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                     interaction = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))
params_good <- list(main_hidden = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                    main_visible = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                    interaction = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))


##plots -----------------------

#computing actual distributions
distn_good <- visible_distn(params = params_good)
distn_degen <- visible_distn(params = params_degen)

get_props <- function(model) {
  props <- apply(model$post_pred, 1, function(x) {
    x %>% 
      data.frame %>%
      rename(v1 = X1, v2 = X2, v3 = X3, v4 = X4) %>%
      right_join(distn_good) %>% 
      group_by(image_id) %>%
      summarise(prop = n()/nrow(x))
  })
  
  props %>% 
    bind_rows(.id = "iteration") %>% 
    mutate(image_id = as.numeric(image_id), iteration = as.numeric(iteration)) %>%
    data.frame
}

props_good <- get_props(models_good)
props_bad <- get_props(models_bad)
props_degen <- get_props(models_degen)

data_props_good <- flat_images_good$visibles %>% 
  data.frame %>%
  rename(v1 = X1, v2 = X2, v3 = X3, v4 = X4) %>%
  right_join(distn_good) %>%
  group_by(image_id) %>%
  summarise(prop = n()/nrow(flat_images_good$visibles))

data_props_degen <- flat_images_degen$visibles %>% 
  data.frame %>%
  rename(v1 = X1, v2 = X2, v3 = X3, v4 = X4) %>%
  right_join(distn_degen) %>%
  group_by(image_id) %>%
  summarise(prop = n()/nrow(flat_images_degen$visibles)) 


props_good %>% rename(good = prop) %>%
  left_join(props_bad %>% rename(bad = prop)) %>%
  gather(model, prop, good, bad) %>%
  ggplot() +
  geom_histogram(aes(prop)) + 
  facet_grid(model~image_id) +
  xlim(c(0,1)) +
  geom_segment(aes(x = prop, xend = prop, y = 0, yend = 500), colour = "red", data = data_props_good) +
  geom_segment(aes(x = prob, xend = prob, y = 0, yend = 500), colour = "blue", data = distn_good)
  

props_degen %>%
  ggplot() +
  geom_histogram(aes(prop)) + 
  geom_segment(aes(x = prop, xend = prop, y = 0, yend = 500), colour = "red", data = data_props_degen) +
  facet_grid(.~image_id) +
  xlim(c(0,1)) +
  geom_segment(aes(x = prob, xend = prob, y = 0, yend = 500), colour = "blue", data = distn_degen)

########################################################
reshape_sample_distn <- function(model) {
  sample_distn <- model$distn
  dim(sample_distn) <- c(dim(sample_distn)[1]*dim(sample_distn)[2], dim(sample_distn)[3])
  sample_distn %>% data.frame() -> sample_distn
  names(sample_distn) <- names(distn_degen)
  
  sample_distn %>%
    group_by(image_id) %>%
    mutate(iter = 1:n()) -> sample_distn
  
  return(sample_distn)
}

sample_good <- reshape_sample_distn(models_good)
sample_bad <- reshape_sample_distn(models_bad)
sample_degen <- reshape_sample_distn(models_degen)


sample_good %>% rename(good = prob) %>%
  left_join(sample_bad %>% rename(bad = prob)) %>%
  gather(model, prob, good, bad) %>%
  ggplot() +
  geom_histogram(aes(prob)) +
  facet_grid(model~image_id) +
  geom_segment(aes(x = prob, xend = prob, y = 0, yend = 500), colour = "blue", data = distn_good) +
  xlim(c(0,1))


