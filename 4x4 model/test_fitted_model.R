#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("sample_functs.R")

#data and params ------------------------
H <- 4
V <- 4

load("written/sample_images.Rdata")
# load("written/fitted_models_adaptive_mh_trunc_distn.Rdata")
load("written/fitted_models_adaptive_mh_full_trunc_distn_marginal_2.Rdata")
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

# get_props <- function(model) {
#   props <- apply(model$post_pred, 1, function(x) {
#     x %>% 
#       data.frame %>%
#       rename(v1 = X1, v2 = X2, v3 = X3, v4 = X4) %>%
#       right_join(distn_good) %>% 
#       group_by(image_id) %>%
#       summarise(prop = n()/nrow(x))
#   })
#   
#   props %>% 
#     bind_rows(.id = "iteration") %>% 
#     mutate(image_id = as.numeric(image_id), iteration = as.numeric(iteration)) %>%
#     data.frame
# }
# 
# props_good <- get_props(models_good)
# props_bad <- get_props(models_bad)
# props_degen <- get_props(models_degen)
# 
# data_props_good <- flat_images_good$visibles %>% 
#   data.frame %>%
#   rename(v1 = X1, v2 = X2, v3 = X3, v4 = X4) %>%
#   right_join(distn_good) %>%
#   group_by(image_id) %>%
#   summarise(prop = n()/nrow(flat_images_good$visibles))
# 
# data_props_degen <- flat_images_degen$visibles %>% 
#   data.frame %>%
#   rename(v1 = X1, v2 = X2, v3 = X3, v4 = X4) %>%
#   right_join(distn_degen) %>%
#   group_by(image_id) %>%
#   summarise(prop = n()/nrow(flat_images_degen$visibles)) 
# 
# 
# props_good %>% rename(good = prop) %>%
#   left_join(props_bad %>% rename(bad = prop)) %>%
#   gather(model, prop, good, bad) %>%
#   ggplot() +
#   geom_histogram(aes(prop)) + 
#   facet_grid(model~image_id) +
#   xlim(c(0,1)) +
#   geom_segment(aes(x = prop, xend = prop, y = 0, yend = 500), colour = "red", data = data_props_good) +
#   geom_segment(aes(x = prob, xend = prob, y = 0, yend = 500), colour = "blue", data = distn_good)
#   
# 
# props_degen %>%
#   ggplot() +
#   geom_histogram(aes(prop)) + 
#   geom_segment(aes(x = prop, xend = prop, y = 0, yend = 500), colour = "red", data = data_props_degen) +
#   facet_grid(.~image_id) +
#   xlim(c(0,1)) +
#   geom_segment(aes(x = prob, xend = prob, y = 0, yend = 500), colour = "blue", data = distn_degen)

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
  left_join(sample_bad %>% rename(bad = prob) %>% ungroup() %>% select(-image_id)) %>%
  left_join(distn_good %>% rename(true_prob = prob) %>% select(-image_id)) %>%
  gather(model, prob, good, bad) %>%
  ggplot() +
  geom_histogram(aes(prob)) +
  facet_grid(model~image_id) +
  geom_segment(aes(x = true_prob, xend = true_prob, y = 0, yend = max(sample_good$iter)), colour = "blue") +
  xlim(c(0,1))

sample_degen %>%
  left_join(distn_degen %>% rename(true_prob = prob) %>% select(-image_id)) %>%
  ggplot() +
  geom_histogram(aes(prob)) +
  facet_grid(~image_id) +
  geom_segment(aes(x = true_prob, xend = true_prob, y = 0, yend = max(sample_good$iter)), colour = "blue") +
  xlim(c(0,1))

sample_good %>% rename(good = prob) %>%
  left_join(sample_bad %>% rename(bad = prob)) %>%
  left_join(distn_good) %>%
  gather(model, est_prob, good, bad) %>%
  mutate(p_val = est_prob > prob) %>%
  group_by(model, image_id) %>%
  summarise(p_val = sum(p_val)/n()) %>%
  spread(model, p_val)

sample_good %>% rename(good = prob) %>%
  left_join(sample_bad %>% rename(bad = prob) %>% ungroup() %>% select(-image_id)) %>%
  left_join(distn_good %>% rename(true_prob = prob) %>% select(-image_id)) %>%
  gather(model, prob, good, bad) %>%
  ggplot() +
  geom_line(aes(iter, prob, colour = model)) +
  geom_abline(aes(intercept = true_prob, slope = 0)) +
  #geom_abline(aes(intercept = prop, slope = 0), data = data_props_good, colour = "red") +
  facet_wrap(~image_id) +
  ylim(c(0,1))

sample_degen %>% 
  left_join(distn_degen %>% rename(true_prob = prob) %>% select(-image_id)) %>%
  ggplot() +
  geom_line(aes(iter, prob)) +
  geom_abline(aes(intercept = true_prob, slope = 0), colour = "blue") +
  #geom_abline(aes(intercept = prop, slope = 0), data = data_props_good, colour = "red") +
  facet_wrap(~image_id) +
  ylim(c(0,1))


possibles <- stats(4, 4)

data.frame(good = colSums(exp(possibles %*% t(models_good$theta))),
           bad = colSums(exp(possibles %*% t(models_bad$theta)))) %>%
  mutate(iter = 1:n()) %>%
  gather(model, normalizer, -iter) %>%
  ggplot() +
  geom_line(aes(iter, normalizer, colour = model))

  


