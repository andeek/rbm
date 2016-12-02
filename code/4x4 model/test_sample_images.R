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
load("written/sample.params.Rdata")
params_degen <- list(main_hidden = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                     main_visible = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                     interaction = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))
params_good <- list(main_hidden = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                    main_visible = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                    interaction = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))


cdf_visible <- function(main_visible, interaction, hiddens, visibles) {
  #main_visible is a 1xH matrix of visible main effects param values
  #interaction is a HxV matrix of interaction effects param values
  #hiddens is an NxH matrix of hidden node values
  stopifnot(nrow(hiddens) == nrow(visibles))
  
  N <- nrow(hiddens)

  a <- exp(2*(matrix(main_visible, nrow = nrow(hiddens), ncol = length(main_visible), byrow = TRUE) + tcrossprod(hiddens, interaction)))
  f <- apply(a, 1:2, function(b) dbinom(1, 1, b/(1 + b)))
  A <- matrix(runif(N), nrow = N, ncol = ncol(visibles))
  
  res <- A*f*(visibles == -1) + (f + A*(1 - f))*(visibles == 1)

  return(res)
}

pit_test_degen <- cdf_visible(params_degen$main_visible, params_degen$interaction, flat_images_degen$hiddens, flat_images_degen$visibles)
pit_test_good <- cdf_visible(params_good$main_visible, params_good$interaction, flat_images_good$hiddens, flat_images_good$visibles)

##look at results -----------------------
pit_test_good %>%
  data.frame() %>% 
  gather(node, cdf, everything()) %>%
  ggplot() +
  geom_histogram(aes(cdf)) +
  facet_wrap(~node)

pit_test_degen %>%
  data.frame() %>% 
  gather(node, cdf, everything()) %>%
  ggplot() +
  geom_histogram(aes(cdf)) +
  facet_wrap(~node)

##true vs sample distribution of images
distn_degen <- visible_distn(params = params_degen)
distn_good <- visible_distn(params = params_good)

distn_good %>%
  left_join(flat_images_good$visibles %>% data.frame() %>% rename(v1 = X1, v2 = X2, v3 = X3, v4 = X4)) %>%
  rename(true_good = prob) %>%
  group_by(image_id, v1, v2, v3, v4, true_good) %>%
  summarise(est_good = n()/nrow(flat_images_good$visibles)) %>% 
  left_join(distn_degen %>%
              left_join(flat_images_degen$visibles %>% data.frame() %>% rename(v1 = X1, v2 = X2, v3 = X3, v4 = X4)) %>%
              rename(true_degen = prob) %>%
              group_by(image_id, v1, v2, v3, v4, true_degen) %>%
              summarise(est_degen = n()/nrow(flat_images_degen$visibles))) %>%
  gather(type, prob, -image_id, -(v1:v4)) %>%
  separate(type, into = c("what", "type")) %>%
  spread(what, prob) %>%
  ggplot() +
  geom_bar(aes(x = image_id, y = est), stat = "identity") +
  geom_point(aes(x = image_id, y = true), colour = "red") +
  facet_wrap(~type)
 
