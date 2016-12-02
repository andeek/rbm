#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("sample_functs.R")

#data and params ------------------------
load("written/sample.params.Rdata")
params_good <- list(main_hidden = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                    main_visible = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                    interaction = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))

distn_good <- visible_distn(params = params_good)
possibles <- stats(4, 4)

load("written/fitted_ml_margin_zeros.Rdata")
theta.df <- theta %>% data.frame()
colnames(theta.df) <- colnames(possibles)

theta.df %>%
  rowwise() %>%
  do(params = list(main_hidden = data.frame(.) %>% select(starts_with("h")) %>% data.matrix(),
                    main_visible = data.frame(.) %>% select(starts_with("v")) %>% data.matrix(),
                    interaction = data.frame(.) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))
  ) %>%
  ungroup() %>%
  mutate(iter = 1:n()) %>%
  group_by(iter) %>%
  do(data.frame(visible_distn(.$params[[1]]))) -> distn_ml
 
distn_good %>% rename(true = prob) %>%
  left_join(distn_ml %>% rename(est = prob) %>% ungroup %>% select(-image_id)) %>% 
  #gather(variable, value, true, est) %>%
  ggplot() +
  geom_line(aes(iter, est), colour = 'blue') +
  geom_line(aes(iter, true)) +
  facet_wrap(~image_id)

#likelihood-----
len <- 20

expand.grid(V1 = seq(-10, 10, length.out = len), V2 = seq(-10, 10, length.out = len)) %>%
  cbind(matrix(rep(rnorm(H + V + H*V - 2, 0, 10), len*len), nrow = len*len, byrow = TRUE)) %>%
  data.frame() %>% 
  group_by(V1, V2) %>%
  do(data.frame(loglik = loglik(as.numeric(.), parms))) -> loglik_grid

loglik_grid %>%
  ggplot() +
  geom_contour(aes(x = V1, y = V2, z = loglik, colour = ..level..))


theta_select <- theta[100, ] %>% matrix(nrow = 1)
colnames(theta_select) <- colnames(theta_good)

