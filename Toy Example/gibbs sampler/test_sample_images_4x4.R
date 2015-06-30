#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
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


cdf_visible <- function(main_visible, interaction, hiddens, visibles) {
  #main_visible is a 1xH matrix of visible main effects param values
  #interaction is a HxV matrix of interaction effects param values
  #hiddens is an NxH matrix of hidden node values
  stopifnot(nrow(hiddens) == nrow(visibles))
  
  N <- nrow(hiddens)
  res <- matrix(NA, nrow = N, ncol = ncol(visibles))
  
  
  for(i in 1:N) {
    a <- exp(2*(sum(main_visible) + rowSums(t(colSums(interaction) * t(matrix(hiddens[i, ], nrow = 1))))))
    A <- runif(1)
    f <- dbinom(0, 1, a/(a + 1))
    res[i, ] <- A*f*(visibles[i, ] == -1) + (f + A*(1 - f))*(visibles[i, ] == 1)
  }

  return(res)
}

pit_test_degen <- cdf_visible(params_degen$main_visible, params_degen$interaction, flat_images_degen$hiddens, flat_images_degen$visibles)
pit_test_good <- cdf_visible(params_good$main_visible, params_good$interaction, flat_images_good$hiddens, flat_images_good$visibles)

##look at results -----------------------
pit_test_good %>%
  data.frame() %>% 
  gather(node, cdf, everything()) %>%
  ggplot() +
  geom_histogram(aes(cdf), binwidth = .01) +
  facet_wrap(~node)

pit_test_degen %>%
  data.frame() %>% 
  gather(node, cdf, everything()) %>%
  ggplot() +
  geom_histogram(aes(cdf), binwidth = .01) +
  facet_wrap(~node)
  




