#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("sample_functs.R")

#data and params ------------------------
H <- 4
V <- 4

#load truncatednormal prior
load("written/fitted_models_adaptive_mh_trunc_distn_shrunk_0.6667.Rdata")
shrunk_bad <- models_bad
shrunk_degen <- models_degen
shrunk_good <- models_good

#load Jing's prior
load("written/fitted_models_distn_0.26.Rdata")
jing_bad <- models_bad
jing_degen <- models_degen
jing_good <- models_good

#rm unneccesary data
rm(models_bad)
rm(models_degen)
rm(models_good)

load("written/sample.params.Rdata")
params_degen <- list(main_hidden = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                     main_visible = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                     interaction = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))
params_good <- list(main_hidden = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                    main_visible = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                    interaction = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))



#computing actual distributions
distn_good <- visible_distn(params = params_good)
distn_degen <- visible_distn(params = params_degen)

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

shrunk_sample_good <- reshape_sample_distn(shrunk_good)
shrunk_sample_bad <- reshape_sample_distn(shrunk_bad)
shrunk_sample_degen <- reshape_sample_distn(shrunk_degen)

jing_sample_good <- reshape_sample_distn(jing_good)
jing_sample_bad <- reshape_sample_distn(jing_bad)
jing_sample_degen <- reshape_sample_distn(jing_degen)

#plots -------------------
jing_sample_good %>% rename(jing = prob) %>%
  left_join(shrunk_sample_good %>% ungroup() %>% select(-image_id) %>% rename(shrunk = prob)) %>%
  left_join(distn_good %>% ungroup() %>% select(-image_id) %>% rename(true = prob)) %>%
  gather(model, prob, jing, shrunk) %>%
  ggplot() +
  geom_line(aes(iter, prob, colour = model)) +
  geom_abline(aes(slope = 0, intercept = true)) +
  facet_wrap(~image_id) +
  ylim(c(0,1))
