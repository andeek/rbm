#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("../../4x4 model/sample_functs.R")

#data and params ------------------------
H <- 4
V <- 4

#marginalized likelihood
load("../../4x4 model/written/fitted_models_adaptive_mh_trunc_distn_marginal_1.Rdata")
marginal_bad_1 <- models_bad
marginal_good_1 <- models_good
marginal_degen_1 <- models_degen

#load Jing's prior
load("../../4x4 model/written/fitted_models_distn_0.26.Rdata")
jing_bad <- models_bad
jing_good <- models_good
jing_degen <- models_degen

#truncated normal
load("../../4x4 model/written/fitted_models_adaptive_mh_trunc_distn_shrunk_jing_match.Rdata")
shrunk_jing_bad <- models_bad
shrunk_jing_good <- models_good
shrunk_jing_degen <- models_degen

#rm unneccesary data
rm(models_bad)
rm(models_degen)
rm(models_good)

load("../../4x4 model/written/sample.params.Rdata")
params_degen <- list(main_hidden = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                     main_visible = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                     interaction = sample.params %>% ungroup() %>% filter(near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))
params_good <- list(main_hidden = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                    main_visible = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                    interaction = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))


load("../../4x4 model/written/sample_images.Rdata")

#computing actual distributions ---------------------
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

marginal_sample_good_1 <- reshape_sample_distn(marginal_good_1)
jing_sample_good <- reshape_sample_distn(jing_good)
shrunk_jing_sample_good <- reshape_sample_distn(shrunk_jing_good) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 

#plots -------------------
marginal_sample_good_1 %>% rename(marginal = prob) %>%
  left_join(jing_sample_good %>% ungroup() %>% select(-image_id) %>% rename(trick = prob)) %>%
  left_join(shrunk_jing_sample_good %>% ungroup() %>% select(-image_id) %>% rename(trunc = prob)) %>%
  ungroup() %>%
  mutate(method = "good", image_id = paste0("image_", image_id)) %>% 
  #select_(.dots = paste0("-v", 1:V)) %>%
  gather(prior, prob, trick, marginal, trunc) %>%
  spread(image_id, prob) -> all_statistics

all_statistics %>%
  gather(statistic, value, -iter, -method, -prior, -starts_with("v")) %>%
  filter(grepl("image", statistic) & !is.na(value)) %>%
  separate(statistic, into = c("junk", "statistic")) %>%
  mutate(statistic = factor(statistic, labels = paste("Image", 1:length(unique(statistic))))) %>%
  left_join(distn_good %>% select(-image_id)) %>% 
  ggplot() +
  geom_line(aes(iter/100, value, colour = prior)) +
  geom_abline(aes(slope = 0, intercept = prob)) +
  facet_grid(~statistic) +
  ylim(c(0,.5)) +
  xlab("MCMC Iteration (in 100s)") + ylab("Posterior Probability") +
  theme(legend.position = "bottom") +
  scale_colour_discrete("Method", labels = c("Marginalized likelihood", "Trick prior", "Truncated Normal prior")) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)) -> p

ggsave("images/image_prediction.pdf",
       plot = p,
       bg = "transparent",
       width = 14.5,
       height = 2.5,
       units = "in")
