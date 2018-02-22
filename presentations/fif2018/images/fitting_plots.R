#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("../../../writing/resources/code/functions.R")

load("../../../writing/resources/data/sample_images.Rdata")

#data and params ------------------------
H <- 4
V <- 4

load("../../../writing/resources/data/params_theta.Rdata")

params_good <- list(main_hidden = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                    main_visible = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                    interaction = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))


#marginalized likelihood
load("../../../writing/resources/data/fitted_models_trunc_marginal_full.Rdata")
marginal_bad <- models_bad
marginal_good <- models_good

#load trick prior
load("../../../writing/resources/data/fitted_models_jing_5.8.Rdata")
trick_bad <- models_bad
trick_good <- models_good

#truncated normal
load("../../../writing/resources/data/fitted_models_trunc_full.Rdata")
trunc_bad <- models_bad
trunc_good <- models_good

#rm unneccesary data
rm(models_bad)
rm(models_good)

#computing actual distributions ---------------------
flat_images_good$visibles %>%
  data.frame() %>% 
  rename(v1 = X1, v2 = X2, v3 = X3, v4 = X4) %>%
  group_by(v1, v2, v3, v4) %>%
  summarise(prob = n()/nrow(flat_images_good$visibles)) -> distn_emp

distn_good <- visible_distn(params = params_good)

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

marginal_sample_good <- reshape_sample_distn(marginal_good)
trick_sample_good <- reshape_sample_distn(trick_good)
trunc_sample_good <- reshape_sample_distn(trunc_good)


#plots -------------------
marginal_sample_good %>% 
  filter(iter > 50) %>%
  group_by(v1, v2, v3, v4, image_id) %>% 
  do(data.frame(acf = acf(.$prob, plot=FALSE)$acf,
                lag = acf(.$prob, plot=FALSE)$lag)) -> marginal_acfs

trunc_sample_good %>% 
  filter(iter > 50) %>%
  group_by(v1, v2, v3, v4, image_id) %>% 
  do(data.frame(acf = acf(.$prob, plot=FALSE)$acf,
                lag = acf(.$prob, plot=FALSE)$lag)) -> trunc_acfs

marginal_acfs %>% ungroup() %>% select(-image_id) %>% rename(marginal = acf) %>%
  left_join(trunc_acfs %>% rename(truncated = acf)) %>%
  ggplot() +
  geom_segment(aes(x=lag, xend=lag, y=0, yend=truncated), colour = "black") +
  geom_segment(aes(x=lag, xend=lag, y=0, yend=marginal), colour = "red") +
  geom_point(aes(x=lag, y=marginal), colour = "red", size = 0.5) +
  facet_grid(~image_id) +
  xlab("Lag") +
  ylab("ACF") +
  theme(legend.position = "bottom") +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)) -> acf_plot

marginal_sample_good %>% rename(marginal = prob) %>% filter(iter > 50) %>%
  left_join(trick_sample_good %>% ungroup() %>% select(-image_id) %>% rename(trick = prob)) %>%
  ungroup() %>%
  mutate(image_id = paste0("image_", image_id)) %>% 
  gather(method, prob, trick, marginal) %>%
  spread(image_id, prob) -> all_statistics

all_statistics %>%
  gather(statistic, value, -iter, -method, -starts_with("v")) %>%
  filter(grepl("image", statistic) & !is.na(value)) %>%
  separate(statistic, into = c("junk", "statistic")) %>%
  mutate(statistic = factor(statistic, labels = paste("Vector", 1:length(unique(statistic))))) %>%
  left_join(distn_emp %>% rename(prob_emp = prob) %>% right_join(distn_good %>% select(-image_id))) %>%
  ggplot() +
  geom_density(aes(value, y=..scaled.., colour = method, fill = method), alpha = .2) +
  geom_vline(aes(xintercept = prob)) +
  geom_vline(aes(xintercept = prob_emp), lty = 2) +
  facet_grid(~statistic, scales = "free_x") + 
  ylab("Scaled Posterior Density") + xlab("Probability of Vector of Visibles") +
  scale_colour_discrete("Method", labels = c("BwTNML", "BwTPLV")) +
  scale_fill_discrete("Method", labels = c("BwTNML", "BwTPLV")) +
  theme(legend.position = "bottom") +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)) -> p.models

sx <- scale_x_continuous()
sx$trans$breaks <- function(range) pretty(range, n = 3)

p <- p.models + sx

ggsave("image_acf.pdf",
       plot = acf_plot,
       bg = "transparent",
       width = 14,
       height = 2.5,
       units = "in")

ggsave("image_prediction.pdf",
       plot = p,
       bg = "transparent",
       width = 14.5,
       height = 2.5,
       units = "in")
