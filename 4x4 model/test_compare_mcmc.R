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
shrunk_good <- models_good
shrunk_degen <- models_degen

load("written/fitted_models_adaptive_mh_trunc_distn_marginal_1.Rdata")
marginal_bad_1 <- models_bad
marginal_good_1 <- models_good
marginal_degen_1 <- models_degen


load("written/fitted_models_adaptive_mh_trunc_distn_marginal_2.Rdata")
marginal_bad_2 <- models_bad
marginal_good_2 <- models_good
marginal_degen_2 <- models_degen

load("written/fitted_models_adaptive_mh_trunc_distn_marginal_3.Rdata")
marginal_bad_3 <- models_bad
marginal_good_3 <- models_good
marginal_degen_3 <- models_degen

#load Jing's prior
load("written/fitted_models_distn_0.26.Rdata")
jing_bad <- models_bad
jing_good <- models_good
jing_degen <- models_degen

load("written/fitted_models_adaptive_mh_trunc_distn_shrunk_jing_match.Rdata")
shrunk_jing_bad <- models_bad
shrunk_jing_good <- models_good
shrunk_jing_degen <- models_degen

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


load("written/sample_images.Rdata")

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

#burn in more on the truncated, there are more iterations than the Jing, which gets complicated with plotting
shrunk_sample_good <- reshape_sample_distn(shrunk_good) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_sample_bad <- reshape_sample_distn(shrunk_bad) %>% filter(iter > 500) %>% mutate(iter = 1:n())
shrunk_sample_degen <- reshape_sample_distn(shrunk_degen) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 

marginal_sample_good_1 <- reshape_sample_distn(marginal_good_1)
marginal_sample_bad_1 <- reshape_sample_distn(marginal_bad_1)
marginal_sample_degen_1 <- reshape_sample_distn(marginal_degen_1)

marginal_sample_good_2 <- reshape_sample_distn(marginal_good_2)
marginal_sample_bad_2 <- reshape_sample_distn(marginal_bad_2)
marginal_sample_degen_2 <- reshape_sample_distn(marginal_degen_2)

marginal_sample_good_3 <- reshape_sample_distn(marginal_good_3)
marginal_sample_bad_3 <- reshape_sample_distn(marginal_bad_3)
marginal_sample_degen_3 <- reshape_sample_distn(marginal_degen_3)

jing_sample_good <- reshape_sample_distn(jing_good)
jing_sample_bad <- reshape_sample_distn(jing_bad)
jing_sample_degen <- reshape_sample_distn(jing_degen)

shrunk_jing_sample_good <- reshape_sample_distn(shrunk_jing_good) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_jing_sample_bad <- reshape_sample_distn(shrunk_jing_bad) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_jing_sample_degen <- reshape_sample_distn(shrunk_jing_degen) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 

#computing loglikelihoods, ratios ---------------------
reshape_loglik <- function(model, visibles) {
  if("theta" %in% names(model)) {
    theta <- model$theta
  } else {
    interaction <- t(apply(model$params$interaction, 1, t))
    theta <- cbind(model$params$main_visible, model$params$main_hidden, interaction)
  }
  #hiddens <- model$hiddens
  N <- visibles %>% nrow()
  H <- dim(model$hiddens)[3]
  V <- ncol(visibles)
  
  loglik <- do.call(rbind, lapply(seq_along(theta[, 1]), function(m) {
    
    possibles <- stats(H, V)
    data <- visibles %>% data.frame
    names(data) <- paste0("v", 1:V)
    
    data <- data %>%
      left_join(possibles %>% data.frame()) %>% data.frame
    
    N <- nrow(data)
    require(dplyr)
    N <- nrow(data)
    data <- data %>% data.frame() 
    names(data) <- colnames(possibles)
    H <- data %>% select(starts_with("h")) %>% ncol
    V <- data %>% select(starts_with("V")) %>% ncol
    
    data %>%
      select(starts_with("v")) %>%
      left_join(possibles %>% data.frame, by = paste0("v", 1:V)) -> data
    
    data %>%
      group_by_(.dots = names(data)) %>%
      summarise(count = n()) %>%
      group_by(count, add = TRUE) %>%
      do(data.frame(A = exp(as.numeric(data.frame(.) %>% select(-count)) %*% theta[m, ]))) %>% 
      mutate(A_count = A*count) %>%
      group_by_(.dots = paste0("v", 1:V)) %>%
      summarise(g_theta = sum(A_count), count = sum(count)/2^(H)) -> inner
 
    data.frame(iter = m, 
               log_normalizer = N*log(sum(exp(possibles %*% theta[m, ]))),
               loglik = -N*log(sum(exp(possibles %*% theta[m, ]))) + sum(inner$count*log(inner$g_theta)))
  }))
  
  loglik
}

shrunk_loglik_good <- reshape_loglik(shrunk_good, flat_images_good$visibles) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_loglik_bad <- reshape_loglik(shrunk_bad, flat_images_good$visibles) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_loglik_degen <- reshape_loglik(shrunk_degen, flat_images_good$visibles) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 

marginal_loglik_good_1 <- reshape_loglik(marginal_good_1, flat_images_good$visibles)
marginal_loglik_bad_1 <- reshape_loglik(marginal_bad_1, flat_images_good$visibles)
marginal_loglik_degen_1 <- reshape_loglik(marginal_degen_1, flat_images_good$visibles)

marginal_loglik_good_2 <- reshape_loglik(marginal_good_2, flat_images_good$visibles)
marginal_loglik_bad_2 <- reshape_loglik(marginal_bad_2, flat_images_good$visibles)
marginal_loglik_degen_2 <- reshape_loglik(marginal_degen_2, flat_images_good$visibles)

marginal_loglik_good_3 <- reshape_loglik(marginal_good_3, flat_images_good$visibles)
marginal_loglik_bad_3 <- reshape_loglik(marginal_bad_3, flat_images_good$visibles)
marginal_loglik_degen_3 <- reshape_loglik(marginal_degen_3, flat_images_good$visibles)

jing_loglik_good <- reshape_loglik(jing_good, flat_images_good$visibles)
jing_loglik_bad <- reshape_loglik(jing_bad, flat_images_good$visibles)
jing_loglik_degen <- reshape_loglik(jing_degen, flat_images_good$visibles)

shrunk_jing_loglik_good <- reshape_loglik(shrunk_jing_good, flat_images_good$visibles) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_jing_loglik_bad <- reshape_loglik(shrunk_jing_bad, flat_images_good$visibles) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_jing_loglik_degen <- reshape_loglik(shrunk_jing_degen, flat_images_good$visibles) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 

reshape_distance <- function(model) {
  if("theta" %in% names(model)) {
    theta <- model$theta
  } else {
    interaction <- t(apply(model$params$interaction, 1, t))
    theta <- cbind(model$params$main_visible, model$params$main_hidden, interaction)
  }
  
  theta %>%
    data.frame() %>%
    mutate_(distance_main = paste0("sqrt(", paste(paste0("X", 1:(V + H), "^2"), collapse = " + "), ")"),
            distance_interaction = paste0("sqrt(", paste(paste0("X", (V + H + 1):(V + H + V*H), "^2"), collapse = " + "), ")"),
            distance_combined = paste0("sqrt(", paste(paste0("X", 1:(V + H + H*V), "^2"), collapse = " + "), ")")) %>%
    mutate(iter = 1:n()) %>%
    select(iter, contains("distance"))
 
}

shrunk_distance_good <- reshape_distance(shrunk_good) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_distance_bad <- reshape_distance(shrunk_bad) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_distance_degen <- reshape_distance(shrunk_degen) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 

marginal_distance_good_1 <- reshape_distance(marginal_good_1)
marginal_distance_bad_1 <- reshape_distance(marginal_bad_1)
marginal_distance_degen_1 <- reshape_distance(marginal_degen_1)

marginal_distance_good_2 <- reshape_distance(marginal_good_2)
marginal_distance_bad_2 <- reshape_distance(marginal_bad_2)
marginal_distance_degen_2 <- reshape_distance(marginal_degen_2)

marginal_distance_good_3 <- reshape_distance(marginal_good_3)
marginal_distance_bad_3 <- reshape_distance(marginal_bad_3)
marginal_distance_degen_3 <- reshape_distance(marginal_degen_3)

jing_distance_good <- reshape_distance(jing_good)
jing_distance_bad <- reshape_distance(jing_bad)
jing_distance_degen <- reshape_distance(jing_degen)

shrunk_jing_distance_good <- reshape_distance(shrunk_jing_good) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_jing_distance_bad <- reshape_distance(shrunk_jing_bad) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_jing_distance_degen <- reshape_distance(shrunk_jing_degen) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 

#plots -------------------
jing_sample_good %>% rename(jing = prob) %>%
  left_join(marginal_sample_good_1 %>% ungroup() %>% select(-image_id) %>% rename(marginal_1 = prob)) %>%
  left_join(marginal_sample_good_2 %>% ungroup() %>% select(-image_id) %>% rename(marginal_2 = prob)) %>%
  left_join(marginal_sample_good_3 %>% ungroup() %>% select(-image_id) %>% rename(marginal_3 = prob)) %>%
  left_join(shrunk_sample_good %>% ungroup() %>% select(-image_id) %>% rename(shrunk = prob)) %>%
  left_join(shrunk_jing_sample_good %>% ungroup() %>% select(-image_id) %>% rename(shrunk_jing = prob)) %>%
  ungroup() %>%
  mutate(method = "good", image_id = paste0("image_", image_id)) %>% 
  rbind_list(jing_sample_bad %>% rename(jing = prob) %>%
              left_join(marginal_sample_bad_1 %>% ungroup() %>% select(-image_id) %>% rename(marginal_1 = prob)) %>%
              left_join(marginal_sample_bad_2 %>% ungroup() %>% select(-image_id) %>% rename(marginal_2 = prob)) %>%
               left_join(marginal_sample_bad_3 %>% ungroup() %>% select(-image_id) %>% rename(marginal_3 = prob)) %>%
              left_join(shrunk_sample_bad %>% ungroup() %>% select(-image_id) %>% rename(shrunk = prob)) %>%
              left_join(shrunk_jing_sample_bad %>% ungroup() %>% select(-image_id) %>% rename(shrunk_jing = prob)) %>%
              ungroup() %>% 
              mutate(method = "bad", image_id = paste0("image_", image_id))) %>%
  select_(.dots = paste0("-v", 1:V)) %>%
  gather(prior, prob, jing, marginal_1, marginal_2, marginal_3, shrunk, shrunk_jing) %>%
  spread(image_id, prob) %>%
  left_join(jing_loglik_good %>% mutate(method = "good") %>%
              rbind_list(jing_loglik_bad %>% mutate(method = "bad")) %>%
              mutate(prior = "jing") %>%
              rbind_list(marginal_loglik_good_1 %>% mutate(method = "good") %>%
                           rbind_list(marginal_loglik_bad_1 %>% mutate(method = "bad")) %>%
                           mutate(prior = "marginal_1")) %>%
              rbind_list(marginal_loglik_good_2 %>% mutate(method = "good") %>%
                           rbind_list(marginal_loglik_bad_2 %>% mutate(method = "bad")) %>%
                           mutate(prior = "marginal_2")) %>%
              rbind_list(marginal_loglik_good_3 %>% mutate(method = "good") %>%
                           rbind_list(marginal_loglik_bad_3 %>% mutate(method = "bad")) %>%
                           mutate(prior = "marginal_3")) %>%
              rbind_list(shrunk_loglik_good %>% mutate(method = "good") %>%
                           rbind_list(shrunk_loglik_bad %>% mutate(method = "bad")) %>%
                           mutate(prior = "shrunk")) %>%
              rbind_list(shrunk_jing_loglik_good %>% mutate(method = "good") %>%
                           rbind_list(shrunk_jing_loglik_bad %>% mutate(method = "bad")) %>%
                           mutate(prior = "shrunk_jing"))) %>%
  left_join(jing_distance_good %>% mutate(method = "good") %>%
              rbind_list(jing_distance_bad %>% mutate(method = "bad")) %>%
              mutate(prior = "jing") %>%
              rbind_list(marginal_distance_good_1 %>% mutate(method = "good") %>%
                           rbind_list(marginal_distance_bad_1 %>% mutate(method = "bad")) %>%
                           mutate(prior = "marginal_1")) %>%
              rbind_list(marginal_distance_good_2 %>% mutate(method = "good") %>%
                           rbind_list(marginal_distance_bad_2 %>% mutate(method = "bad")) %>%
                           mutate(prior = "marginal_2")) %>%
              rbind_list(marginal_distance_good_3 %>% mutate(method = "good") %>%
                           rbind_list(marginal_distance_bad_3 %>% mutate(method = "bad")) %>%
                           mutate(prior = "marginal_3")) %>%
              rbind_list(shrunk_distance_good %>% mutate(method = "good") %>%
                           rbind_list(shrunk_distance_bad %>% mutate(method = "bad")) %>%
                           mutate(prior = "shrunk")) %>%
              rbind_list(shrunk_jing_distance_good %>% mutate(method = "good") %>%
                           rbind_list(shrunk_jing_distance_bad %>% mutate(method = "bad")) %>%
                           mutate(prior = "shrunk_jing"))
  ) -> all_statistics

jing_sample_degen %>% rename(jing = prob) %>%
  left_join(marginal_sample_degen_1 %>% ungroup() %>% select(-image_id) %>% rename(marginal_1 = prob)) %>%
  left_join(marginal_sample_degen_2 %>% ungroup() %>% select(-image_id) %>% rename(marginal_2 = prob)) %>%
  left_join(marginal_sample_degen_3 %>% ungroup() %>% select(-image_id) %>% rename(marginal_3 = prob)) %>%
  left_join(shrunk_sample_degen %>% ungroup() %>% select(-image_id) %>% rename(shrunk = prob)) %>%
  left_join(shrunk_jing_sample_degen %>% ungroup() %>% select(-image_id) %>% rename(shrunk_jing = prob)) %>%
  ungroup() %>%
  mutate(image_id = paste0("image_", image_id)) %>% 
  select_(.dots = paste0("-v", 1:V)) %>%
  gather(prior, prob, jing, marginal_1, marginal_2, marginal_3, shrunk, shrunk_jing) %>%
  spread(image_id, prob) %>%
  left_join(jing_loglik_degen %>% 
              mutate(prior = "jing") %>%
              rbind_list(marginal_loglik_degen_1 %>%
                           mutate(prior = "marginal_1")) %>%
              rbind_list(marginal_loglik_degen_2 %>%
                           mutate(prior = "marginal_2")) %>%
              rbind_list(marginal_loglik_degen_3 %>%
                           mutate(prior = "marginal_3")) %>%
              rbind_list(shrunk_loglik_degen %>%
                           mutate(prior = "shrunk")) %>%
              rbind_list(shrunk_jing_loglik_degen %>%
                           mutate(prior = "shrunk_jing"))) %>%
  left_join(jing_distance_degen %>% 
              mutate(prior = "jing") %>%
              rbind_list(marginal_distance_degen_1 %>%
                           mutate(prior = "marginal_1")) %>%
              rbind_list(marginal_distance_degen_2 %>%
                           mutate(prior = "marginal_2")) %>%
              rbind_list(marginal_distance_degen_3 %>%
                           mutate(prior = "marginal_3")) %>%
              rbind_list(shrunk_distance_degen %>%
                           mutate(prior = "shrunk")) %>%
              rbind_list(shrunk_jing_distance_degen %>%
                           mutate(prior = "shrunk_jing"))
  ) -> all_statistics_degen


all_statistics %>%
  gather(statistic, value, -iter, -method, -prior) %>%
  filter(prior != "shrunk_jing") %>%
  filter(grepl("image", statistic)) %>%
  separate(statistic, into = c("junk", "statistic")) %>%
  mutate(statistic = as.numeric(statistic)) %>%
  ggplot() +
  geom_line(aes(iter, value, colour = prior)) +
  geom_abline(aes(slope = 0, intercept = prob), data = distn_good %>% rename(statistic = image_id)) +
  facet_grid(method~statistic) +
  theme_bw(base_family = "serif") +
  ylim(c(0,1))

all_statistics_degen %>%
  gather(statistic, value, -iter, -prior) %>%
  filter(prior != "shrunk_jing") %>%
  mutate(marginal = factor(grepl("marginal", prior), labels = c("Latents", "Marginalized"))) %>%
  filter(grepl("image", statistic)) %>%
  separate(statistic, into = c("junk", "statistic")) %>%
  mutate(statistic = as.numeric(statistic)) %>%
  ggplot() +
  geom_line(aes(iter, value, colour = prior)) +
  geom_abline(aes(slope = 0, intercept = prob), data = distn_degen %>% rename(statistic = image_id)) +
  facet_grid(marginal~statistic) +
  theme_bw(base_family = "serif") +
  ylim(c(0,1))

all_statistics %>%
  gather(statistic, value, -iter, -method, -prior) %>%
  #filter(prior != "shrunk_jing") %>%
  filter(grepl("loglik", statistic)) %>%
  ggplot() +
  geom_line(aes(iter, value, colour = prior)) +
  theme_bw(base_family = "serif") +
  facet_grid(method~statistic)

all_statistics_degen %>%
  gather(statistic, value, -iter, -prior) %>%
  filter(prior != "marginal_3") %>%
  filter(grepl("log", statistic)) %>%
  ggplot() +
  geom_line(aes(iter, value, colour = prior)) +
  theme_bw(base_family = "serif") +
  facet_grid(.~statistic)

all_statistics %>%
  gather(statistic, value, -iter, -method, -prior) %>%
  #filter(prior != "shrunk_jing") %>%
  filter(grepl("distance", statistic)) %>%
  ggplot() +
  geom_line(aes(iter, value, colour = prior)) +
  theme_bw(base_family = "serif") +
  facet_grid(method~statistic)

all_statistics_degen %>%
  gather(statistic, value, -iter, -prior) %>%
  filter(prior != "marginal_3") %>%
  filter(grepl("distance", statistic)) %>%
  ggplot() +
  geom_line(aes(iter, value, colour = prior)) +
  theme_bw(base_family = "serif") +
  facet_grid(.~statistic)

all_statistics %>%
  group_by(method, prior) %>%
  do(cor = cor(data.matrix(data.frame(.) %>% select(-method, -prior, -iter)))) -> correlations

correlations %>% 
  rowwise() %>% 
  do(data.frame(data.frame(.$cor, method = .$method, prior = .$prior, var1 = rownames(.$cor)))) %>%
  gather(var2, cor, -var1, -method, -prior) %>%
  mutate(var2 = as.character(var2)) %>% 
  filter(prior != "shrunk_jing") %>%
  ggplot() +
  geom_point(aes(var1, var2, colour = cor, size = abs(cor))) +
  scale_colour_gradient2(limits=c(-1,1)) +
  coord_fixed() +
  theme_bw(base_family = "serif") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  facet_grid(prior ~ method)

all_statistics %>% 
  gather(statistic, value, -method, -prior, -iter) %>%
  filter(grepl("distance", statistic)) %>%
  filter(prior != "shrunk_jing") %>%
  ggplot() +
  geom_boxplot(aes(statistic, value, colour = prior)) +
  facet_wrap(~method)

all_statistics %>% 
  gather(statistic, value, -method, -prior, -iter) %>%
  filter(grepl("log", statistic)) %>%
  filter(prior != "shrunk_jing") %>%
  ggplot() +
  geom_boxplot(aes(statistic, value, colour = prior)) +
  facet_wrap(~method)

all_statistics %>%
  group_by(method, prior) %>%
  summarise_each(funs(mean), -iter) %>%
  gather(metric, mean_value, -method, -prior) %>%
  mutate(combined_id = paste(method, prior, sep = "-")) %>%
  select(-method, -prior) %>%
  spread(combined_id, "mean_value") %>% data.frame



