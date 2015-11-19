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

jing_sample_good <- reshape_sample_distn(jing_good)
jing_sample_bad <- reshape_sample_distn(jing_bad)
jing_sample_degen <- reshape_sample_distn(jing_degen)

#computing loglikelihoods, ratios ---------------------
reshape_loglik <- function(model, visibles) {
  if("theta" %in% names(model)) {
    theta <- model$theta
  } else {
    interaction <- t(apply(model$params$interaction, 1, t))
    theta <- cbind(model$params$main_visible, model$params$main_hidden, interaction)
  }
  hiddens <- model$hiddens
  N <- visibles %>% nrow()
  H <- dim(hiddens)[3]
  V <- ncol(visibles)
  
  loglik <- do.call(rbind, lapply(seq_along(theta[, 1]), function(m) {
    
    data <- cbind(visibles, hiddens[m, , ])
    for(i in 1:V) {
      for(j in 1:H) {
        data <- cbind(data, visibles[i, ]*hiddens[m, j, ])
      }
    }
    possibles <- stats(H, V) 
    W <- exp(possibles %*% theta[m, ])
    
    data.frame(iter = m, 
               log_numerator = sum(theta[m, ]*colSums(data)), 
               log_normalizer = N*log(sum(W)),
               reference_loglik = sum(rep(0,length(theta[m, ]))*colSums(data)) - N*log(sum(exp(possibles %*% rep(0,length(theta[m, ]))))))
  }))
  
  loglik %>%
    mutate(loglik = log_numerator - log_normalizer) %>%
    select(-log_numerator) %>%
    mutate(log_ratio = loglik - reference_loglik)
}

shrunk_loglik_good <- reshape_loglik(shrunk_good, flat_images_good$visibles) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_loglik_bad <- reshape_loglik(shrunk_bad, flat_images_good$visibles) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_loglik_degen <- reshape_loglik(shrunk_degen, flat_images_degen$visibles) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 

jing_loglik_good <- reshape_loglik(jing_good, flat_images_good$visibles)
jing_loglik_bad <- reshape_loglik(jing_bad, flat_images_good$visibles)
jing_loglik_degen <- reshape_loglik(jing_degen, flat_images_degen$visibles)

reshape_distance <- function(model) {
  if("theta" %in% names(model)) {
    theta <- model$theta
  } else {
    interaction <- t(apply(model$params$interaction, 1, t))
    theta <- cbind(model$params$main_visible, model$params$main_hidden, interaction)
  }
  
  theta %>%
    data.frame() %>%
    mutate_(main_distance = paste0("sqrt(", paste(paste0("X", 1:(V + H), "^2"), collapse = " + "), ")"),
            interaction_distance = paste0("sqrt(", paste(paste0("X", (V + H + 1):(V + H + V*H), "^2"), collapse = " + "), ")"),
            combined_distance = paste0("sqrt(", paste(paste0("X", 1:(V + H + H*V), "^2"), collapse = " + "), ")")) %>%
    mutate(iter = 1:n()) %>%
    select(iter, contains("distance"))
    
    
  
  
}

shrunk_distance_good <- reshape_distance(shrunk_good) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_distance_bad <- reshape_distance(shrunk_bad) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 
shrunk_distance_degen <- reshape_distance(shrunk_degen) %>% filter(iter > 500) %>% mutate(iter = 1:n()) 

jing_distance_good <- reshape_distance(jing_good)
jing_distance_bad <- reshape_distance(jing_bad)
jing_distance_degen <- reshape_distance(jing_degen)


#plots -------------------
jing_sample_good %>% rename(jing = prob) %>%
  left_join(shrunk_sample_good %>% ungroup() %>% select(-image_id) %>% rename(shrunk = prob)) %>%
  ungroup() %>%
  mutate(method = "good", image_id = paste0("image_", image_id)) %>% 
  rbind_list(jing_sample_bad %>% rename(jing = prob) %>%
              left_join(shrunk_sample_bad %>% ungroup() %>% select(-image_id) %>% rename(shrunk = prob)) %>%
              ungroup() %>% 
              mutate(method = "bad", image_id = paste0("image_", image_id))) %>%
  select_(.dots = paste0("-v", 1:V)) %>%
  gather(prior, prob, jing, shrunk) %>%
  spread(image_id, prob) %>%
  left_join(jing_loglik_good %>% mutate(method = "good") %>%
              rbind_list(jing_loglik_bad %>% mutate(method = "bad")) %>%
              mutate(prior = "jing") %>%
              rbind_list(shrunk_loglik_good %>% mutate(method = "good") %>%
                           rbind_list(shrunk_loglik_bad %>% mutate(method = "bad")) %>%
                           mutate(prior = "shrunk"))) %>%
  left_join(jing_distance_good %>% mutate(method = "good") %>%
              rbind_list(jing_distance_bad %>% mutate(method = "bad")) %>%
              mutate(prior = "jing") %>%
              rbind_list(shrunk_distance_good %>% mutate(method = "good") %>%
                           rbind_list(shrunk_distance_bad %>% mutate(method = "bad")) %>%
                           mutate(prior = "shrunk"))
  ) %>%
  select(-reference_loglik) -> all_statistics

all_statistics %>%
  gather(statistic, value, -iter, -method, -prior) %>%
  filter(grepl("image", statistic)) %>%
  separate(statistic, into = c("junk", "statistic")) %>%
  mutate(statistic = as.numeric(statistic)) %>%
  ggplot() +
  geom_line(aes(iter, value, colour = prior)) +
  geom_abline(aes(slope = 0, intercept = prob), data = distn_good %>% rename(statistic = image_id)) +
  facet_grid(method~statistic) +
  ylim(c(0,1))

all_statistics %>%
  gather(statistic, value, -iter, -method, -prior) %>%
  filter(!grepl("image", statistic)) %>%
  ggplot() +
  geom_line(aes(iter, value, colour = prior)) +
  facet_wrap(method~statistic, scales = "free_y", nrow = 2)

all_statistics %>%
  gather(statistic, value, -iter, -method, -prior) %>%
  filter(method == "bad") %>%
  ggplot() +
  geom_line(aes(iter, value, colour = prior)) +
  facet_wrap(~statistic, scales = "free_y")

all_statistics %>%
  group_by(method, prior) %>%
  do(cor = cor(data.matrix(data.frame(.) %>% select(-method, -prior, -iter)))) -> correlations

i <- 4
correlations$cor[[i]] %>%
  data.frame() %>%
  mutate(var1 = rownames(correlations$cor[[i]])) %>%
  gather(var2, cor, -var1) %>% 
  mutate(var2 = as.character(var2)) %>% 
  ggplot() +
  geom_tile(aes(var1, var2, fill = cor)) +
  scale_fill_gradient2(limits=c(-1,1)) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle(paste0("Method: ", correlations$method[i], ", Prior: ", correlations$prior[i]))



