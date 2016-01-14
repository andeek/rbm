#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("sample_functs.R")

#data and params ------------------------
H <- 4
V <- 4

load("written/sample_images.Rdata")

load("written/fitted_models_adaptive_mh_full_trunc_distn_marginal_1.Rdata")
marginal_bad_1 <- models_bad
marginal_good_1 <- models_good
load("written/fitted_models_adaptive_mh_full_trunc_distn_marginal_2.Rdata")
marginal_bad_2 <- models_bad
marginal_good_2 <- models_good
load("written/fitted_models_adaptive_mh_full_trunc_distn_marginal_3.Rdata")
marginal_bad_3 <- models_bad
marginal_good_3 <- models_good


#rm unneccesary data
rm(models_bad)
rm(models_good)

load("written/sample.params.Rdata")
params <- sample.params %>% 
  mutate(label = ifelse(near_hull, "degen", "good")) %>% 
  ungroup() %>%
  select(label, starts_with("v"), starts_with("h"), starts_with("theta"), -H, -V) %>%
  gather(theta, value, -label)


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

# marginal_loglik_good_1 <- reshape_loglik(marginal_good_1, flat_images_good$visibles)
# marginal_loglik_bad_1 <- reshape_loglik(marginal_bad_1, flat_images_good$visibles)
# 
# marginal_loglik_good_2 <- reshape_loglik(marginal_good_2, flat_images_good$visibles)
# marginal_loglik_bad_2 <- reshape_loglik(marginal_bad_2, flat_images_good$visibles)
# 
# marginal_loglik_good_3 <- reshape_loglik(marginal_good_3, flat_images_good$visibles)
# marginal_loglik_bad_3 <- reshape_loglik(marginal_bad_3, flat_images_good$visibles)
# 
# save(marginal_loglik_good_1, marginal_loglik_bad_1, marginal_loglik_good_2, marginal_loglik_bad_2, marginal_loglik_good_3, marginal_bad_3, file = "written/test_compare_marginal_inits.RData")
load("written/test_compare_marginal_inits.RData")

##plots --------------------------
colnames(marginal_good_3$theta) <- colnames(marginal_bad_3$theta) <- colnames(marginal_good_1$theta) <- colnames(marginal_good_2$theta) <- colnames(marginal_bad_1$theta) <- colnames(marginal_bad_2$theta) <- unique(params$theta)

marginal_good_1$theta %>%
  data.frame %>%
  mutate(iter = 1:n()) %>%
  gather(theta, value, -iter) %>% 
  mutate(inits = "1") %>% 
  rbind_list(marginal_good_2$theta %>%
               data.frame %>%
               mutate(iter = 1:n()) %>%
               gather(theta, value, -iter) %>% 
               mutate(inits = "2")) %>%
  rbind_list(marginal_good_3$theta %>%
               data.frame %>%
               mutate(iter = 1:n()) %>%
               gather(theta, value, -iter) %>% 
               mutate(inits = "3")) %>%
  ggplot() +
  geom_line(aes(iter, value, colour = inits)) +
  geom_abline(aes(slope = 0, intercept = value), data = params %>% filter(label == "good")) +
  facet_wrap(~theta)


marginal_loglik_bad_1 %>% mutate(inits = "1", prior = "bad") %>%
  rbind_list(marginal_loglik_bad_2 %>% mutate(inits = "2", prior = "bad")) %>%
  rbind_list(marginal_loglik_bad_3 %>% mutate(inits = "3", prior = "bad")) %>%
  rbind_list(marginal_loglik_good_1 %>% mutate(inits = "1", prior = "good")) %>%
  rbind_list(marginal_loglik_good_2 %>% mutate(inits = "2", prior = "good")) %>%
  rbind_list(marginal_loglik_good_3 %>% mutate(inits = "3", prior = "good")) %>%
  gather(metric, value, -iter, -prior, -inits) %>%
  ggplot() +
  geom_line(aes(iter, value, colour = inits)) +
  facet_grid(metric ~ prior)
  
  