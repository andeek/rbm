#libraries --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
source("sample_functs.R")

#data and params ------------------------
H <- 4
V <- 4

load("written/sample_images.Rdata")
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

## all models
model_files <- paste0("written/", list.files("written/")[grep("fitted_models_distn_", list.files("written/"))])

load_model <- function(file) {
  load(file)
  sample_good <- reshape_sample_distn(models_good)
  sample_bad <- reshape_sample_distn(models_bad)
  sample_degen <- reshape_sample_distn(models_degen)
  sample_good$const <- sample_bad$const <- sample_degen$const <- as.numeric(gsub(".Rdata", "", strsplit(file, "_")[[1]][4]))
  sample_good$N <- sample_bad$N <- sample_degen$N <- dim(models_good$hiddens)[2]
  sample_good$model <- "good"
  sample_bad$model <- "bad"
  sample_degen$model <- "degen"
  bind_rows(sample_good, sample_bad, sample_degen)
}

data.frame(file = model_files) %>%
  mutate(file = as.character(file)) %>%
  rowwise() %>%
  do(models = load_model(.$file)) -> all_models

rbind_all(all_models$models) -> all_models
  
all_models %>%
  left_join(rbind_all(list(distn_good %>% mutate(model = "good") %>% rename(true_prob = prob), 
                           distn_good %>% mutate(model = "bad") %>% rename(true_prob = prob), 
                           distn_degen %>% mutate(model = "degen") %>% rename(true_prob = prob)))) %>%
  mutate(diff_sq = (true_prob - prob)^2) %>%
  group_by(image_id, const, model) %>%
  summarise(mse = mean(diff_sq)) %>%
  filter(model %in% c("good", "bad")) %>%
  ggplot(aes(const, mse, colour = model)) +
  geom_point() +
  geom_line() +
  facet_wrap(~image_id) +
  xlab(paste0("Constant*", unique(all_models$N)))


all_models %>%
  left_join(rbind_all(list(distn_good %>% mutate(model = "good") %>% rename(true_prob = prob), 
                           distn_good %>% mutate(model = "bad") %>% rename(true_prob = prob), 
                           distn_degen %>% mutate(model = "degen") %>% rename(true_prob = prob)))) %>%
  mutate(diff_sq = (true_prob - prob)^2) %>%
  group_by(image_id, const, model) %>%
  summarise(mse = mean(diff_sq)) %>%
  filter(model %in% c("good", "bad")) %>%
  filter(const < .5) %>%
  ggplot(aes(const, mse, colour = model)) +
  geom_point() +
  geom_line() +
  facet_wrap(~image_id, scales = "free") +
  xlab(paste0("Constant*", unique(all_models$N)))

all_models %>%
  left_join(rbind_all(list(distn_good %>% mutate(model = "good") %>% rename(true_prob = prob), 
                           distn_good %>% mutate(model = "bad") %>% rename(true_prob = prob), 
                           distn_degen %>% mutate(model = "degen") %>% rename(true_prob = prob)))) %>%
  mutate(diff_sq = (true_prob - prob)^2) %>%
  group_by(image_id, const, model) %>%
  summarise(mse = mean(diff_sq)) %>%
  group_by(const, model) %>%
  summarise(total_mse = sum(mse)) %>%
  ggplot(aes(const, total_mse, colour = model)) +
  geom_point() +
  geom_line() +
  xlab(paste0("Constant*", unique(all_models$N)))

all_models %>%
  filter(const >= .24 & const <= .26 & model == "good") %>%
  ggplot() +
  geom_histogram(aes(prob)) +
  geom_segment(aes(x = prob, xend = prob, y = 0, yend = 500), data = distn_good, colour = "blue") +
  facet_grid(const ~ image_id) 

all_models %>%
  filter(const >= .24 & const <= .26 & model == "bad") %>%
  ggplot() +
  geom_histogram(aes(prob)) +
  geom_segment(aes(x = prob, xend = prob, y = 0, yend = 500), data = distn_good, colour = "blue") +
  facet_grid(const ~ image_id) 
  
all_models %>%
  filter(const >= .24 & const <= .26 & model == "degen") %>%
  ggplot() +
  geom_histogram(aes(prob)) +
  geom_segment(aes(x = prob, xend = prob, y = 0, yend = 500), data = distn_degen, colour = "blue") +
  facet_grid(const ~ image_id) 
  



  