library(dplyr)
library(tidyr)

load("written/sample_images.Rdata")

#data and params ------------------------
H <- 4
V <- 4

load("written/fitted_models_trunc_marginal_full.Rdata")
marginal <- models_good

load("written/fitted_models_trunc_full.Rdata")
trunc <- models_good

rm(models_good)
rm(models_bad)

source("functs_sample.R")
load("written/params_theta.Rdata")

params_good <- list(main_hidden = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("h"), -H) %>% data.matrix(),
                    main_visible = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("v"), -V) %>% data.matrix(),
                    interaction = sample.params %>% ungroup() %>% filter(!near_hull) %>% select(starts_with("theta")) %>% data.matrix() %>% matrix(4))
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

marginal_sample <- reshape_sample_distn(marginal)
trunc_sample <- reshape_sample_distn(trunc)



marginal_sample %>% 
  group_by(v1, v2, v3, v4, image_id) %>% 
  do(data.frame(acf = acf(.$prob, plot=FALSE)$acf,
     lag = acf(.$prob, plot=FALSE)$lag)) -> marginal_acfs

trunc_sample %>% 
  group_by(v1, v2, v3, v4, image_id) %>% 
  do(data.frame(acf = acf(.$prob, plot=FALSE)$acf,
                lag = acf(.$prob, plot=FALSE)$lag)) -> trunc_acfs

marginal_acfs %>% ungroup() %>% select(-image_id) %>% rename(marginal = acf) %>%
  left_join(trunc_acfs %>% rename(truncated = acf)) %>%
  ggplot() +
  geom_segment(aes(x=lag, xend=lag, y=0, yend=truncated), colour = "black") +
  geom_segment(aes(x=lag, xend=lag, y=0, yend=marginal), colour = "red") +
  # geom_hline(aes(yintercept = qnorm((1 + .95)/2)/sqrt(2000)), colour="red", lty=2) +
  facet_wrap(~image_id) +
  xlab("Lag") +
  ylab("ACF")




  


