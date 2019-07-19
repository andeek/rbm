# libraries -----------------------
library(ggplot2) # plots
library(dplyr) # manipulate data frames
library(tidyr) # tidy data
library(purrr) # map functions to lists
library(geometry) # for hull functions

# source additional functions
source("functions.R")

# plot theme
theme_set(theme_bw(base_family="serif"))

# more functions ----------------------
# sample uniformly on a sphere
sample_sphere_unif <- function(space, rep, center, r) {
  vec <- matrix(rnorm(space*rep), nrow = rep)
  r*vec/t(matrix(rep(sqrt(rowSums(vec^2)), each = space), nrow = space))
}

# sample thetas --------------
n <- 100
V <- 9
H <- 5

expand.grid(H = H, V = V, r1 = seq(0.001, 3, length.out = 20), r2 = seq(0.001, 3, length.out = 20)) %>%
  mutate(r1 = r1*(H + V), r2 = r2*(H*V)) %>%
  mutate(N = H + V) %>%
  group_by(H, V, r1, r2, N) %>%
  do(samp = cbind(sample_sphere_unif(.$N, n, 0, .$r1), sample_sphere_unif(.$H*.$V, n, 0, .$r2))) -> grid_sample

expand.grid(H = H, V = V, r1 = seq(0.001, 3, length.out = 20), r2 = seq(0.001, 3, length.out = 20)) %>%
  mutate(r1 = r1*(H + V), r2 = r2*(H*V)) %>%
  mutate(N = H + V) %>%
  group_by(H, V, r1, r2, N) %>%
  do(stat = stats(.$H, .$V, "negative")) -> grid_sample_stat

grid_sample <- inner_join(grid_sample, grid_sample_stat)

grid_sample %>%
  group_by(H, V, N, r1, r2) %>%
  do(g_theta = t(rbind(t(.$samp[[1]]), expected_value(theta = t(.$samp[[1]]), stats = data.matrix(.$stat[[1]]))))) -> grid

res <- inner_join(grid, grid_sample)

# reshape data functions ----------------
elpr <- function(theta, stats) {
  exp(tcrossprod(stats, theta)) %>% data.frame() %>%
    cbind(stats) %>%
    group_by_(.dots = colnames(stats)[grepl("v", colnames(stats)) & !grepl("theta", colnames(stats))]) %>%
    summarise_each(funs(sum), contains("X")) %>%
    ungroup() -> marginalized
  
  marginalized %>%
    summarise_each(funs(max), contains("X")) %>%
    select(contains("X")) %>%
    data.matrix() -> max_marg
  
  marginalized %>%
    summarise_each(funs(min), contains("X")) %>%
    select(contains("X")) %>%
    data.matrix() -> min_marg
  
  t(log(max_marg/min_marg))
}

# instability
res %>%
  group_by(H, V, N, r1, r2) %>%
  do(data.frame(elpr = elpr(.$samp[[1]], .$stat[[1]]))) %>%
  ungroup() %>%
  group_by(H, V, N, r1, r2) %>% 
  summarise(mean_elpr = mean(elpr)) %>%
  mutate(scaled_mean_elpr = mean_elpr/V) %>%
  ungroup() -> elpr_summary

res %>%
  group_by(H, V, N, r1, r2) %>%
  do(data.frame(max_prob = apply(.$samp[[1]], 1, function(row) {
    params <- list(main_visible = row[1:V] %>% data.matrix(),
                   main_hidden = row[(V + 1):(V + H)] %>% data.matrix(),
                   interaction = row[(V + H + 1):length(row)] %>% data.matrix() %>% matrix(H))
    
    distn <- visible_distn(params) %>% mutate(out_vals = pmap(lst_(paste0("v", 1:V)), c)) %>% select(-(starts_with("v")))
    
    compare <- expand.grid(out1 = distn$image_id, out2 = distn$image_id) %>%
      left_join(distn, by = c("out1" = "image_id")) %>%
      left_join(distn, by = c("out2" = "image_id")) %>%
      mutate(test = map2_dbl(out_vals.x, out_vals.y, function(a, b) sum(a != b))) %>%
      filter(test == 1) %>%
      mutate(delta = log(prob.x/prob.y))
    
    max(compare$delta)
    })
    )) %>%
  ungroup() %>%
  group_by(H, V, N, r1, r2) %>% 
  summarise(mean_max_prob = mean(max_prob)) %>%
  ungroup() -> mean_max_prob_diff_summary


## ----rbm-plots----
elpr_summary %>%
  mutate(r1 = round(r1/(H + V), 8), r2 = round(r2/(H * V), 8)) %>% 
  ggplot() +
  geom_tile(aes(x = r1, y = r2, fill = scaled_mean_elpr)) +
  geom_contour(aes(x = r1, y = r2, z = scaled_mean_elpr), colour = "black", bins = 8) +
  geom_abline(aes(intercept = 0, slope = 1), alpha = .5, lty = 2) +
  scale_fill_gradient(expression(paste("Mean ", frac(ELPR(theta), N[italic(V)]))), low = "white", high = "grey20") +
  xlab(expression(paste(paste("||", theta[main], "||", "/(", N[italic(H)]," + ", N[italic(V)], ")")))) +
  ylab(expression(paste(paste("||", theta[interaction], "||", "/(", N[italic(H)], "*", N[italic(V)], ")")))) +
  theme(aspect.ratio = 1, legend.position = "bottom")

# fix labels
mean_max_prob_diff_summary %>%
  mutate(r1 = round(r1/(H + V), 8), r2 = round(r2/(H * V), 8)) %>% 
  ggplot() +
  geom_tile(aes(x = r1, y = r2, fill = mean_max_prob)) +
  geom_contour(aes(x = r1, y = r2, z = mean_max_prob), colour = "black", bins = 8) +
  geom_abline(aes(intercept = 0, slope = 1), alpha = .5, lty = 2) +
  scale_fill_gradient(expression(paste("Mean ", Delta[N](theta))), low = "white", high = "grey20") +
  xlab(expression(paste(paste("||", theta[main], "||", "/(", N[italic(H)]," + ", N[italic(V)], ")")))) +
  ylab(expression(paste(paste("||", theta[interaction], "||", "/(", N[italic(H)], "*", N[italic(V)], ")")))) +
  theme(aspect.ratio = 1, legend.position = "bottom")

