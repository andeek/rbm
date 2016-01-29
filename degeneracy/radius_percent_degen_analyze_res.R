
## analyze results -------------------------------------
## plot function ----------------------------------
plot_data <- function(res) {
  
  plot.data <- data.frame()
  for(i in 1:nrow(res)) {
    tmp <- res$samp[[i]]
    H <- res[i,]$H
    V <- res[i,]$V
    
    tmp %>% 
      data.frame() %>% 
      rowwise() %>% 
      mutate_(ss_ratio = paste0("(", paste(paste0(names(.)[(H + V + 1):(H + V + H*V)], "^2"), collapse = " + "), ")/sum(", paste(paste0(names(.)[1:(H+V)],"^2"), collapse = " + "), ")")) %>%    
      ungroup() -> ratio
    
    inner_join(ratio, res$outside[[i]] %>% ungroup()) %>%
      select(ss_ratio, near_hull) %>%
      mutate(H = H, V = V, n_param = H + V + H*V) %>%
      rbind(plot.data) -> plot.data
  }
  return(plot.data)
}

make_plot <- function(plot.data) { 
  plot.data %>%
    ggplot() + 
    geom_boxplot(aes(as.factor(n_param), ss_ratio, colour = near_hull)) +
    xlab("Number of params") +
    ylab("Ratio of sum of squares of cross terms to main effects") +
    scale_colour_discrete("Within .05 of hull") -> plot1
  
  plot.data %>% 
    ggplot() + 
    geom_bar(aes(as.factor(n_param), y = ..count.., fill = near_hull), position = "fill") +
    xlab("Number of params") +
    scale_fill_discrete("Within .05 of hull") -> plot2
  
  return(list(p1 = plot1, p2 = plot2))
}

plot.dat <- lapply(res, plot_data)
plots <- lapply(plot.dat, make_plot)


library(gridExtra)
for (i in 1:length(plots)) {
  grid.arrange(plots[[i]][[1]], plots[[i]][[2]], main = paste0("Radius multiplier: ", r_center[i]))
}


names(plot.dat) <- gsub("[.]", "_", r_center)
do.call(rbind, plot.dat) %>%
  ungroup() %>%
  mutate(name = rownames(do.call(rbind, plot.dat))) %>%
  separate(name, into = c("multiplier", "num"), "[.]") %>% 
  mutate(multiplier = as.numeric(gsub("_", ".", multiplier))) %>%
  group_by(n_param, multiplier) %>%
  summarise(frac_degen = sum(near_hull)/n()) %>%
  ggplot() +
  geom_point(aes(multiplier, frac_degen, colour = factor(n_param))) +
  geom_line(aes(multiplier, frac_degen, group = factor(n_param), colour = factor(n_param))) +
  #geom_smooth(aes(multiplier, frac_degen, group = factor(n_param), colour = factor(n_param), fill = factor(n_param)), alpha = 0) +
  xlab("Radius Multiplier") + ylab("Fraction degenerate") +
  scale_colour_discrete("Num. of Parameters") +
  scale_fill_discrete("Num. of Parameters") +
  theme_bw()


