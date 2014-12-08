library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

calc_hull <- function(H, V, type="binary") {
  #only possible values of t are 0 or 1
  if(type == "binary") t <- data.frame(c(0,1))[,rep(1, H + V)]
  else if(type == "negative") t <- data.frame(c(-1,1))[,rep(1, H + V)]
  else if(type == "experiment") t <- data.frame(c(1,3))[,rep(1, H + V)]
  else stop("Type must be 'binary', 'negative', or 'experiment'")
  
  names(t) <- c(paste0("h", 1:H), paste0("v", 1:V))  
  
  #all possibilities of statistics in (H + V)-space
  t.grid <- expand.grid(t)

  for(i in 1:H){
    for(j in (H+1):(H+V)){
      t.grid <- cbind(t.grid, t.grid[,i]*t.grid[,j])
      names(t.grid)[ncol(t.grid)] <- paste0("theta", i, j - H)
    }
  }
  
  #calculate convex hull
  library(geometry)
  C <- convhulln(t.grid, options="FA")
  return(list(possible_t=t.grid, c_hull=C))
}

test_cases <- expand.grid(data.frame(1:4)[,rep(1,2)]) %>% 
  rename(H = X1.4, V = X1.4.1) %>%
  filter(H <= V) %>% #remove those cases with less visibles than hiddens. Is this necessary?
  mutate(n_param = H*V + H + V) %>%
  mutate(max_facets = (2^(H+V))^(floor(n_param/2))) %>%
  filter(n_param <= 11) #calc_hull can't handle any higher dimensions currently

bin_res <- dlply(test_cases, .(n_param), function(x) calc_hull(x$V, x$H, "binary"), .progress="text")
neg_res <- dlply(test_cases, .(n_param), function(x) calc_hull(x$V, x$H, "negative"), .progress="text")

volume_plot <- ldply(bin_res, function(x) x$c_hull$vol) %>% 
  mutate(frac_vol = V1/(1^n_param)) %>%
  inner_join(ldply(neg_res, function(x) x$c_hull$vol) %>% 
               mutate(frac_vol = V1/(2^n_param)),
             by="n_param") %>%
  rename(vol.bin = V1.x, vol.neg = V1.y, frac_vol.bin = frac_vol.x, frac_vol.neg = frac_vol.y) %>%
  gather(vars, value, -n_param) %>%
  separate(vars, c("type", "encoding"), "\\.") %>%
  spread(type, value) %>%
  ggplot() +
  geom_point(aes(x=n_param, y=frac_vol, colour=encoding)) +
  geom_line(aes(x=n_param, y=frac_vol, colour=encoding, group=encoding)) +
  theme_bw() +
  ylab("Fraction of unrestricted volume") +
  xlab("Number of parameters") +
  scale_colour_discrete("Encoding", labels=c("Binary (1,0)", "Negative (1,-1)"))
  
  


