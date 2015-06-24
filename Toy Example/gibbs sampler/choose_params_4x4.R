#libraries -----------------------
library(dplyr)
library(tidyr)
library(ggplot2)

#function -----------------------
#this function takes the results from the simulation study and rearranges for ease of use
rearrange_data <- function(res) {
  
  data <- data.frame()
  for(i in 1:nrow(res)) {
    tmp <- res$samp[[i]]
    
    H <- res[i,]$H
    V <- res[i,]$V
    
    if(H == 4 & V == 4) { #only looking at 4x4 model
      tmp %>% 
        data.frame() %>% 
        rowwise() %>% 
        mutate_(ss_ratio = paste0("(", paste(paste0(names(.)[(H + V + 1):(H + V + H*V)], "^2"), collapse = " + "), ")/sum(", paste(paste0(names(.)[1:(H+V)],"^2"), collapse = " + "), ")")) %>%    
        ungroup() -> ratio
      
      data <- rbind(data, inner_join(ratio, res$outside[[i]] %>% ungroup()) %>%
        mutate(H = H, V = V, n_param = H + V + H*V))
    }
  }
  return(data)
}

#do it -------------------------
#uses the function to rearrange the data
n <- 100
r_center <- seq(0.2, 3, by = 0.2)
r_exp <- .05

plot.dat <- list()
for(i in seq(0, 1, by = 0.1)) {
  load(paste0("../expected value/apps/results_", i, ".RData")) 
  plot_dat <- lapply(res, rearrange_data)
  names(plot_dat) <- gsub("[.]", "_", r_center)
  
  do.call(rbind, plot_dat) %>% 
    mutate(name = rownames(do.call(rbind, plot_dat))) %>%
    separate(name, into = c("multiplier", "num"), "[.]") %>% 
    mutate(multiplier = as.numeric(gsub("_", ".", multiplier))) %>%
    mutate(n_param_f = factor(n_param)) %>%
    select(-starts_with("exp")) -> plot_dat
  
  plot.dat[[i*10 + 1]] <- plot_dat
}

names(plot.dat) <- gsub("[.]", "_", seq(0, 1, by = 0.1))
do.call(rbind, plot.dat) %>%
  mutate(name2 = rownames(do.call(rbind, plot.dat))) %>%
  separate(name2, into = c("power", "num2"), "[.]") %>% 
  mutate(power = as.numeric(gsub("_", ".", power))) -> plot.dat

#sample --------------------------------
#samples the parameters to get one set from the "degenerate" and one from "non-degenerate"

set.seed(1000) #reproducible

plot.dat %>%
  group_by(near_hull) %>%
  sample_n(1) -> sample.params

save(sample.params, file = "written/sample.params.Rdata")

