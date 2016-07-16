#libraries -----------------------
library(dplyr)
library(tidyr)
library(ggplot2)


#look at functional relationship ----------
expand.grid(V = seq(1, 256, 5), H = seq(1, 100, 5)) %>%
  mutate(r_main = 1/(H + V), r_int = 1/(H*V)) -> params

params %>%
  ggplot() +
  geom_line(aes(V, r_main, colour = H, group = H)) 

params %>%
  ggplot() +
  geom_tile(aes(V, H, fill = r_main/r_int)) +
  geom_contour(aes(V, H, z = r_main/r_int), colour = "black") +
  scale_fill_gradient(low = "yellow", high = "red")

expand.grid(V = 4, H = 1:10) %>%
  mutate(r_main = 1/(H + V), r_int = 1/(H*V)) %>%
  gather(which, r, r_main, r_int) %>%
  ggplot() +
  geom_line(aes(H, r, group = which, colour = which))

  