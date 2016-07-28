library(ElemStatLearn)
library(ggplot2)

image <- cbind(expand.grid(x = 1:16, y = -1:-16), zip = zip.test[50,-1])
ggplot(image) +
  geom_tile(aes(x = x, y = y, fill = zip), colour = "#0066cc") +
  scale_fill_gradient(low = "white", high = "black") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        aspect.ratio = 1) 
