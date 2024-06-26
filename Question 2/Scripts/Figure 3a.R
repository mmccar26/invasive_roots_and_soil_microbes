#===load libraries======
library(ggplot2)
library(tidyverse)
library(vegan)
library(plot3Drgl)
library(scatterplot3d)
library(plot3D)

#=====import data======
#relative pathname
plt <- file.path(".", "Data", "dataForOrdination.csv")

##Load data
ord<-read_csv(plt)%>%
  mutate(plot_cover = case_when(plotStatus == "Highly Invaded" ~ "red",
                                plotStatus == "Moderately Invaded" ~ "yellow",
                                plotStatus == "Not Invaded" ~ "#56B4E9"))

mod1<-plot3d(x = ord$PC1, y = ord$PC2, z = ord$PC3, col = ord$plot_cover, 
       xlab = "PC1",
       ylab = "PC2",
       zlab = "PC3",
       type = "s", 
      legend3d("topright"), size = 1)
bg3d(color = "white") 
bbox3d(col = "lightgrey") 
par3d(cex=1.0)


rgl.snapshot('mod1.png', fmt = 'png')
