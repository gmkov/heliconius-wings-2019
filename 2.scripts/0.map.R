##### Altitude and life-history shape the evolution of Heliconius wings ####
#### Gabriela Montejo-Kovacevich 2019 ####
####### packages #######
rm(list=ls())
dev.off()
library(maps)
library(mapdata)
library(RColorBrewer)
library(sp)
library(scatterplot3d)
library(dplyr)
library(sf)
library(cowplot)
library(googleway)
library(ggspatial)
library(ggrepel)
library(ggplot2)
library(googleway)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library("ggplot2")
theme_set(theme_bw())
library("sf")

######## data ######
setwd("../heliconius-wings-2019/")
master <- read.csv("1.data/analyses/master.analyses.csv")

###### localities #####
loc <- summarise(group_by(master, locality),
                 n=n(),
                 lat=mean(latitude),
                 lon=mean(longitude),
                 alt=mean(altitude))

###### map with points with altitude gradient ######
world <- ne_countries(scale = "medium", returnclass = "sf")
theme_set(theme_bw())
ggplot(data = world) +
  geom_sf(color="grey") + ylab("Latitude") + xlab("Longitude")+
  geom_jitter(position = "jitter",data=loc, aes(x=lon, y=lat, color=alt), size=2.5, inherit.aes = FALSE, alpha=.8)+
  coord_sf( xlim=c(-88,-50), ylim=c(-10,15), expand = FALSE)+
  theme(panel.grid.major = element_line(color = "white", linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "white"))+
  annotation_scale(location = "bl", width_hint = 0.08, line_width = 1)+
  scale_color_viridis(name="Altitude\n(m.a.s.l)")+
  scale_y_continuous(breaks=c(-10,0,10))+
  scale_x_continuous(breaks=c(-85, -70,  -50))

# fig 1 main text
#ggsave("3.plots/map_dig1.png", width = 6, height = 4, dpi = 300)


