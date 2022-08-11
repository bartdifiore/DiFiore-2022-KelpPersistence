library(ggplot2)
library(ggspatial)
library(tidyverse)
library(ggsn)
library(sf)



#-----------------------------------------------------------------
## Make the map
#-----------------------------------------------------------------

#us <- raster::getData("GADM", country = c("United States"), level = 1)

us <- sf::st_read("data/spatial/shapefiles/gadm36_USA_1.shp")

clipper_small <- st_polygon(list(rbind(c(-120.7, 34.65),
                                       c(-119.25, 34.65), 
                                       c(-119.25, 33.8), 
                                       c(-120.7, 33.8), 
                                       c(-120.7, 34.65)))) %>% st_sfc() %>% st_set_crs(4326)


shore_small <- us %>% sf::st_as_sf() %>% sf::st_intersection(clipper_small) %>% sf::st_transform(4326) %>% sf::st_union()

depth <- marmap::read.bathy("data/spatial/sbc.xyz")
coords <- st_coordinates(clipper_small)
depth <- marmap::subsetBathy(depth, x = coords[,1], y = coords[,2], locator = F)
depth.df <- marmap::fortify.bathy(depth)


blues <- colorRampPalette(colors = c("#94AAC3", "#F9FAFB")) #Low #94AAC3, high #F9FAFB
browns <- colorRampPalette(colors = c("#ACD0A5", "#C3A76B"))

persistence <- read.csv("data/intermediary/persistence_metrics.csv") %>%
  group_by(site) %>%
  summarize(perturbations = mean(perturbations))

sites <- st_read("data/spatial/shapefiles/sampledsites.shp") %>%
  st_set_crs("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs ") %>% st_transform(4326) %>%
  group_by(site) %>%
  summarise(st_union(geometry)) %>%
  st_centroid() %>%
  left_join(persistence) %>%
  filter(site != "AHND")

ggplot()+
  geom_sf(data = sites, aes(color = site))

# map no legend
zoom_map <- ggplot()+
  geom_tile(data = filter(depth.df, z < 0), aes(x = x, y = y, fill = z), show.legend = F)+
  scale_fill_gradientn(colours = blues(10))+
  geom_contour(data = filter(depth.df, z < 10), aes(x = x, y = y, z = z), color = "black", binwidth = 100, alpha = 0.25)+
  geom_sf(data = shore_small, fill = "#596778", lwd = 0.01)+
  geom_sf(data = sites, aes(color = perturbations, size = perturbations), alpha = 0.75)+
  scale_x_continuous(breaks = -1*c(120.5, 119.5))+
  scale_y_continuous(breaks = c(33.8, 34.2, 34.6))+
  scale_color_gradient(low = "#fcc5c5", high = "#960f0f", limits = c(1, 7), breaks=seq(1, 7, by=1))+
  scale_size_continuous(limits=c(1,7), breaks=seq(1, 7, by=1), range = c(1.5, 12))+
  labs(x = "", y = "", size = "", color = "")+
  coord_sf(xlim = c(-120.7, -119.25), ylim = c(34, 34.65), expand = F)+
  annotation_scale(location = "tl", style = "ticks",  pad_x = unit(1, "cm"), pad_y = unit(0.25, "cm"), line_width = 2, text_cex = 1)+
  annotation_north_arrow(location = "tr", style = north_arrow_nautical, height = unit(0.75, "cm"), width = unit(0.75, "cm"), pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"))+
  theme(legend.key=element_blank(),legend.background = element_blank(), axis.text = element_text(size = 14), legend.position = "bottom")+
  guides(color= guide_legend(nrow = 1), size=guide_legend(nrow = 1))

ggsave("figures/map.png", device = "png", width = 9, height = 5.5)
