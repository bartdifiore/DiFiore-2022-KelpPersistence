library(rgeos) 
library(rgdal)
library(sp)
library(ggpubr)
library(tidyverse)
library(dplyr)


################################################
## BB Survey's Spatial Data
################################################

# Read in metadata
meta <- read.csv("data/cleaned/metadata.csv") %>%
  select(-notes)

# clean up metadata
coords <- meta %>%
  select(lat_decimal, lon_decimal, site, transect) %>%
  drop_na() %>%
  mutate(lon_decimal = -1*lon_decimal, 
         lat_decimal = lat_decimal)


# convert to a spatial data.frame
coordinates(coords) <- ~lon_decimal + lat_decimal

proj4string(coords) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

new.projection <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")

coords <- spTransform(coords, new.projection)
plot(coords)


# build data.frame to recalculate the OAKS lat x lon as we only had one waypoint
oaks.fix <- data.frame(site = "OAKS", transect = c(1,2,3), lon_decimal = c(758458.9, 758458.9, 758458.9), lat_decimal = c(3817637,3817637+20, 3817637+40)) %>%
  select(lat_decimal, lon_decimal, site, transect)

coordinates(oaks.fix) <- ~lon_decimal + lat_decimal
proj4string(oaks.fix) <- new.projection

# add in new OAKS lat x lon
coords <- coords[coords$site != "OAKS", ] 

coords <- rbind(coords, oaks.fix)

# build data.frame to recalculate the TANK, INNP, and REFU lat x lon as we estimated some transect start positions as we dove multiple transects at a time

fix <- data.frame(site = c("INNP", "TANK", "REFU"), 
                  transect = c(3,2,2), 
                  lon_decimal = c(coords@coords[coords$site == "INNP" & coords$transect == 3][1], coords@coords[coords$site == "TANK" & coords$transect == 1][1]+60, coords@coords[coords$site == "REFU" & coords$transect == 1][1]), 
                  lat_decimal = c(coords@coords[coords$site == "INNP" & coords$transect == 3][2]-20, coords@coords[coords$site == "TANK" & coords$transect == 1][2], coords@coords[coords$site == "REFU" & coords$transect == 1][2]+35)) %>%
                    mutate(id = paste(site, transect, sep = "")) %>%
                    select(lat_decimal, lon_decimal, site, transect, id)

coordinates(fix) <- ~lon_decimal + lat_decimal
proj4string(fix) <- new.projection

coords$id <- paste(coords$site, coords$transect, sep ="")
coords <- coords[coords$id != "INNP3", ]

coords <- rbind(coords, fix)

plot(coords)


###################################################
## SBC LTER Site Spatial Data
###################################################

lter <- read.csv("data/spatial/csvs/lter_waypoints.csv") %>%
  filter(treatment != "annual" & treatment != "continual") %>%
  mutate(id = paste(site, transect, sep = ""))

coordinates(lter) <- ~long + lat

proj4string(lter) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

new.projection <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")

lter <- spTransform(lter, new.projection)
plot(lter)

lter <- lter[, c("site", "transect", "id")]


###################################################
## Combine all sites and extract persistence data
###################################################

coords <- rbind(coords, lter)
plot(coords)

# read in the persistence datafile
dist <- readOGR("data/spatial/shapefiles/localpatchdisturbance.shp")

# extract persistence information for each transect 
temp <- over(coords, dist)
coords@data <- data.frame(site = coords$site, transect = coords$transect, id = coords$id, patch = temp$patch)

#If transects fell outside of patch boundaries assign them to the closest patch.
coords[coords$id == "STEP1"] <- c("STEP", "1", "STEP1", "151")
coords[coords$id == "LEAD3"] <- c("LEAD", "3", "LEAD3", "152")
coords[coords$id == "NAPL8"] <- c("NAPL", "8", "NAPL8", "181")
coords[coords$id == "AQUE3"] <- c("AQUE", "3", "AQUE3", "209")

coords$transect <- as.numeric(coords$transect)

writeOGR(coords, dsn = "data/spatial/shapefiles",  layer = "sampledsites", driver = "ESRI Shapefile", overwrite_layer = T)

dat <- read.csv("data/final_biomass/allsites_combined_biomass_2018.csv") 

dat <- dat %>% left_join(coords@data) 

write.csv(dat, "data/intermediary/species_withpatches.csv", row.names = F)


#########################################################
## Study design summary
#########################################################

coords@data %>% 
  group_by(patch) %>%
  summarize(f = length(transect)) %>%
  View()


fig <- coords@data %>% 
  group_by(patch) %>%
  summarize(f = length(transect)) %>%
  filter(f >= 3) %>%
  ggplot()+
  geom_point(aes(x = patch, y = f))+
  annotate("text", x = 10, y = 5, label = "n = 20 patches", size = 10)+
  coord_flip()+
  labs(x = "Patch ID", y = "Number of transects")+
  theme_pubclean()

ggsave("figures/study_replication.png", fig, device = "png")



########################################################
## Site map
########################################################



# build map for orientation purposes

#get bathy data
library(marmap)
bath <- getNOAA.bathy(lon1 = -120.5, lon2 = -119.5,
                      lat1 = 34.3, lat2 = 34.65, resolution = 1)

#get cal state outline
cal <- readOGR("~/GitHub/Rennick-2019-urchin-grazing/data/spatial/")
plot(cal)


d <- par(las = 1, mgp = c(3, 0.75, 0))
blues <- colorRampPalette(c('#ece7f2','#a6bddb','#2b8cbe'))
plot(bath, image = TRUE, land = TRUE, lwd = 0,  xlim = c(-120.5,-119.5), ylim = c(34.3,34.65), xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex.axis = 1.5, bpal = blues(100)) # bathy layer
plot(cal, col = "#FFEB9B", add =T) # cal state
plot(spTransform(dist[dist$patch %in% temp$patch, ], "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"), add = T, col = "forestgreen")
plot(spTransform(coords, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"), pch = 4, add = T)
par(d)




# zoom to IVEE for example


d <- par(las = 1, mgp = c(3, 0.75, 0))
plot(spTransform(dist[dist$patch %in% temp$patch, ], "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"), col = "forestgreen",  xlim = c(-119.82,-119.77), ylim = c(34.40,34.42))
plot(cal, col = "#FFEB9B", add =T) # cal state
plot(spTransform(coords, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"), pch = 4, add = T, col = "red", cex = 2)
par(d)















                       
     