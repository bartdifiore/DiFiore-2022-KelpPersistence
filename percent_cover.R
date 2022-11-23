#----------------------------
## Get data 
#----------------------------

library(ncdf4)
library(raster)
library(rgeos) 
library(rgdal)
library(sp)

nc_data <- nc_open('~/Downloads/LandsatKelpBiomass_2022_Q3_withmetadata.nc')
print(nc_data)

lon <- ncvar_get(nc_data, "longitude")
lat <- ncvar_get(nc_data, "latitude")
year <- ncvar_get(nc_data, "year")
quarter <- ncvar_get(nc_data, "quarter")

cover.array  <- ncvar_get(nc_data, "area")
fillvalue <- ncatt_get(nc_data, "area", "_FillValue")
fillvalue


cover.array[cover.array == fillvalue$value] <- NA


nc_close(nc_data)

#########################################################
## Pixel scale kelp canopy metrics
#########################################################
library(tidyverse)

df <- data.frame(Lon = lon, Lat = lat, cover.array)

col_names <- paste("Y", year, quarter, sep = "_")

names(df) <- c("Lon", "Lat", col_names)

df <- df %>%
  select(Lon, Lat, Y_2008_1:Y_2018_4) %>%
  filter(Lat > 34.0 & Lat < 34.5)

#########################################
## Buffer the transects
#########################################

dat <- df

coordinates(dat) <- ~ Lon + Lat

proj4string(dat) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"


new.projection <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")

dat <- spTransform(dat, new.projection)

boundary <- readOGR("data/spatial/shapefiles/boundary2.shp")

# fix coordinates to be at the center of each of the landsat pixels rather than the upper left corner... 

dat@coords <- cbind(dat@coords[,1]+15, dat@coords[,2]-15)


# Essentially, this approach takes the points where we have community data: 
# 
#         1. Builds a 60m x 60m square buffer around the points.
#         2. Builds a .shp file out of the buffer
#         3. extracts the time series of kelp biomass for each square region
#         4. Examine the alignment of the buffered squares and the "pixel" data as distributed by Tom

####################################

# 1. Get the community data points, and build a 60m x 60m square buffer. 

#       Buffer workflow from https://stackoverflow.com/questions/28665918/create-square-polygons-from-single-centre-coordinates-and-area-in-r and https://www.neonscience.org/field-data-polygons-centroids


# Get the sites where we have community data

coords <- readOGR("data/spatial/shapefiles/sampledsites.shp")
lines <- readOGR("data/spatial/shapefiles/transects-as-lines.shp")
lines <- lines[lines$site != "AHND", ]
coords <- coords[coords$site != "AHND", ]

ww <- 30

# sq <- gBuffer(coords, byid = T, id = coords$id, width = ww, quadsegs = 1, capStyle = "SQUARE")
sq <- gBuffer(lines, byid = T, id = lines$id, width = ww, quadsegs = 1, capStyle = "SQUARE")
px <- gBuffer(dat, byid = T, id = rownames(dat), width = 15, quadsegs = 1, capStyle = "SQUARE")

# Check to see if it worked
plot(sq[sq$id == "MOHK1", ])
plot(lines[lines$id == "MOHK1", ], add = T, col = "red")
plot(dat, add = T)
plot(px, add = T)


# Check to see if it worked
plot(sq[sq$id == "AQUE1" | sq$id == "AQUE2", ])
plot(lines[sq$id == "AQUE1" | sq$id == "AQUE2", ], add = T, col = "red")
#plot(patches, add = T)
plot(dat, add = T)
plot(px, add = T)

# Check to see if it worked
plot(sq[sq$id == "ABUR1" | sq$id == "ABUR2", ])
plot(lines[sq$id == "ABUR1" | sq$id == "ABUR2", ], add = T, col = "red")
#plot(patches, add = T)
plot(dat, add = T)
plot(px, add = T)

writeOGR(sq, dsn = "data/spatial/shapefiles", layer = "squares", driver = "ESRI Shapefile", overwrite_layer = T)

writeOGR(px, dsn = "data/spatial/shapefiles", layer = "pixels", driver = "ESRI Shapefile", overwrite_layer = T)


# Extract the time series for each point

overlay <- over(sq, dat, returnList = T)

over.df <- do.call(rbind, overlay)

over.df$id <- rownames(over.df)

ts <- over.df %>% separate(id, into = c("id", "junk"), sep = "[.]") %>%
  dplyr::select(-junk) %>%
  mutate(id2 = id) %>%
  separate(id2, 
           into = c("site", "transect"), 
           sep = "(?<=[A-Za-z])(?=[0-9])"
  )
rownames(ts) <- NULL

write.csv(ts, "data/spatial/csvs/percover_timeseries-bysite.csv", row.names = F)









