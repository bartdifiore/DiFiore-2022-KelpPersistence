#########################################################
## Pixel scale kelp canopy metrics
#########################################################

source("code/1_setup.R")

df <- read_csv(("data/Large_files/satellite_kelp/kelp_months.csv"), col_types = cols(.default = col_character()))

patches <- read.table("data/spatial/csvs/So_Cal_kelp_patches.csv", header = T, sep = ",")

df <- df %>% mutate(Lat = as.numeric(Lat), 
              Lon = as.numeric(Lon)) %>%
  left_join(patches)
#NOTE: THIS INTRODUCES 6 NEW ROWS TO THE DATA FRAME AND I CANNOT DETERMINE WHY!!!! NEED TO INSPECT BEFORE FINAL PUBLICATION

df[is.na(df)] <- NA
df$patch <- as.factor(df$Patch)


coords <- readOGR("data/spatial/shapefiles/sampledsites.shp")

aoi <- coords@data %>% filter(site != "AHND") %>% dplyr::select(patch) %>% distinct()
aoi <- as.vector(aoi$patch) #build vector of sites where we collected data


df <- df %>%
  filter(patch %in% aoi) # subset to the patches that we are interested in! 

#########################################
## Buffer the transects
#########################################

dat <- df
coordinates(dat) <- ~ Lon + Lat

proj4string(dat) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"


new.projection <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")

dat <- spTransform(dat, new.projection)

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


# Extract the time series for each point

overlay <- over(sq, dat, returnList = T)

over.df <- do.call(rbind, overlay)

over.df$id <- rownames(over.df)
  
ts <- over.df %>% separate(id, into = c("id", "junk"), sep = "[.]") %>%
  dplyr::select(-junk, -patch) %>%
  mutate(id2 = id) %>%
  separate(id2, 
           into = c("site", "transect"), 
           sep = "(?<=[A-Za-z])(?=[0-9])"
  )
rownames(ts) <- NULL

write.csv(ts, "data/spatial/csvs/timeseries-bysite.csv", row.names = F)

