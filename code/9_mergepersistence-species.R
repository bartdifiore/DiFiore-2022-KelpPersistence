source("code/1_setup.R")
source("code/7_dist_function.R")

library(rgeos) 
library(rgdal)
library(sp)
library(ggpubr)
library(tidyverse)


############################################################
## read in data
############################################################

coords <- readOGR("data/spatial/shapefiles/sampledsites.shp")@data %>%
  dplyr::select(id, patch) %>%
  filter(id != "AHND1" | id != "AHND2") # this provides the patch level designations, which is the scale that the waves were measured at..


ts <- read.csv("data/spatial/csvs/timeseries-bysite.csv") %>%
  mutate(across(Y1984.1:Y2018.12,as.numeric)) %>%
  mutate(across(Y1984.1:Y2018.12, ~na_if(., "NaN"))) %>%
  left_join(coords)

#############################################################

# Estimate quarterly averages
# Quarter 1 is Jan - Mar, (winter)
# Quarter 2 is Apr -Jun, (spring)
# Quarter 3 is July - Sep, (summer)
# Quarter 4 is Oct - Dec. (fall)

formerge <- data.frame(month = 1:12, quarter = rep(c("winter", "spring", "summer", "fall"), each = 3), quart.id = rep(c(0, 25, 5, 75), each = 3))


wide <- ts %>% 
  group_by(id) %>%
  summarize_at(vars(matches("y")), funs(mean), na.rm=T) %>%
  gather(key = time, value = measurement, -id) %>% #this will gather into long form...
  separate(time, into = c("junk", "year", "month"), sep = "Y|[.]") %>%
  filter(year >= 2008) %>%
  mutate(junk = NULL, 
         month = as.integer(month)) %>%
  left_join(formerge) %>%
  mutate(year.quarter = as.numeric(paste(year, quart.id, sep = ".")))%>%
  group_by(year.quarter, id) %>%
  summarize(biomass = mean(measurement, na.rm = T)) %>%
  spread(id, biomass)


long <- ts %>% 
  group_by(id) %>%
  summarize(across(Y1984.1:Y2018.12, ~mean(., na.rm = T))) %>%
  mutate(across(Y1984.1:Y2018.12, ~na_if(., "NaN"))) %>%
  pivot_longer(cols = Y1984.1:Y2018.12, names_to = "time", values_to = "measurement") %>%
  tidyr::separate(time, into = c("junk", "year", "month"), sep = "Y|[.]") %>%
  mutate(junk = NULL, 
         month = as.integer(month)) %>%
  left_join(formerge) %>%
  group_by(year, quarter, id) %>%
  summarize(biomass = mean(measurement, na.rm = T)) %>%
  mutate(biomass = na_if(biomass, "NaN"))


#############################################################
## Calculate summary stats for last 10 years at each id
#############################################################

# We are interesting in understanding how variable kelp is at each site. Therefore, we will calculate mean, median, and CV of kelp biomass at each id. 

sum <- long %>%
  filter(year >= 2008) %>%
  group_by(id) %>%
  summarize(mean.canopy = mean(biomass, na.rm = T), sd.canopy = sd(biomass, na.rm = T), median.canopy = median(biomass, na.rm = T), cv.canopy = sd(biomass, na.rm = T)/mean(biomass, na.rm = T)) %>% 
  arrange(id)


# We are also interesting in understanding how recent kelp biomass incluences community structure. Therefore, this included the mean biomass of kelp over the previous year. 

sum <- long %>%
  filter(year == 2017 & quarter %in% c("summer", "fall") | year == 2018 & quarter %in% c("winter", "spring")) %>%
  group_by(id) %>%
  summarize(mean.canopy_previousyear = mean(biomass, na.rm = T)) %>% 
  arrange(id) %>% 
  right_join(sum)

###########################################################################
## Alternative metrics of variability
###########################################################################
wide.temp <- as.data.frame(wide[,-1])

sum$pv <- apply(wide.temp, 2, PV) # AQUE contains an NaN... need to figure out why
hist(sum$pv)


sum$d <- apply(wide.temp, 2, Dc)
hist(sum$d)



##########################################################################
## Calculate number of extinction events and time since last extinction
##########################################################################

sum$extinctions <- bd.apply(wide[,-1], extinctions, row.ids = wide$year.quarter)


timesincelast <- vector()


for(i in 1:dim(wide.temp)[2]){
  timesincelast[i] <- time.since(col = wide.temp[,i], time.vec = wide$year.quarter)
}

sum$timesincelastextinct <- timesincelast


##############################################################################
## All disturbance events regardless of season
##############################################################################

#number of times kelp has disappeared (< 80% reduction) for more than 2 quarters.
sum$perturbations <- bd.apply(wide[,-1], perturbations, row.ids = wide$year.quarter)


##############################################################################
## Merge into single data file
##############################################################################

persistence <- sum

psych::pairs.panels(persistence %>% select(-id), ellipses = F)

psych::pairs.panels(persistence %>% select(mean.canopy, cv.canopy, perturbations, timesincelastextinct, extinctions), ellipses = F)

psych::pairs.panels(persistence %>% select(mean.canopy, cv.canopy, d, pv), ellipses = F)


###########################################################################
## Merge with community data by patch 
###########################################################################

dat <- read.csv("data/intermediary/species_withpatches.csv") %>%
  mutate(patch = as.factor(patch)) %>%
  filter(!site %in% c("SCDI", "SCTW", "AHND") ) %>%
  mutate(group = case_when(coarse_grouping == "SESSILE INVERT" & !sp_code %in% c("PACA", "CHOV") ~ "endoSI", 
                           coarse_grouping == "GIANT KELP" ~ "kelp",
                           coarse_grouping == "UNDERSTORY ALGAE" ~ "ua")) %>%
  select(c(year:common_name, id, patch, group)) %>%
  filter(group %in% c("endoSI", "kelp", "ua"))

dat <- dat %>% left_join(persistence, by = "id") 

# There are 17 sites in the final dataset sampled in 2018!!!!

write.csv(dat, "data/intermediary/species_withpersistence-pixelscale.csv", row.names = F)

write.csv(wide, "data/intermediary/kelpcanopytimeseries_wide.csv", row.names = F)

###########################################################################
## Clean up the substrate data and merge with the LTER core sites
###########################################################################


bbs <- read.csv("data/intermediary/substrate/BBenthic_substrate.csv") %>%
  dplyr::select(-DATE, -QUAD, -SIDE)


lts <- read.csv("data/intermediary/substrate/Annual_Substrate_All_Years.csv") %>%
  dplyr::select(-DATE, -QUAD, -SIDE) %>%
  filter(YEAR == 2018, SITE != "SCDI" | SITE != "SCTW")

df <- bind_rows(bbs, lts) %>% group_by(SITE, TRANSECT, SUBSTRATE_TYPE) %>%
  summarize(percent_cover = mean(PERCENT_COVER, na.rm = T))

write.csv(df, "data/intermediary/substrate/combined_substrate.csv", row.names = F)

###########################################################################
## Plot kelp biomass dynamics
###########################################################################

forplot <- sum %>%
  separate(id, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  group_by(site) %>%
  summarize(pv = mean(pv))


long %>% 
  mutate(quarter.id = case_when(quarter == "fall" ~ ".75", 
                                quarter == "winter" ~ ".0", 
                                quarter == "spring" ~ ".25", 
                                quarter == "summer" ~ ".5"), 
         year.quarter = as.numeric(paste(year, quarter.id, sep = ""))) %>%
  filter(year >= 2008) %>%
  separate(id, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  left_join(forplot) %>%
  mutate(site = forcats::fct_reorder(site, pv)) %>%
  ggplot(aes(x = year.quarter, y = biomass))+
  geom_line(aes(color = site, group = transect))+
  facet_wrap(~site)+
  theme_classic()
ggsave("figures/kelptimeseries.png")





