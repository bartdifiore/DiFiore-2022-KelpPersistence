source("code/1_setup.R")
source("code/2_dist_function.R")
source("code/theme.R")


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

# Number of extinctions and time when the extinction started
mat <- as.matrix(wide[,-1])

temp <- matrix(nrow = dim(mat)[1], ncol = dim(mat)[2])

for(k in 1:dim(mat)[2]){
  for(i in 1:(dim(mat)[1]-2)){
    temp[i,k] <- ifelse(mat[i,k] > 0
                        & mat[i+1,k] == 0 
                        & mat[i+2,k] == 0, 1, 0)
  }
}

sum$extinctions <- apply(temp, MARGIN = 2, sum, na.rm = T)
map_extinctions <- data.frame(temp)
names(map_extinctions) <- names(wide[-1])
map_extinctions$year.quarter <- wide$year.quarter

# Time since kelp last went extinct

##### This code will show the last extinction event (kelp present then absent) BUT the last period when kelp was gone (regardless of if it was present before) may be the more important
time_since <- map_extinctions %>% 
  pivot_longer(cols = ABUR1:WOOD3, names_to = "id", values_to = "extinction_event") %>%
  group_by(id) %>%
  filter(extinction_event == 1) %>% 
  filter(year.quarter == max(year.quarter)) %>%
  mutate(timesincelastextinct = 2018.75 - year.quarter) %>% 
  select(id, timesincelastextinct)

sum <- sum %>% 
  left_join(time_since) %>%
  mutate(timesincelastextinct = replace_na(timesincelastextinct, 10))

# Time since kelp was zero for 6 months

mat <- as.matrix(wide[,-1])

temp <- matrix(nrow = dim(mat)[1], ncol = dim(mat)[2])

for(k in 1:dim(mat)[2]){
  for(i in 1:(dim(mat)[1]-1)){
    temp[i+1,k] <- ifelse(mat[i,k] == 0
                        & mat[i+1,k] == 0, 1, 0)
  }
}

map_absent <- data.frame(temp)
names(map_absent) <- names(wide[-1])
map_absent$year.quarter <- wide$year.quarter
timesince_absent <- map_absent %>%
  pivot_longer(cols = c(ABUR1:WOOD3), names_to = "id", values_to = "absences") %>%
  group_by(id) %>%
  filter(absences == 1) %>% 
  filter(year.quarter == max(year.quarter)) %>%
  mutate(timesincelastabsent = 2018.75 - year.quarter) %>% 
  select(id, timesincelastabsent)
  
sum <- sum %>% 
  left_join(timesince_absent) %>%
  mutate(timesincelastabsent = replace_na(timesincelastabsent, 10))

##############################################################################
## All disturbance events
##############################################################################

#number of times kelp has disappeared (< 80% reduction) for more than 2 quarters.

mat <- as.matrix(wide[,-1])

temp <- matrix(nrow = dim(mat)[1], ncol = dim(mat)[2])

threshold <- apply(mat, MARGIN = 2, FUN = function(x){max(x)*0.1})

for(k in 1:dim(mat)[2]){
  for(i in 1:(dim(mat)[1]-2)){
    temp[i,k] <- ifelse(mat[i,k] >= threshold[k]
                        & ((mat[i,k] - mat[i+1,k]) / mat[i,k]) >= 0.80 
                        & ((mat[i,k] - mat[i+2,k]) / mat[i,k]) >= 0.80, 1, 0)
  }
}

(1000 - 1100) / 1000

for(k in 1:dim(mat)[2]){
  for(i in 1:(dim(mat)[1]-1)){
    temp[i,k] <- ifelse(temp[i,k] == 1
                        & temp[i+1,k] == 1
                        & is.na(temp[i+1,k]) == F, NA, temp[i,k] )
  }
} # This for loop corrects for two disturbances in a row. Basically, I only want the first one when the decline was initiated, not a subsequent one in the next quarter.

sum$perturbations <- apply(temp, MARGIN = 2, sum, na.rm = T)
map_perturbations <- data.frame(temp)
names(map_perturbations) <- names(wide[-1])
map_perturbations$year.quarter <- wide$year.quarter


##############################################################################
## Proportional metrics
##############################################################################

percover <- read.csv("data/spatial/csvs/percover_timeseries-bysite.csv") %>%
  pivot_longer(cols = c(Y_2008_1:Y_2018_4), names_to = "time", values_to = "percent_cover") %>%
  group_by(id, site, transect, time) %>%
  summarize(num_pixels = n(), 
            percent_cover = sum(percent_cover, na.rm = T)/(num_pixels*900)) %>%
  separate(time, into = c("junk", "year", "quarter"), sep = "_") %>%
  select(-junk) %>%
  mutate(year.quarter = as.numeric(paste(year, quarter, sep = ".")))


ggplot(percover, aes(x = year.quarter, y = percent_cover))+
  geom_line(aes(color = as.character(transect)))+
  facet_wrap(~site)


# Proportion of time percent cover was zero

zero_percover <- percover %>%
  filter(percent_cover == 0) %>%
  summarize(prop_zero = n()/40)



# Proporiton of time kelp was > 80% areal coverage

fifty_percover <- percover %>%
  filter(percent_cover >= .5) %>%
  summarize(prop_50 = n()/40)



##############################################################################
## Merge into single data file
##############################################################################

persistence <- sum %>% 
  separate(id, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  mutate(transect = as.numeric(transect)) %>%
  left_join(zero_percover) %>%
  left_join(fifty_percover)

psych::pairs.panels(persistence %>% select(-id), ellipses = F)

psych::pairs.panels(persistence %>% select(mean.canopy, cv.canopy, perturbations, timesincelastextinct, extinctions), ellipses = F)

psych::pairs.panels(persistence %>% select(mean.canopy, cv.canopy, d, pv), ellipses = F)

psych::pairs.panels(persistence %>% select(perturbations, timesincelastextinct, extinctions, prop_zero, prop_50))

write.csv(persistence, "data/intermediary/persistence_metrics.csv", row.names = F)

###########################################################################
## Merge with community data by patch 
###########################################################################

# merge non-core sites with data from core sites in 2018

#LTER data
lt <- read.csv("data/raw/Annual_All_Species_Biomass_at_transect_20210108.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  filter(YEAR == 2018) %>%
  mutate(survey = "core") %>%
  rename_all(tolower)


sp.meta <- lt %>% 
  select(sp_code, scientific_name:coarse_grouping) %>%
  distinct()



# Noncore site data
df <- read.csv("data/raw/BartsBenthic_All_Species_Biomass_at_transect_20220424.csv", stringsAsFactors = F,na.strings ="-99999") %>% 
  mutate(survey = "noncore") %>%
  rename_all(tolower)



dat <- bind_rows(df, lt) %>%
  select(site:sp_code, dry_gm2, survey) %>%
  complete(sp_code, 
           nesting(site, transect, survey), 
           fill = list(dry_gm2 = 0)) %>%
  left_join(sp.meta) %>% 
  mutate(group = case_when(coarse_grouping == "SESSILE INVERT" & !sp_code %in% c("PACA", "CHOV") ~ "epiSI", 
                           coarse_grouping == "GIANT KELP" ~ "kelp",
                           coarse_grouping == "UNDERSTORY ALGAE" ~ "ua", 
                           coarse_grouping == "FISH" ~ "fish", 
                           coarse_grouping == "MOBILE INVERT" ~ "mobile_invert", 
                           coarse_grouping == "SESSILE INVERT" & sp_code %in% c("PACA", "CHOV") ~ "endoSI")) %>%
  select(c(sp_code:common_name, group)) %>%
  left_join(persistence) %>% 
  filter(!site %in% c("AHND", "SCDI", "SCTW"))


# find and filter out species that were never observed

noobs <- dat %>% group_by(sp_code) %>% 
  summarize(total = sum(dry_gm2)) %>% 
  filter(total == 0)

noobs.v <- as.vector(noobs$sp_code)

dat2 <- dat %>% 
  filter(!sp_code %in% noobs.v) %>%
  mutate(across(everything(), ~replace_na(.x, -99999)))


# There are 17 sites in the final dataset sampled in 2018!!!!

write.csv(dat2, "data/intermediary/species_withpersistence-pixelscale.csv", row.names = F)

write.csv(wide, "data/intermediary/kelpcanopytimeseries_wide.csv", row.names = F)

###########################################################################
## Clean up the substrate data and merge with the LTER core sites
###########################################################################


bbs <- read.csv("data/raw/substrate/BBenthic_substrate.csv") %>%
  dplyr::select(-DATE, -QUAD, -SIDE)


lts <- read.csv("data/raw/substrate/Annual_Substrate_All_Years.csv") %>%
  dplyr::select(-DATE, -QUAD, -SIDE) %>%
  filter(YEAR == 2018, SITE != "SCDI" | SITE != "SCTW")

df <- bind_rows(bbs, lts) %>% group_by(SITE, TRANSECT, SUBSTRATE_TYPE) %>%
  summarize(percent_cover = mean(PERCENT_COVER, na.rm = T))

write.csv(df, "data/intermediary/combined_substrate.csv", row.names = F)

###########################################################################
## Plot kelp biomass dynamics
###########################################################################


kelp_ts <- wide %>%
  pivot_longer(cols = c(ABUR1:WOOD3), names_to = "id", values_to = "biomass") %>%
  separate(id, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])", remove = F) %>%
  left_join(sum)

ggplot(kelp_ts, aes(x = year.quarter, y = biomass))+
  geom_line(aes(color = site, group = transect))+
  facet_wrap(~site)+
  theme_classic()
ggsave("figures/kelptimeseries.png")

temp <- map_extinctions %>% 
  pivot_longer(cols = c(ABUR1:WOOD3), names_to = "id", values_to = "extinction_events") %>%
  filter(extinction_events == 1)

ggplot(kelp_ts, aes(x = year.quarter, y = biomass))+
  geom_line(aes(group = id, color = id), show.legend = F)+
  geom_vline(data = temp, aes(xintercept = year.quarter[extinction_events == 1], group = id), color = "red", lty = 4)+
  facet_wrap(~id)+
  theme_bd()

ggsave("figures/kelptimeseries_wextinctions.png", width = 20, height = 20)

temp2 <- map_perturbations %>% 
  pivot_longer(cols = c(ABUR1:WOOD3), names_to = "id", values_to = "perturbation_events") %>%
  filter(perturbation_events == 1)

ggplot(kelp_ts, aes(x = year.quarter, y = biomass))+
  geom_line(aes(group = id, color = id), show.legend = F)+
  geom_vline(data = temp2, aes(xintercept = year.quarter[perturbation_events == 1], group = id), color = "red", lty = 4)+
  facet_wrap(~id)+
  theme_bd()

ggsave("figures/kelptimeseries_wperturbations.png", width = 20, height = 20)



###########################################################################
## Figures for Table S1 - Summary of persistence metrics
###########################################################################


ggplot(kelp_ts, aes(x = year.quarter, y = biomass))+
  geom_line(aes(color = site, group = transect))+
  facet_wrap(~forcats::fct_reorder(id, cv.canopy))+
  theme_classic()

kelp_ts %>% 
  filter(id == "OAKS3") %>% 
  ggplot(aes(x = year.quarter, y = biomass))+
  geom_line()+
  labs(x = "", y = "Biomass")+
  theme_classic()

ggsave("figures/tables1_1.png", width = 3, height = 2)


persistence %>%
  ggplot(aes(x = mean.canopy, y = sd.canopy))+
  geom_point()+ 
  labs(x = "Mean kelp biomass", y = "SD kelp biomass")+
  theme_classic()

ggsave("figures/tables1_2.png", width = 3, height = 2)


persistence %>%
  ggplot(aes(x = median.canopy))+
  geom_histogram(color = "white")+
  labs(x = "Median kelp biomass")+
  theme_classic()

ggsave("figures/tables1_3.png", width = 3, height = 2)


persistence %>%
  ggplot(aes(x = d))+
  geom_histogram(color = "white", bins = 15)+
  labs(x = "d")+
  theme_classic()

ggsave("figures/tables1_4.png", width = 3, height = 2)

persistence %>%
  ggplot(aes(x = pv))+
  geom_histogram(color = "white", bins = 15)+
  labs(x = "pv")+
  theme_classic()

ggsave("figures/tables1_5.png", width = 3, height = 2)

kelp_ts %>%
  filter(id == "CARP6") %>%
ggplot(aes(x = year.quarter, y = biomass))+
  geom_line()+
  geom_vline(data = filter(temp, id == "CARP6"), aes(xintercept = year.quarter[extinction_events == 1], group = id), color = "red", lty = 4)+
  labs(x = "", y = "Kelp biomass")+
  theme_classic()

ggsave("figures/tables1_6.png", width = 3, height = 2)


persistence %>%
  ggplot(aes(x = perturbations))+
  geom_histogram(color = "white", bins = 8)+
  labs(x = "Num. of perturbations")+
  theme_classic()

ggsave("figures/tables1_7.png", width = 3, height = 2)


kelp_ts %>%
  filter(id == "INNP3") %>%
  ggplot(aes(x = year.quarter, y = biomass))+
  geom_line()+
  geom_vline(data = filter(temp, id == "CARP6"), aes(xintercept = year.quarter[extinction_events == 1], group = id), color = "red", lty = 4)+
  labs(x = "", y = "Kelp biomass")+
  theme_classic()

ggsave("figures/tables1_8.png", width = 3, height = 2)


persistence %>%
  ggplot(aes(x = timesincelastabsent))+
  geom_histogram(color = "white", bins = 15)+
  labs(x = "Time since last absent")+
  theme_classic()

ggsave("figures/tables1_9.png", width = 3, height = 2)



persistence %>%
  ggplot(aes(x = prop_zero))+
  geom_histogram(color = "white", bins = 10)+
  labs(x = "Proportion of time kelp absent")+
  theme_classic()

ggsave("figures/tables1_10.png", width = 3, height = 2)


persistence %>%
  ggplot(aes(x = prop_50))+
  geom_histogram(color = "white", bins = 10)+
  labs(x = "Proportion of time kelp >50% of max")+
  theme_classic()

ggsave("figures/tables1_11.png", width = 3, height = 2)





















