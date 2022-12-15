source("code/1_setup.R")

d1 <- read.csv("data/Large_files/Raine_modeloutput/EndModelRunsAllSites.csv") %>%
  pivot_longer(cols = c(ABUR1:WOOD3), names_to = "id", values_to = "predicted_kelpbiomass")

d2 <- read.csv("data/Large_files/Raine_modeloutput/EndModelRunsAllSitesObsKelp.csv")%>%
  pivot_longer(cols = c(ABUR1:WOOD3), names_to = "id", values_to = "predicted_kelpbiomass")



d1 %>%
  separate(id, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  ggplot(aes(x = tmpt, y = predicted_kelpbiomass))+
  geom_line(aes(color = group, linetype = as.factor(transect)))+
  facet_wrap(~site)




dt <- read.csv("data/Large_files/Raine_modeloutput/EndModelRunsAllSitesObsKelp.csv")%>%
  pivot_longer(cols = c(ABUR1:WOOD3), names_to = "id", values_to = "predicted_biomass") %>%
  select(group, month, year.quarter, id, predicted_biomass) %>%
  separate(id, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  filter(year.quarter == 2018.50) %>%
  group_by(site, transect, group) %>%
  summarize(predicted_percover = mean(predicted_biomass)) %>%
  mutate(group = case_when(group == "algae" ~ "ua", 
                           group == "invert" ~ "epiSI"), 
         transect = as.numeric(transect))


# get kelp maximums for the previous decade to relativize the observations

kelp <- read.csv("data/intermediary/kelpcanopytimeseries_wide.csv") %>% 
  pivot_longer(cols = ABUR1:WOOD3) %>%
  group_by(name) %>%
  filter(value == max(value)) %>%
  select(name, value) %>%
  separate(name, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  mutate(transect = as.numeric(transect)) %>%
  rename(max_kelp = value)

df <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv") %>% 
  filter(group %in% c("ua", "epiSI")) %>%
  group_by(site, transect, group, perturbations, timesincelastextinct, extinctions, d, pv, median.canopy, mean.canopy_previousyear) %>%
  summarize(biomass = sum(dry_gm2, na.rm = T)) %>% #Estimate total biomass of macroalgae
  left_join(kelp) %>%
  mutate(rel_biomass = biomass / max_kelp) %>%
  left_join(dt)%>%
  rename(obs_biomass = biomass)

plot(rel_biomass ~ predicted_biomass, df)


ggplot(df, aes(x = predicted_biomass, y = rel_biomass))+
  geom_point()+
  geom_abline(slope =1, intercept = 0, lty = 4)+
  facet_wrap(~group, scales = "free")+
  coord_cartesian(xlim = c(0, 0.75))


df %>%
  group_by(site, transect, group) %>%
  pivot_longer(cols = c(rel_biomass, predicted_biomass)) %>%
  ggplot(aes(x = group, y = value)) +
  geom_point(aes(color = name, group = group))+
  facet_wrap(~site+transect)



ggplot(df, aes(x = rel_biomass))+
  geom_density(aes(fill = group), alpha = 0.5)
  
  
  
  
  
  
  
forplot <- df %>% 
  ungroup() %>%
  select(site, transect, group, rel_biomass, predicted_biomass) %>%
  group_by(site, transect) %>%
  pivot_wider(id_cols = c(site, transect), names_from = group, values_from = c(rel_biomass, predicted_biomass)) %>%
  mutate(rel_ratio = rel_biomass_ua/rel_biomass_epiSI, 
         pred_ratio = predicted_biomass_ua/predicted_biomass_epiSI) %>%
  pivot_longer(cols = c(rel_ratio, pred_ratio)) %>%
  select(site, transect, name, value) %>%
  group_by(name) %>%
  mutate(id = paste(site, transect, sep = "_"))

ggplot(forplot, aes( x= id, y = value))+
  geom_point(aes(color = name))+
  coord_cartesian(ylim = c(-10, 10))


hist(forplot$value[forplot$name == "pred_ratio"])  
hist(forplot$value[forplot$name == "rel_ratio" & forplot$value < 500], breaks = 100)
forplot$value[forplot$name == "rel_ratio" & forplot$value < 1] 









df <- read.csv("data/raw/BartsBenthic_All_Species_Biomass_at_transect_20220424.csv", na.strings = "-99999") %>%
  as_tibble() %>%
  rename_all(tolower) %>%
  filter(coarse_grouping == "UNDERSTORY ALGAE" | coarse_grouping == "SESSILE INVERT") %>%
  group_by(site, transect, coarse_grouping) %>%
  summarize(obs_percover = sum(percent_cover, na.rm = T)) %>%
  mutate(coarse_grouping = case_when(coarse_grouping == "UNDERSTORY ALGAE" ~ "ua", 
                           coarse_grouping == "SESSILE INVERT" ~ "epiSI"), 
         obs_percover  = obs_percover/100) %>%
  rename(group = coarse_grouping) %>%
  pivot_wider(values_from = obs_percover, names_from = group) %>%
  mutate(obs_ratio = epiSI/ua)%>%
  select(-c(ua, epiSI))

dt <- read.csv("data/Large_files/Raine_modeloutput/EndModelRunsAllSitesObsKelp.csv")%>%
  pivot_longer(cols = c(ABUR1:WOOD3), names_to = "id", values_to = "predicted_biomass") %>%
  select(group, month, year.quarter, id, predicted_biomass) %>%
  separate(id, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  filter(year.quarter == 2018.50) %>%
  group_by(site, transect, group) %>%
  summarize(predicted_percover = mean(predicted_biomass)) %>%
  mutate(group = case_when(group == "algae" ~ "ua", 
                           group == "invert" ~ "epiSI"), 
         transect = as.numeric(transect)) %>%
  pivot_wider(values_from = predicted_percover, names_from = group) %>%
  mutate(predicted_ratio = epiSI/ua) %>%
  select(-c(ua, epiSI))

df <- df %>% 
  left_join(dt)


df %>%
  mutate(id = paste(site, transect, sep = "-")) %>%
  ungroup() %>%
  select(-c(site, transect)) %>%
  pivot_longer(cols = obs_ratio:predicted_ratio)%>%
  ggplot(aes(x = id, y = value))+
  geom_point(aes(color = name))

ggplot(df, aes(y = obs_ratio, x = predicted_ratio))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, lty = 4)+
  coord_cartesian(ylim = c(0,1), xlim = c(0, 1))+
  theme_classic()

plot(obs_percover ~ predicted_percover, df)
abline(0, 1)

df %>%
  group_by(site, transect) 












df <- read.csv("data/raw/BartsBenthic_All_Species_Biomass_at_transect_20220424.csv", na.strings = "-99999") %>%
  as_tibble() %>%
  rename_all(tolower) %>%
  filter(coarse_grouping == "UNDERSTORY ALGAE" | coarse_grouping == "SESSILE INVERT") %>%
  group_by(site, transect, coarse_grouping) %>%
  summarize(obs_percover = sum(percent_cover, na.rm = T)) %>%
  mutate(coarse_grouping = case_when(coarse_grouping == "UNDERSTORY ALGAE" ~ "ua", 
                                     coarse_grouping == "SESSILE INVERT" ~ "epiSI"), 
         obs_percover  = obs_percover/100) %>%
  rename(group = coarse_grouping)

dt <- read.csv("data/Large_files/Raine_modeloutput/EndModelRunsAllSitesObsKelp.csv")%>%
  pivot_longer(cols = c(ABUR1:WOOD3), names_to = "id", values_to = "predicted_biomass") %>%
  select(group, month, year.quarter, id, predicted_biomass) %>%
  separate(id, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  filter(year.quarter == 2018.50) %>%
  group_by(site, transect, group) %>%
  summarize(predicted_percover = mean(predicted_biomass)) %>%
  mutate(group = case_when(group == "algae" ~ "ua", 
                           group == "invert" ~ "epiSI"), 
         transect = as.numeric(transect))

df %>%
  left_join(dt) %>%
  mutate(id = paste(site, transect, sep = "-")) %>%
  ungroup() %>%
  select(-c(site, transect)) %>%
  pivot_longer(cols = c(predicted_percover, obs_percover)) %>%
  ggplot(aes(x = id, y = value))+
  geom_point(aes(color = name))+
  facet_wrap(~group)
  









