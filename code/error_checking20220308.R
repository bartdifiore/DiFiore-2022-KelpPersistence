df <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv")



df %>%
  select(year, site, transect, perturbations, cv.canopy) %>%
  distinct() %>%
  arrange(year, site, transect)

# Errors to fix!!!!
    # STEP transect 1 is missing data for the persistence metrics ??? maybe the buffered region fell outside of the patch? When I view in Qgis it appears that the coordinates for STEP 1 may be incorrect. I will need to cross check with the original datafile. 

  # The extinction forloop is really funky and not working correctly... need to fix it. Seem scrap code below.

long %>%
  filter(year >= 2008) %>%
  filter(id == "AQUE3") %>% View()


df <- read.csv("data/spatial/csvs/timeseries-bysite.csv") %>%
  mutate(across(Y1984.1:Y2018.12,as.numeric)) %>%
  mutate(across(Y1984.1:Y2018.12, ~na_if(., "NaN"))) %>%
  left_join(coords)
  
View(ts)




ts %>% 
  group_by(id) %>%
  summarize(across(Y1984.1:Y2018.12, ~mean(., na.rm = T))) %>%
  mutate(across(Y1984.1:Y2018.12, ~na_if(., "NaN"))) %>%
  pivot_longer(cols = Y1984.1:Y2018.12, names_to = "time", values_to = "measurement") %>%
  tidyr::separate(time, into = c("junk", "year", "month"), sep = "Y|[.]") %>%
  mutate(junk = NULL, 
         month = as.integer(month)) %>%
  left_join(formerge) %>%
  group_by(year, quarter, id) %>%
  filter(year >= 2008) %>%
  filter(id == "AQUE3") %>% View()
  summarize(biomass = mean(measurement, na.rm = T)) %>% View()

  
  
  
  
  
  
df %>%
  group_by(coarse_grouping, site, transect, pv, d) %>%
  summarize(total_biomass = sum(dry_gm2, na.rm = T)) %>%
  ggplot(aes(x = pv, y = total_biomass))+
  geom_point(aes(color = site))+
  facet_wrap(~coarse_grouping, scales = "free")
  
  
df %>%
  group_by(coarse_grouping, site, transect, pv, d) %>%
  summarize(total_biomass = sum(dry_gm2, na.rm = T)) %>%
  ggplot(aes(x = d, y = total_biomass))+
  geom_point(aes(color = site))+
  facet_wrap(~coarse_grouping, scales = "free")

df %>%
  group_by(coarse_grouping, site, transect, pv, d) %>%
  summarize(total_biomass = sum(dry_gm2, na.rm = T)) %>%
  filter(coarse_grouping %in% c("SESSILE INVERT", "UNDERSTORY ALGAE")) %>%
  ggplot(aes(x = d, y = total_biomass))+
  geom_point(aes(color = site))+
  geom_smooth(method = "lm")+
  facet_wrap(~coarse_grouping, scales = "free")


library(lme4)
library(lmerTest)

formod <- df %>%
  group_by(coarse_grouping, site, transect, pv, d, mean.canopy, timesincelastextinct) %>%
  summarize(total_biomass = sum(dry_gm2, na.rm = T))

mod1 <- lmer(total_biomass ~ d + I(d^2) + mean.canopy + timesincelastextinct + (1|site), data = formod %>% filter(coarse_grouping == "UNDERSTORY ALGAE"))
summary(mod1)  

mod2 <- lmer(total_biomass ~ d + I(d^2) + mean.canopy + timesincelastextinct + (1|site), data = formod %>% filter(coarse_grouping == "SESSILE INVERT"))
summary(mod2)  
  




# I found that BB benthic data has duplicate entries for species. So I'm trying to figure out why

temp <- read.csv("~/Downloads/benthic_kelppersistence_20210908.csv")



