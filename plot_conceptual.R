
library(tidyverse)


dist <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv") %>% 
  select(site, transect, perturbations, timesincelastextinct, extinctions, d, pv, median.canopy, mean.canopy_previousyear) %>%
  distinct()

df <- read.csv("data/Large_files/Raine_modeloutput/FullModelRunsAllSites.csv") %>% 
  pivot_longer(cols = ABUR1:WOOD3) %>%
  group_by(name) %>%
  separate(name, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  mutate(transect = as.numeric(transect)) %>%
  left_join(dist) %>%
  mutate(id = paste(site, transect, sep = "_"))

filt <- df %>% 
  filter(extinctions == max(extinctions, na.rm= T) | extinctions == min(extinctions, na.rm = T)) %>%
  distinct(id, extinctions)


df %>% 
  filter(id %in% as.vector(filt$id)) %>%
  ggplot(aes(x = year.quarter, y = value))+
  geom_line(aes(color = group))+
  facet_wrap(~id)


