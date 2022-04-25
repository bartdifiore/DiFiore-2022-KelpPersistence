# Castorani et al. 2021

#  "We quantified the combined effects of sea urchin density and sand cover on macroalgal biomass using a Gamma GLM with a log link function. We then made predictions from this model using measurements of sea urchin density and sand cover to create a continuous macroalgal habitat quality index (rescaled to a maximum of one), and used a linear model to compare values of this variable among sites."

# The goal is to replicate how Castorani et al. 2021 estimate habitat quality, and include habitat quality in the model as a preditor variable. For example:
      # response ~ historic variability + habitat quality + current kelp biomass + timesincelast disturbance 
# I suspect that the lack of evidence for historic variabiltiy impacting community structure are differences in habitat quality. My hope is that my accounting for those differences I can pull out the historic effects. 

# However, as I think through this is doesn't make sense to construct the habitat quality variable based on how macroalgae (both kelp and understory) responds to sand and urchins, and then use this predictor to model understory algae. So maybe either including sand and urchins in the model is sufficient???

library(tidyverse)

urc <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv") %>% 
  filter(group == "mobile_invert", common_name %in% c("Purple Urchin", "Red Urchin")) %>%
  group_by(site, transect) %>%
  summarize(total_urc_biomass = sum(dry_gm2, na.rm = T))

sand <- read.csv("data/intermediary/combined_substrate.csv") %>%
  rename_all(tolower) %>%
  filter(!site %in% c("AHND", "SCDI", "SCTW")) %>%
  group_by(site, transect) %>%
  mutate(rocknot = ifelse(substrate_type %in% c("B", "BL", "BM", "BS", "C"), "rock", "sand")) %>%
  group_by(site, transect, rocknot) %>%
  summarize(sand_cover = sum(percent_cover)) %>%
  filter(rocknot == "sand") %>%
  dplyr::select(-rocknot)

df <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv") %>% 
  filter(group %in% c("kelp", "ua")) %>%
  group_by(site, transect) %>%
  summarize(total_macro_biomass = sum(dry_gm2, na.rm = T)) %>% #Estimate total biomass of macroalgae
  left_join(urc) %>%
  left_join(sand) %>%
  filter(!site %in% c("AHND", "SCDI", "SCTW")) %>%
  mutate(id = paste(site, transect, sep = "_"))

# What could total_macro_biomass be given observed urchin biomass and the % cover of sand


mod1 <- glm(total_macro_biomass ~ total_urc_biomass + sand_cover, df, family = Gamma(link = "log"))
summary(mod1)


p1 <- ggeffects::ggpredict(mod1)  
plot(p1, add.data = T)


  
df %>%
  mutate(id = fct_reorder(as.factor(id), total_urc_biomass)) %>%
  ggplot(aes(x = id, y = total_urc_biomass))+
  geom_point()
  
ggplot(df, aes(x = total_urc_biomass, y = total_macro_biomass))+
  geom_point()
ggplot(df, aes(x = sand_cover, y = total_macro_biomass))+
  geom_point()




predictions <- df %>% 
  select(site, transect)
predictions$hq <- predict(mod1, type = "response")
predictions$quality <- cut(predictions$hq, quantile(predictions$hq, prob = c(0, 0.33, 0.66, 1)), labels = c("low", "medium", "high"), include = T)

predictions %>% 
  arrange(hq) %>%
  ggplot(aes(x = as.numeric(as.factor(paste(site, transect))), y = hq))+
  geom_point(aes(color = quality))

mod2 <- aov(hq ~ site, predictions)
summary(mod2)

TukeyHSD(mod2)


split(predictions$hq, cut(predictions$hq, quantile(predictions$hq, prob = c(0, 0.33, 0.66, 1))))





summary(lm(total_urc_biomass ~ sand_cover, df))




