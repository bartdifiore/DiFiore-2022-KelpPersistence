source("code/1_setup.R")


df <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv")



df %>% 
  filter(group %in% c("epiSI", "ua", "kelp")) %>%
  group_by(site, transect, group, pv, d, perturbations) %>%
  summarize(biomass = sum(dry_gm2, na.rm = T)) %>%
  ggplot(aes(y = biomass, x = perturbations))+
  geom_jitter(shape = 1)+
  geom_smooth(method = "lm")+
  facet_wrap(~group, scales = "free")+
  labs(x = "Persistence metric (# of prolonged periods of kelp loss in the last 10 years)", y = expression(paste("Biomass (g"~m^{-2}~")")))+
  theme_pubr()


df %>% 
  filter(group %in% c("epiSI", "ua", "kelp")) %>%
  group_by(site, transect, group, pv, d, perturbations) %>%
  summarize(biomass = sum(dry_gm2, na.rm = T)) %>%
  ggplot(aes(y = biomass, x = pv))+
  geom_jitter(shape = 1)+
  geom_smooth(method = "lm")+
  facet_wrap(~group, scales = "free")+
  labs(x = "Persistence metric (# of prolonged periods of kelp loss in the last 10 years)", y = expression(paste("Biomass (g"~m^{-2}~")")))+
  theme_pubr()


df %>% 
  filter(group %in% c("epiSI", "ua", "kelp")) %>%
  group_by(site, transect, group, pv, d, perturbations) %>%
  summarize(biomass = sum(dry_gm2, na.rm = T)) %>%
  ggplot(aes(y = biomass, x = d))+
  geom_jitter(shape = 1)+
  geom_smooth(method = "lm")+
  facet_wrap(~group, scales = "free")+
  labs(x = "Persistence metric (# of prolonged periods of kelp loss in the last 10 years)", y = expression(paste("Biomass (g"~m^{-2}~")")))+
  theme_pubr()


# modeling

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

kelp <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv") %>%
  filter(group %in% c("kelp")) %>%
  mutate(current_kelpbiomass = dry_gm2) %>%
  select(site, transect, current_kelpbiomass)

df <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv") %>% 
  filter(group %in% c("ua", "epiSI")) %>%
  group_by(site, transect, group, perturbations, timesincelastextinct, extinctions, d, pv, median.canopy, mean.canopy_previousyear) %>%
  summarize(biomass = sum(dry_gm2, na.rm = T)) %>% #Estimate total biomass of macroalgae
  left_join(urc) %>%
  left_join(sand) %>%
  left_join(kelp) %>%
  mutate(id = paste(site, transect, sep = "_"))


mod1 <- glmmTMB::glmmTMB(biomass ~ scale(perturbations) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + (1|site), df[df$group == "ua", ], family = Gamma(link = "log"))
summary(mod1)
hist(residuals(mod1))
car::qqPlot(residuals(mod1))

plot(ggeffects::ggpredict(mod1), add.data = T)


mod1 <- glm(biomass ~ perturbations + sand_cover + total_urc_biomass, df[df$group == "ua", ], family = Gamma(link = "log"))
summary(mod1)


df %>%
  filter(group == "ua") %>%
  ggplot(aes(x = perturbations, y = biomass))+
  geom_point(aes(color = site))+
  facet_wrap(~site)


mod1 <- rstanarm::stan_glmer(biomass ~ scale(perturbations) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + scale(timesincelastextinct) + (1|site), df[df$group == "epiSI", ], family = Gamma(link = "log"))


print(summary(mod1), digits = 3)

coef(summary(mod1))
coef(mod1)
plot(ggeffects::ggpredict(mod1), add.data = T)
temp <- rstantools::posterior_predict(mod1, re.form = NA)


mod2 <- rstanarm::stan_glmer(biomass ~ scale(perturbations) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + scale(timesincelastextinct) + (1|site), df[df$group == "ua", ], family = Gamma(link = "log"))


print(summary(mod2), digits = 3)


# general conclusions after first pass of the analysis: 

# 1. Understory algae biomass is highly dependent on the time since kelp last went extinct, urchin biomass, and current kelp cover. This may be because UA life history is relatively fast (annual or semi-annual species). The number of extinctions or perturbations does not seem to have a large impact. 
# 
# 2. Sessile invert biomass is highly dependent on the number of perturbations (declines as disturbance increases) and sand cover. However, sessile invert biomass is not dependent on the time since the last disturbance, likely due to slower life history strategies relative to UA's.




p1 <- mod1 %>%
  gather_draws(`scale(perturbations)`,`scale(sand_cover)`, `scale(total_urc_biomass)`, `scale(current_kelpbiomass)`, `scale(timesincelastextinct)`) %>%
  median_qi(.width = c(0.75, 0.95)) %>%
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper))+
  geom_pointinterval()+
  geom_vline(xintercept = 0, lty = 4)+
  scale_y_discrete(labels = rev(c("Urchin biomass", "Time since\nlast extinction", "Sand cover", "Perturbations", "Current kelp biomass")))+
  labs(x = "Effect on sessile invertebrate biomass", y = "")+
  theme_classic()


p2 <- mod2 %>%
  gather_draws(`scale(perturbations)`,`scale(sand_cover)`, `scale(total_urc_biomass)`, `scale(current_kelpbiomass)`, `scale(timesincelastextinct)`) %>%
  median_qi(.width = c(0.75, 0.95)) %>%
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper))+
  geom_pointinterval()+
  geom_vline(xintercept = 0, lty = 4)+
  scale_y_discrete(labels = rev(c("Urchin biomass", "Time since\nlast extinction", "Sand cover", "Perturbations", "Current kelp biomass")))+
  labs(x = "Effect on understory macroalgae biomass", y = "")+
  theme_classic()


p3 <- df %>%
  filter(group == "epiSI") %>%
  filter(!site %in% c("AHND", "SCDI", "SCTW")) %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(perturbations = modelr::seq_range(perturbations, n = 100), .model = mod1) %>%
  add_predicted_draws(mod1, re_formula = NA) %>%
  ggplot(aes(x = perturbations, y = biomass)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray", color = "black") +
  geom_jitter(data = filter(df, group == "epiSI"), aes(x = perturbations, y = biomass))+
  # coord_cartesian(ylim = c(0, 250))+
  theme_classic()+
  labs(x = "Perturbations", y = "Sessile invertebrate biomass")
  

p4 <- df %>%
  filter(group == "ua") %>%
  filter(!site %in% c("AHND", "SCDI", "SCTW")) %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(timesincelastextinct = modelr::seq_range(timesincelastextinct, n = 100), .model = mod2) %>%
  add_predicted_draws(mod2, re_formula = NA) %>%
  ggplot(aes(x = timesincelastextinct, y = biomass)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray", color = "black") +
  geom_jitter(data = filter(df, group == "ua"), aes(x = timesincelastextinct, y = biomass))+
  theme_classic()+
  labs(x = "Time since last kelp canopy extinction", y = "Understory algae biomass")

fourpanel <- cowplot::plot_grid(p1, p3, p2, p4)
fourpanel

ggsave("figures/4panel_bayesresults.png", fourpanel, width = 12, height = 8)







mod1 <- rstanarm::stan_glmer(biomass ~ scale(d) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + scale(timesincelastextinct) + (1|site), df[df$group == "epiSI", ], family = Gamma(link = "log"))


print(summary(mod1), digits = 3)


mod2 <- rstanarm::stan_glmer(biomass ~ scale(perturbations)*scale(timesincelastextinct) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + (1|site), df[df$group == "ua", ], family = Gamma(link = "log"))


print(summary(mod2), digits = 3)









