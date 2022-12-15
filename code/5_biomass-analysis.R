source("code/1_setup.R")
source("code/theme.R")

ua_col = "chocolate1"
si_col = "royalblue1"
kelp_col = "darkgreen"


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

df <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv", na.strings = "-99999") %>% 
  filter(group %in% c("ua", "epiSI")) %>%
  group_by(site, transect, group, perturbations, timesincelastextinct, timesincelastabsent, extinctions, d, pv, median.canopy, mean.canopy_previousyear, prop_zero, prop_50) %>%
  summarize(biomass = sum(dry_gm2, na.rm = T)) %>% #Estimate total biomass of macroalgae
  left_join(urc) %>%
  left_join(sand) %>%
  left_join(kelp) %>%
  mutate(id = paste(site, transect, sep = "_"))



df %>%
  filter(group == "ua") %>%
  ggplot(aes(x = perturbations, y = biomass))+
  geom_point(aes(color = site))+
  facet_wrap(~site)


mod1 <- rstanarm::stan_glmer(biomass ~ scale(perturbations) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + scale(timesincelastabsent) + scale(prop_zero) + (1|site), df[df$group == "epiSI", ], family = Gamma(link = "log"))


rstanarm::prior_summary(mod1)
print(summary(mod1), digits = 3)

coef(summary(mod1))
coef(mod1)
plot(ggeffects::ggpredict(mod1), add.data = T)
temp <- rstantools::posterior_predict(mod1, re.form = NA)


mod2 <- rstanarm::stan_glmer(biomass ~ scale(perturbations) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + scale(timesincelastabsent) + scale(prop_zero) + (1|site), df[df$group == "ua", ], family = Gamma(link = "log"))


print(summary(mod2), digits = 3)


# general conclusions after first pass of the analysis: 

# 1. Understory algae biomass is highly dependent on the time since kelp last went extinct, urchin biomass, and current kelp cover. This may be because UA life history is relatively fast (annual or semi-annual species). The number of extinctions or perturbations does not seem to have a large impact. 
# 
# 2. Sessile invert biomass is highly dependent on the number of perturbations (declines as disturbance increases) and sand cover. However, sessile invert biomass is not dependent on the time since the last disturbance, likely due to slower life history strategies relative to UA's.




p1 <- mod1 %>%
  gather_draws(`scale(perturbations)`,`scale(sand_cover)`, `scale(total_urc_biomass)`, `scale(current_kelpbiomass)`, `scale(timesincelastabsent)`, `scale(prop_zero)`) %>%
  median_qi(.width = c(0.75, 0.95)) %>%
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper))+
  geom_pointinterval()+
  geom_vline(xintercept = 0, lty = 4)+
  scale_y_discrete(labels = rev(c("Urchin biomass", "Time since\nlast absent", "Sand cover", "Prop. zero", "Perturbations", "Current kelp biomass")))+
  labs(x = "Effect on sessile invertebrate biomass", y = "")+
  theme_classic()


p2 <- mod2 %>%
  gather_draws(`scale(perturbations)`,`scale(sand_cover)`, `scale(total_urc_biomass)`, `scale(current_kelpbiomass)`, `scale(timesincelastabsent)`, `scale(prop_zero)`) %>%
  median_qi(.width = c(0.75, 0.95)) %>%
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper))+
  geom_pointinterval()+
  geom_vline(xintercept = 0, lty = 4)+
  scale_y_discrete(labels = rev(c("Urchin biomass", "Time since\nlast absent", "Sand cover", "Prop. zero", "Perturbations", "Current kelp biomass")))+
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


#--------------------------
## One panel version
#--------------------------

si <- mod1 %>%
  gather_draws(`scale(perturbations)`,`scale(sand_cover)`, `scale(total_urc_biomass)`, `scale(current_kelpbiomass)`, `scale(timesincelastabsent)`, `scale(prop_zero)`) %>%
  median_qi(.width = c(0.75, 0.95)) %>%
  mutate(group = "Sessile\ninvertebrates")

ua <- mod2 %>%
  gather_draws(`scale(perturbations)`,`scale(sand_cover)`, `scale(total_urc_biomass)`, `scale(current_kelpbiomass)`, `scale(timesincelastabsent)`, `scale(prop_zero)`) %>%
  median_qi(.width = c(0.75, 0.95)) %>%
  mutate(group = "Understory\nalgae")


ua %>% 
  bind_rows(si) %>%
  mutate(predictors = forcats::fct_rev(case_when(.variable == "scale(current_kelpbiomass)" ~ "D Current kelp\nbiomass",
                                .variable == "scale(perturbations)" ~ "B Perterbations",
                                .variable == "scale(prop_zero)" ~ "A Proportion time\nno kelp",
                                .variable == "scale(sand_cover)" ~ "E Sand cover",
                                .variable == "scale(timesincelastabsent)" ~ "C Time since\nlast absent",
                                .variable == "scale(total_urc_biomass)" ~ "F Urchin biomass"
                                ))) %>%
  ggplot(aes(y = predictors, x = .value, xmin = .lower, xmax = .upper))+
  geom_pointinterval(aes(color = group), position = position_dodge(width = 0.5))+
  scale_color_manual(values = c(si_col, ua_col))+
  geom_vline(xintercept = 0, lty = 4)+
  scale_y_discrete(labels = rev(c("Proportion time\nno kelp", "Perturbations", "Time since\nlast absent", "Current kelp\nbiomass", "Sand cover", "Urchin biomass")))+
  labs(x = "Effect on biomass", y = "", color = "")+
  theme_bd()+
  theme(legend.position = c(0.8,0.2))

ggsave("figures/coef_plot.png", width = 7, height = 6)


#----------------------------------
## Time Since last absention
#----------------------------------
dist <- read.csv("data/intermediary/persistence_metrics.csv")

rd_ts  <- read.csv("data/Large_files/Raine_modeloutput/FullModelRunsAllSites.csv") %>%
  pivot_longer(cols = ABUR1:WOOD3, names_to = "id", values_to = "biomass") %>%
  group_by(id) %>%
  separate(id, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])", remove = F) %>%
  mutate(transect = as.numeric(transect)) %>%
  left_join(dist) %>%
  separate(year.quarter, into = c("year", "quarter"), sep = "[.]") %>%
  mutate(quarter = case_when(quarter == 1 ~ ".0", 
                             quarter == 2 ~ ".25", 
                             quarter == 3 ~ ".5", 
                             quarter == 4 ~ ".75")) %>%
  mutate(year.quarter = as.numeric(paste(year, quarter, sep = "")))


rd_predictions <- rd_ts %>% 
  filter(year.quarter == 2018.5, group != "kelp")

p1 <- ggplot(rd_predictions, aes(x = timesincelastabsent, y = biomass))+
  geom_point(aes(color = group))+
  geom_smooth(method = "lm", aes(color = group))+
  scale_color_manual(values = c(ua_col, si_col))+
  labs(x = "", y = "Theoretical prediction\nrelative abundance", color = "")+
  theme_bd()+
  theme(legend.position = "top")

lm1 <- lm(biomass ~ timesincelastabsent*group, rd_predictions)
summary(lm1)

p2 <- df %>%
  filter(group == "ua") %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(timesincelastabsent = modelr::seq_range(timesincelastabsent, n = 100), .model = mod2) %>%
  add_predicted_draws(mod2, re_formula = NA) %>%
  ggplot(aes(x = timesincelastabsent, y = biomass)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray80", color = ua_col) +
  geom_point(data = filter(df, group == "ua"), aes(x = timesincelastabsent, y = biomass), color = ua_col)+
  theme_bd()+
  labs(x = "", y = "Understory algae biomass")

p3 <- df %>%
  filter(group == "epiSI") %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(timesincelastabsent = modelr::seq_range(timesincelastabsent, n = 100), .model = mod1) %>%
  add_predicted_draws(mod1, re_formula = NA) %>%
  ggplot(aes(x = timesincelastabsent, y = biomass)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray80", color = si_col) +
  geom_point(data = filter(df, group == "epiSI"), aes(x = timesincelastabsent, y = biomass), color = si_col)+
  theme_bd()+
  labs(x = "Time since last kelp canopy extinction", y = "Sessile invertebrate biomass")

cowplot::plot_grid(p1, p2, p3, nrow = 3)

ggsave("figures/time_since.png", width = 6, height = 12 )



#----------------------------------
## Proportion of time no kelp
#----------------------------------


p4 <- ggplot(rd_predictions, aes(x = prop_zero, y = biomass))+
  geom_point(aes(color = group))+
  geom_smooth(method = "lm", aes(color = group))+
  scale_color_manual(values = c(ua_col, si_col))+
  labs(x = "", y = "Theoretical prediction\nrelative abundance", color = "")+
  theme_bd()+
  theme(legend.position = "top")

lm1 <- lm(biomass ~ prop_zero*group, rd_predictions)
summary(lm1)

p5 <- df %>%
  filter(group == "ua") %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(prop_zero = modelr::seq_range(prop_zero, n = 100), .model = mod2) %>%
  add_predicted_draws(mod2, re_formula = NA) %>%
  ggplot(aes(x = prop_zero, y = biomass)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray80", color = ua_col) +
  geom_point(data = filter(df, group == "ua"), aes(x = prop_zero, y = biomass), color = ua_col)+
  theme_bd()+
  labs(x = "", y = "Understory algae biomass")

p6 <- df %>%
  filter(group == "epiSI") %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(prop_zero = modelr::seq_range(prop_zero, n = 100), .model = mod1) %>%
  add_predicted_draws(mod1, re_formula = NA) %>%
  ggplot(aes(x = prop_zero, y = biomass)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray80", color = si_col) +
  geom_point(data = filter(df, group == "epiSI"), aes(x = prop_zero, y = biomass), color = si_col)+
  theme_bd()+
  labs(x = "Proportion of time kelp was absent", y = "Sessile invertebrate biomass")

cowplot::plot_grid(p4, p5, p6, nrow = 3)

ggsave("figures/prop_zero.png", width = 6, height = 12 )


#----------------------------------
## Number of perturbations
#----------------------------------


p7 <- ggplot(rd_predictions, aes(x = perturbations, y = biomass))+
  geom_point(aes(color = group))+
  geom_smooth(method = "lm", aes(color = group))+
  scale_color_manual(values = c(ua_col, si_col))+
  labs(x = "", y = "Theoretical prediction\nrelative abundance", color = "")+
  theme_bd()+
  theme(legend.position = "top")

lm1 <- lm(biomass ~ perturbations*group, rd_predictions)
summary(lm1)

p8 <- df %>%
  filter(group == "ua") %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(perturbations = modelr::seq_range(perturbations, n = 100), .model = mod2) %>%
  add_predicted_draws(mod2, re_formula = NA) %>%
  ggplot(aes(x = perturbations, y = biomass)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray80", color = ua_col) +
  geom_point(data = filter(df, group == "ua"), aes(x = perturbations, y = biomass), color = ua_col)+
  theme_bd()+
  labs(x = "", y = "Understory algae biomass")

p9 <- df %>%
  filter(group == "epiSI") %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(perturbations = modelr::seq_range(perturbations, n = 100), .model = mod1) %>%
  add_predicted_draws(mod1, re_formula = NA) %>%
  ggplot(aes(x = perturbations, y = biomass)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray80", color = si_col) +
  geom_point(data = filter(df, group == "epiSI"), aes(x = perturbations, y = biomass), color = si_col)+
  theme_bd()+
  labs(x = "Number of disturbances", y = "Sessile invertebrate biomass")

cowplot::plot_grid(p7, p8, p9, nrow = 3)

ggsave("figures/perturbations.png", width = 6, height = 12 )

#----------------------------------
## Six panel version
#----------------------------------

cowplot::plot_grid(p1, p2, p3, p4+labs(y = ""), p5+labs(y = ""), p6+labs(y = ""), nrow = 3, ncol = 2, byrow = F, labels = "auto")

ggsave("figures/six_panel_fig3.png", width = 12, height = 12 )


#-----------------------------------
## JUNK
#-----------------------------------

temp <- ua %>% 
  bind_rows(si)

mod_glmm <- glmmTMB(biomass ~ scale(perturbations) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + scale(timesincelastabsent) + (1|site), df[df$group == "epiSI", ], family = Gamma(link = "log"))

mod1 <- rstanarm::stan_glmer(biomass ~ scale(d) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + scale(timesincelastextinct) + (1|site), df[df$group == "epiSI", ], family = Gamma(link = "log"))


print(summary(mod1), digits = 3)


mod2 <- rstanarm::stan_glmer(biomass ~ scale(extinctions) + scale(timesincelastextinct) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + (1|site), df[df$group == "ua", ], family = Gamma(link = "log"))

print(summary(mod2), digits = 3)

mod3 <- rstanarm::stan_glmer(biomass ~ scale(extinctions) + scale(timesincelastextinct) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + (1|site), df[df$group == "epiSI", ], family = Gamma(link = "log"))

print(summary(mod3), digits = 3)



mod3 %>%
  gather_draws(`scale(extinctions)`,`scale(sand_cover)`, `scale(total_urc_biomass)`, `scale(current_kelpbiomass)`, `scale(timesincelastextinct)`) %>%
  median_qi(.width = c(0.75, 0.95)) %>%
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper))+
  geom_pointinterval()+
  geom_vline(xintercept = 0, lty = 4)+
  scale_y_discrete(labels = rev(c("Urchin biomass", "Time since\nlast extinction", "Sand cover", "extinctions", "Current kelp biomass")))+
  labs(x = "Effect on sessile invertebrate biomass", y = "")+
  theme_classic()





