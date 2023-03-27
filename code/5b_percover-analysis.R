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
    
    df <- read.csv("data/intermediary/percover_withpersitence.csv") %>% 
      left_join(urc) %>%
      left_join(sand) %>%
      left_join(kelp) %>%
      mutate(id = paste(site, transect, sep = "_")) %>%
      filter(!site %in% c("AHND", "SCDI", "SCTW")) %>%
      mutate(percent_cover = ifelse(percent_cover == 0, 0.001, percent_cover), 
             group = as.factor(group)) # In order to model w/ a Gamma distribution there can't be any zeros. There is only one zero in the dataset (epiSI at CARP2), so I added 0.001 to the zero value.


# Build the model and run it
mod1 <- rstanarm::stan_glmer(percent_cover ~ (scale(perturbations) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + scale(timesincelastabsent) + scale(prop_zero))*group + (1|site), df, family = Gamma(link = "log"))

# Run model diagnostics and check model priors
rstanarm::prior_summary(mod1) # Used default priors. Rstanarm adjusts priors to make them weakly informative

launch_shinystan(mod1, ppd = FALSE) # Check model convergence, Rhat, and chain mixing

print(summary(mod1), digits = 3) # Print summary of the model

#--------------------------
## Coefficient plot
#--------------------------

# Build out the coefficient plot, with some very inefficient code to deal with the variable names and variable name order in the figure
mod1 %>%
  gather_draws(`scale(perturbations)`, `scale(sand_cover)`, `scale(total_urc_biomass)`, `scale(current_kelpbiomass)`, `scale(timesincelastabsent)`, `scale(prop_zero)`, `scale(perturbations):groupua`, `scale(sand_cover):groupua`, `scale(total_urc_biomass):groupua`, `scale(current_kelpbiomass):groupua`, `scale(timesincelastabsent):groupua`, `scale(prop_zero):groupua`) %>%
  separate(.variable, into = c(".variable", "group"), sep = ":") %>%
  replace_na(list(group = "Sessile\ninvertebrates")) %>%
  mutate(group = ifelse(group == "groupua", "Understory\nalgae", group)) %>%
  group_by(.variable, group) %>%
  median_qi(.width = c(0.75, 0.95)) %>%
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
  labs(x = "Effect on percent_cover", y = "", color = "")+
  theme_bd()+
  theme(legend.position = c(0.8,0.2), legend.background = element_blank())

ggsave("figures/coef_plot-percover.png", width = 7, height = 6)


#--------------------------
## Frequentist approach
#--------------------------

# Check that the model would look similar if fit with a frequentist approach

mod1.lmer <- glmmTMB::glmmTMB(percent_cover ~ (scale(perturbations) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + scale(timesincelastabsent) + scale(prop_zero))*group + (1|site), df, family = Gamma(link = "log"))

summary(mod1.lmer) # In general, all qualitative interpretations will be the same, which makes sense because of the use of weakly informative priors. 

#------------------------------------
## Components of multipanel plot
#------------------------------------

# Read in the theoretical predictions and clean up the data
dist <- read.csv("data/intermediary/persistence_metrics.csv")

rd_ts  <- read.csv("data/Large_files/Raine_modeloutput/FullModelRunsAllSites.csv") %>%
  pivot_longer(cols = ABUR1:WOOD3, names_to = "id", values_to = "percent_cover") %>%
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


#----------------------------------
## Time Since last absent
#----------------------------------

p1 <- ggplot(rd_predictions, aes(x = timesincelastabsent, y = percent_cover))+
  geom_point(aes(color = group))+
  geom_smooth(method = "lm", aes(color = group), se = F)+
  scale_color_manual(values = c(ua_col, si_col))+
  labs(x = "", y = "Theoretical prediction\nrelative abundance", color = "")+
  theme_bd()+
  theme(legend.position = "top")

lm1 <- lm(percent_cover ~ timesincelastabsent*group, rd_predictions)
summary(lm1)

p2 <- df %>%
  filter(group == "ua") %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(timesincelastabsent = modelr::seq_range(timesincelastabsent, n = 100), .model = mod1) %>%
  add_predicted_draws(mod1, re_formula = NA) %>%
  ggplot(aes(x = timesincelastabsent, y = percent_cover)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray80", color = ua_col) +
  geom_point(data = filter(df, group == "ua"), aes(x = timesincelastabsent, y = percent_cover), color = ua_col)+
  theme_bd()+
  labs(x = "", y = "Understory algae percent cover")

p3 <- df %>%
  filter(group == "epiSI") %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(timesincelastabsent = modelr::seq_range(timesincelastabsent, n = 100), .model = mod1) %>%
  add_predicted_draws(mod1, re_formula = NA) %>%
  ggplot(aes(x = timesincelastabsent, y = percent_cover)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray80", color = si_col) +
  geom_point(data = filter(df, group == "epiSI"), aes(x = timesincelastabsent, y = percent_cover), color = si_col)+
  theme_bd()+
  labs(x = "Time since last kelp canopy extinction", y = "Sessile invertebrate percent cover")

cowplot::plot_grid(p1, p2, p3, nrow = 3)

ggsave("figures/time_since-percover.png", width = 6, height = 12 )



#----------------------------------
## Proportion of time no kelp
#----------------------------------


p4 <- ggplot(rd_predictions, aes(x = prop_zero, y = percent_cover))+
  geom_point(aes(color = group))+
  geom_smooth(method = "lm", aes(color = group), se = F)+
  scale_color_manual(values = c(ua_col, si_col))+
  labs(x = "", y = "Theoretical prediction\nrelative abundance", color = "")+
  theme_bd()+
  theme(legend.position = "top")

lm1 <- lm(percent_cover ~ prop_zero*group, rd_predictions)
summary(lm1)

p5 <- df %>%
  filter(group == "ua") %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(prop_zero = modelr::seq_range(prop_zero, n = 100), .model = mod1) %>%
  add_predicted_draws(mod1, re_formula = NA) %>%
  ggplot(aes(x = prop_zero, y = percent_cover)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray80", color = ua_col) +
  geom_point(data = filter(df, group == "ua"), aes(x = prop_zero, y = percent_cover), color = ua_col)+
  theme_bd()+
  labs(x = "", y = "Understory algae percent cover")

p6 <- df %>%
  filter(group == "epiSI") %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(prop_zero = modelr::seq_range(prop_zero, n = 100), .model = mod1) %>%
  add_predicted_draws(mod1, re_formula = NA) %>%
  ggplot(aes(x = prop_zero, y = percent_cover)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray80", color = si_col) +
  geom_point(data = filter(df, group == "epiSI"), aes(x = prop_zero, y = percent_cover), color = si_col)+
  theme_bd()+
  labs(x = "Proportion of time kelp was absent", y = "Sessile invertebrate percent cover")

cowplot::plot_grid(p4, p5, p6, nrow = 3)

ggsave("figures/prop_zero-percover.png", width = 6, height = 12 )


#----------------------------------
## Number of perturbations
#----------------------------------


p7 <- ggplot(rd_predictions, aes(x = perturbations, y = percent_cover))+
  geom_point(aes(color = group))+
  geom_smooth(method = "lm", aes(color = group), se = F)+
  scale_color_manual(values = c(ua_col, si_col))+
  labs(x = "", y = "Theoretical prediction\nrelative abundance", color = "")+
  theme_bd()+
  theme(legend.position = "top")

lm1 <- lm(percent_cover ~ perturbations*group, rd_predictions)
summary(lm1)

p8 <- df %>%
  filter(group == "ua") %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(perturbations = modelr::seq_range(perturbations, n = 100), .model = mod1) %>%
  add_predicted_draws(mod1, re_formula = NA) %>%
  ggplot(aes(x = perturbations, y = percent_cover)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray80", color = ua_col) +
  geom_point(data = filter(df, group == "ua"), aes(x = perturbations, y = percent_cover), color = ua_col)+
  theme_bd()+
  labs(x = "", y = "Understory algae percent cover")

p9 <- df %>%
  filter(group == "epiSI") %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(perturbations = modelr::seq_range(perturbations, n = 100), .model = mod1) %>%
  add_predicted_draws(mod1, re_formula = NA) %>%
  ggplot(aes(x = perturbations, y = percent_cover)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray80", color = si_col) +
  geom_point(data = filter(df, group == "epiSI"), aes(x = perturbations, y = percent_cover), color = si_col)+
  theme_bd()+
  labs(x = "Number of disturbances", y = "Sessile invertebrate percent cover")

cowplot::plot_grid(p7, p8, p9, nrow = 3)

ggsave("figures/perturbations-percover.png", width = 6, height = 12 )

#----------------------------------
## Six panel version
#----------------------------------

cowplot::plot_grid(p1, p2, p3, p7+labs(y = ""), p8+labs(y = ""), p9+labs(y = ""), nrow = 3, ncol = 2, byrow = F, labels = "auto")

ggsave("figures/six_panel_fig3-percover.png", width = 12, height = 12 )








