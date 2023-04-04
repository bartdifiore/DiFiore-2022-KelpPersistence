
rd_predictions <- rd_predictions %>%
  mutate(group = case_when(group == "invert" ~ "epiSI", 
                           group == "algae" ~ "ua"))

mod2 <- rstanarm::stan_glmer(percent_cover ~ (scale(perturbations) + scale(timesincelastabsent) + scale(prop_zero))*group + (1|site), rd_predictions, family = Gamma(link = "log"))

summary(mod2, digits = 3)

mod2 %>%
  gather_draws(`scale(perturbations)`, `scale(timesincelastabsent)`, `scale(prop_zero)`, `scale(perturbations):groupua`, `scale(timesincelastabsent):groupua`, `scale(prop_zero):groupua`) %>%
  separate(.variable, into = c(".variable", "group"), sep = ":") %>%
  replace_na(list(group = "Sessile\ninvertebrates")) %>%
  mutate(group = ifelse(group == "groupua", "Understory\nalgae", group)) %>%
  group_by(.variable, group) %>%
  median_qi(.width = c(0.75, 0.95)) %>%
  # mutate(predictors = forcats::fct_rev(case_when(.variable == "scale(current_kelpbiomass)" ~ "D Current kelp\nbiomass",
  #                                                .variable == "scale(perturbations)" ~ "B Perterbations",
  #                                                .variable == "scale(prop_zero)" ~ "A Proportion time\nno kelp",
  #                                                .variable == "scale(sand_cover)" ~ "E Sand cover",
  #                                                .variable == "scale(timesincelastabsent)" ~ "C Time since\nlast absent",
  #                                                .variable == "scale(total_urc_biomass)" ~ "F Urchin biomass"
  # ))) %>%
  ggplot(aes(y = .variable, x = .value))+
  geom_point(aes(color = group), position = position_dodge(width = 0.5))+
  scale_color_manual(values = c(si_col, ua_col))+
  geom_vline(xintercept = 0, lty = 4)+
  # scale_y_discrete(labels = rev(c("Proportion time\nno kelp", "Perturbations", "Time since\nlast absent", "Current kelp\nbiomass", "Sand cover", "Urchin biomass")))+
  labs(x = "Effect on percent cover", y = "", color = "")+
  theme_bd()+
  theme(legend.position = c(0.8,0.2), legend.background = element_blank())
