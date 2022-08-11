source("code/theme.R")

sim <- read.csv("data/intermediary/stochastic_sims.csv")


sim %>%
  pivot_longer(cols = c(inverts, algae)) %>%
  ggplot(aes(x = n_dist, y = value))+
  geom_point(aes(color = as.factor(t_since)))+
  facet_wrap(~name)


sim %>%
  pivot_longer(cols = c(inverts, algae)) %>%
  ggplot(aes(x = t_since, y = value))+
  geom_point(aes(color = as.factor(n_dist)))+
  facet_wrap(~name)

p1 <- sim %>% 
  select(-algae) %>%
  ggplot(aes(x = n_dist, y = inverts, alpha = t_since))+
  geom_point(color = "#bf7a21", size = 3.5)+
  coord_cartesian(ylim = c(0,1))+
  labs(x = "Disturbance frequency", y = "Relative sessile\ninvertebrate abundance", alpha = "Time since")+
  theme_bd()+
  theme(legend.position = c(0.3, 0.8), legend.direction = "horizontal")+
  guides(alpha = guide_legend(title.position = "top"))


p2 <- df %>%
  filter(group == "epiSI") %>%
  filter(!site %in% c("AHND", "SCDI", "SCTW")) %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(perturbations = modelr::seq_range(perturbations, n = 100), .model = mod1) %>%
  add_predicted_draws(mod1, re_formula = NA) %>%
  ggplot(aes(x = perturbations, y = biomass)) +
  stat_lineribbon(aes(y = .prediction), .width = c(0.95), fill = "gray90", color = "#bf7a21") +
  geom_jitter(data = filter(df, group == "epiSI"), aes(x = perturbations, y = biomass), color = "#bf7a21", alpha = 0.5, size = 3.5)+
  labs(x = "Frequency of significant kelp declines", y = "Sessile invertebrate biomass")+
  theme_bd()

p3 <- sim %>% 
  select(-inverts) %>%
  ggplot(aes(x = t_since, y = algae, alpha = n_dist))+
  geom_point(color = "#0f7046", size = 3.5)+
  coord_cartesian(ylim = c(0,1))+
  labs(x = "Time since disturbance", y = "Relative understory\nalgae abundance", alpha = "Num. dist")+
  theme_bd()+
  theme(legend.position = c(0.3, 0.8), legend.direction = "horizontal")+
  guides(alpha = guide_legend(title.position = "top"))

p4 <- df %>%
  filter(group == "ua") %>%
  filter(!site %in% c("AHND", "SCDI", "SCTW")) %>%
  filter(id != "STEP_1") %>%
  modelr::data_grid(timesincelastextinct = modelr::seq_range(timesincelastextinct, n = 100), .model = mod2) %>%
  add_predicted_draws(mod2, re_formula = NA) %>%
  ggplot(aes(x = timesincelastextinct, y = biomass)) +
  stat_lineribbon(aes(y = .prediction), .width = c( 0.95), fill = "gray90", color = "#0f7046") +
  geom_jitter(data = filter(df, group == "ua"), aes(x = timesincelastextinct, y = biomass), color = "#0f7046", alpha = 0.5, size = 3.5)+
  labs(x = "Time since last kelp canopy extinction", y = "Understory algae biomass")+
  theme_bd()

cowplot::plot_grid(p1,p2,p3,p4, nrow = 2, align = "hv", labels = "AUTO")

ggsave("figures/Raine_ESA.png", device = "png")















