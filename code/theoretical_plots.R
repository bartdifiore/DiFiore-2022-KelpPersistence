source("code/theme.R")

ua_col = "chocolate1"
  si_col = "royalblue1"
    kelp_col = "darkgreen"

#------------------------------------------------------------------------------------------------------
## How long does it take for the benthic community to recover as a function of disturbance frequency?
#------------------------------------------------------------------------------------------------------

recovery <- read.csv("data/intermediary/ModelRecovTimes.csv") %>%
  select(-X) %>%
  pivot_longer(cols = rtimes1:rtimes25, names_to = "scouring", values_to = "recovery_time") %>%
  mutate(scouring = case_when(scouring == "rtimes1" ~ 1,
                              scouring == "rtimes25" ~ 0.25,
                              scouring == "rtimes50" ~ 0.50,
                              scouring == "rtimes75" ~ 0.75), 
         recovery_time_years = recovery_time/365)

ggplot(recovery, aes(x = ndist, y = recovery_time_years))+
  geom_line(aes(color = as.factor(scouring)))+
  theme_bd()

recovery_plot <- ggplot(recovery, aes(x = ndist, y = recovery_time_years))+
  geom_line(aes(color = as.factor(scouring)), linewidth = 1.5)+
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 3))+
  labs(x = "Number of consecutive annual disturbances", y = "Benthic community recover time (years)", color = "Scouring level")+
  scale_y_continuous(breaks = seq(0.5, 3.5, by = 0.5))+
  scale_x_continuous(breaks = seq(0, 10, by = 2))+
  scale_color_manual(values = c('#bae4b3','#74c476','#31a354','#006d2c'))+
  theme_bd()+
  theme(legend.position = c(0.8, 0.4))


#------------------------------------------------------------------------------------------------------
## How long does it take for the benthic community to recover as a function of disturbance frequency?
#------------------------------------------------------------------------------------------------------

main_sims <- read.csv("data/Large_files/Raine_modeloutput/MainSims.csv")

# t_since = years since last disturbance, n_dist = total number of disturbances in simulation, inverts = invert abundance (fractional cover) at end of simulation, algae = algal abundance (fraction cover) at end of simulation, kelp = kelp abundance (fraction of carrying capacity) at end of simulation

main_sims %>% 
  select(-kelp) %>%
  pivot_longer(cols = inverts:algae) %>%
  ggplot(aes(x = n_dist, y = value))+
  geom_line(aes(color = name, alpha = as.factor(t_since)))


main_sims %>% 
  select(-kelp) %>%
  mutate(t_since = as.factor(t_since), 
         n_dist = as.factor(n_dist)) %>%
  pivot_longer(cols = inverts:algae) %>%
  ggplot(aes(x = t_since, y = value))+
  geom_line(aes(color = name, alpha = interaction(n_dist, name)), lwd = 1.5)

p1 <- main_sims %>% 
  select(-kelp) %>%
  filter(n_dist == 1 | n_dist == 8) %>%
  mutate(n_dist = as.factor(n_dist)) %>%
  pivot_longer(cols = inverts:algae) %>%
  ggplot(aes(x = t_since, y = value, group = interaction(n_dist, name)))+
  geom_line(aes(color = name, linetype = n_dist), lwd = 1.5)+
  scale_color_manual(values = c(ua_col, si_col))+
  scale_x_continuous(breaks = seq(0, 10, by = 2))+
  labs(x = "Time since the last disturbance (years)", y = "Estimated percent cover", linetype = "Number of\ndisturbances", color = "Guild")+
  theme_bd()+
  theme(legend.position = c(0.8,0.5), legend.background = element_blank())


twopanel <- cowplot::plot_grid(p1, recovery_plot, nrow = 2)

ggsave("figures/theoretical_p2.png", twopanel, width = 5, height = 8)














