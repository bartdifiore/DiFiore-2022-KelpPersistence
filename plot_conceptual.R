
source("code/1_setup.R")
source("code/4_mergepersistence-species.R")

ua_col = "chocolate1"
si_col = "royalblue1"
kelp_col = "darkgreen"

dist <- read.csv("data/intermediary/persistence_metrics.csv")

df <- read.csv("data/Large_files/Raine_modeloutput/FullModelRunsAllSites.csv") %>%
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

filt <- df %>% 
  filter(perturbations == max(perturbations, na.rm= T) | perturbations == min(perturbations, na.rm = T)) %>%
  distinct(id, perturbations)


temp2 <- temp2 %>% left_join(select(dist, id, perturbations))

ggplot(df, aes(x = year.quarter, y = biomass))+
  geom_line(aes(color = group))+
  scale_color_manual(values = c(ua_col, si_col, kelp_col))+
  geom_vline(data = temp2, aes(xintercept = year.quarter[perturbation_events == 1], group = forcats::fct_reorder(id, perturbations)), color = "red", lty = 4)+
  facet_wrap(~forcats::fct_reorder(id, perturbations))+
  theme_bd()

ggsave("figures/theoretical_timeseries_wperturbations.png", width = 20, height = 20)


df %>%
  filter(id %in% c("ABUR1", "REFU1")) %>%
  ggplot(aes(x = year.quarter, y = biomass))+
  geom_line(aes(color = group))+
  geom_line(aes(color = group))+
  scale_color_manual(values = c(ua_col, si_col, kelp_col))+
  facet_wrap(~forcats::fct_reorder(id, perturbations))+
  theme_bd()

ggsave("figures/forconceptual.svg", width = 9.7, height = 14/3)


# individual plot for each panel

ids <- unique(df$id)

for(i in 1:length(ids)){
  temp_plot <- df %>% 
    filter(id == ids[i]) %>%
  ggplot(aes(x = year.quarter, y = biomass))+
    geom_line(aes(color = group), show.legend = F)+
    scale_color_manual(values = c(ua_col, si_col, kelp_col))+
    scale_x_continuous(breaks = c(2010, 2014, 2018))+
    labs(x = "", y = "")+
    ggtitle(ids[i])+
    theme_bd()
  
  ggsave(temp_plot, file = paste0("figures/panels_plots/plot_", ids[i], ".svg"), width = 3, height = 2, units = "in")
}











