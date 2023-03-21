source("code/1_setup.R")
source("code/theme.R")

ua_col = "chocolate1"
si_col = "royalblue1"
kelp_col = "darkgreen"
  
bb_upc <- read.csv("data/raw/BartsBenthic_Cover.csv", na.strings = "-99999") %>% 
  janitor::clean_names()

sbc_upc <- read.csv("data/raw/Annual_Cover_All_Years_20220809.csv", na.strings = "-99999") %>%
  janitor::clean_names() %>%
  as_tibble() %>%
  filter(year == 2018)

df <- bind_rows(sbc_upc, bb_upc) %>%
  mutate(sp_code = ifelse(sp_code == "NA", "NA_sp", sp_code)) %>%
  filter(!sp_code %in% c("PACA", "CHOV", "PHTO", "ZOMA", "DMH", "MH")) %>% # Filter out seagrasses and boring clams and dead macrocystis holdfasts and living macrocystis holdfasts
  mutate(group = case_when(group == "INVERT" ~ "epiSI", 
                           group == "ALGAE" ~ "ua")) %>%
  group_by(site, transect, sp_code, scientific_name, common_name, taxon_genus, group) %>%
  summarize(percent_cover = mean(percent_cover, na.rm = T)) %>%# average across the different quads on a transect at a site for each species
  mutate(id = paste(site, transect, sep = ""), .before = sp_code)

# Find top 10 by biomass across all sites

filt <- df %>% 
  group_by(sp_code, group, scientific_name) %>%
  summarize(total_percover = sum(percent_cover, na.rm = T)) %>%
  ungroup() %>%
  group_by(group) %>% 
  slice_max(n = 10, total_percover) %>%
  filter(group %in% c("epiSI", "ua")) %>%
  select(sp_code, total_percover, scientific_name)

write.csv(filt, "data/intermediary/toptenbypercover.csv", row.names = F)

filt_vec <- as.vector(filt$sp_code)


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

persistence <- read.csv("data/intermediary/persistence_metrics.csv")

sp <- df %>% 
  filter(sp_code %in% filt_vec) %>%
  left_join(urc) %>%
  left_join(sand) %>%
  left_join(kelp) %>%
  left_join(persistence) %>%
  mutate(id = paste(site, transect, sep = "_"))



# Just based on correlation

# Prop zero
out <- vector()
for(i in 1:length(filt_vec)){
  out[i] <- cor.test(sp$percent_cover[sp$sp_code == filt_vec[i]], sp$prop_zero[sp$sp_code == filt_vec[i]], method = "spearman", na.rm = T)$estimate
}
correlations <- data.frame(filt, correlation = out)

p1 <- correlations %>% 
  filter(group == "ua") %>%
  ggplot(aes(y = forcats::fct_reorder(sp_code, correlation), x = correlation))+
  geom_bar(stat = "identity", fill = ua_col)+
  labs(y = "Species code", x = expression(rho))+
  theme_bd()

p2 <- correlations %>% 
  filter(group == "epiSI") %>%
  ggplot(aes(y = forcats::fct_reorder(sp_code, correlation), x = correlation))+
  geom_bar(stat = "identity", fill = si_col)+
  labs(y = "", x = expression(rho))+
  theme_bd()

panel_1 <- cowplot::plot_grid(p1, p2)


# Perturbations

out <- vector()
for(i in 1:length(filt_vec)){
  out[i] <- cor.test(sp$percent_cover[sp$sp_code == filt_vec[i]], sp$perturbations[sp$sp_code == filt_vec[i]], method = "spearman", na.rm = T)$estimate
}
correlations <- data.frame(filt, correlation = out)

p1 <- correlations %>% 
  filter(group == "ua") %>%
  ggplot(aes(y = forcats::fct_reorder(sp_code, correlation), x = correlation))+
  geom_bar(stat = "identity", fill = ua_col)+
  labs(y = "Species code", x = expression(rho))+
  theme_bd()

p2 <- correlations %>% 
  filter(group == "epiSI") %>%
  ggplot(aes(y = forcats::fct_reorder(sp_code, correlation), x = correlation))+
  geom_bar(stat = "identity", fill = si_col)+
  labs(y = "", x = expression(rho))+
  theme_bd()

panel_2 <- cowplot::plot_grid(p1, p2)



# Time since
out <- vector()
for(i in 1:length(filt_vec)){
  out[i] <- cor.test(sp$percent_cover[sp$sp_code == filt_vec[i]], sp$timesincelastabsent[sp$sp_code == filt_vec[i]], method = "spearman", na.rm = T)$estimate
}
correlations <- data.frame(filt, correlation = out)

p1 <- correlations %>% 
  filter(group == "ua") %>%
  ggplot(aes(y = forcats::fct_reorder(sp_code, correlation), x = correlation))+
  geom_bar(stat = "identity", fill = ua_col)+
  labs(y = "Species code", x = expression(rho))+
  theme_bd()

p2 <- correlations %>% 
  filter(group == "epiSI") %>%
  ggplot(aes(y = forcats::fct_reorder(sp_code, correlation), x = correlation))+
  geom_bar(stat = "identity", fill = si_col)+
  labs(y = "", x = expression(rho))+
  theme_bd()

panel_3 <- cowplot::plot_grid(p1, p2)

cowplot::plot_grid(NULL, panel_1, panel_2, panel_3, nrow = 2, ncol = 2)
ggsave("figures/species_correlations-percover.svg")
    
    