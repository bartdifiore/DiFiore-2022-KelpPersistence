source("code/1_setup.R")
source("code/theme.R")

ua_col = "chocolate1"
si_col = "royalblue1"
kelp_col = "darkgreen"

df <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv", na.strings = "-99999", stringsAsFactors = F)

# Find top 10 by biomass across all sites

filt <- df %>% 
  group_by(sp_code, group, scientific_name) %>%
  summarize(total_biomass = sum(dry_gm2, na.rm = T)) %>%
  ungroup() %>%
  group_by(group) %>% 
  slice_max(n = 10, total_biomass) %>%
  filter(group %in% c("epiSI", "ua")) %>%
  select(sp_code, total_biomass, scientific_name)

filt_vec <- as.vector(filt$sp_code)

write.csv(filt, "data/intermediary/toptenbybiomass.csv", row.names = F)


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

sp <- df %>% 
  filter(sp_code %in% filt_vec) %>%
  left_join(urc) %>%
  left_join(sand) %>%
  left_join(kelp) %>%
  mutate(id = paste(site, transect, sep = "_"))

fit_rstanarm <- function(sp_code){
  fit_mod <- rstanarm::stan_glmer(I(dry_gm2+0.01) ~ scale(perturbations) + scale(sand_cover) + scale(total_urc_biomass) + scale(current_kelpbiomass) + scale(timesincelastextinct) + scale(prop_zero) + (1|site), sp[sp$sp_code == sp_code, ], family = Gamma(link = "log"))
  
  posterior_summary <- fit_mod %>%
    gather_draws(`scale(perturbations)`, `scale(timesincelastextinct)`, `scale(prop_zero)`) %>%
    median_qi(.width = c(0.75, 0.95)) %>%
    mutate(sp_code = filt[i])
  
  return(posterior_summary)
  
}


out <- list()
for(i in 1:length(filt_vec)){
  out[[i]] <- fit_rstanarm(sp_code = filt_vec[i])  
}

mod_sum <- do.call(rbind, out) %>%
  mutate(predictors = case_when(.variable == "scale(perturbations)" ~ "Perterbations",
                                .variable == "scale(prop_zero)" ~ "Proportion time\nno kelp",
                                .variable == "scale(timesincelastextinct)" ~ "Time since\nlast extinction"
  ))

mod_sum %>%
  filter(predictors == "Time since\nlast extinction") %>%
  ggplot(aes(y = sp_code, x = .value))+
  geom_pointinterval(aes(xmin = .lower, xmax = .upper, color = predictors), position = position_dodge(width = 0.5))+
  geom_vline(xintercept = 0, lty = 3)+
  theme_classic()





# Just based on correlation

# Prop zero
out <- vector()
for(i in 1:length(filt_vec)){
  out[i] <- cor.test(sp$dry_gm2[sp$sp_code == filt_vec[i]], sp$prop_zero[sp$sp_code == filt_vec[i]], method = "spearman", na.rm = T)$estimate
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
  out[i] <- cor.test(sp$dry_gm2[sp$sp_code == filt_vec[i]], sp$perturbations[sp$sp_code == filt_vec[i]], method = "spearman", na.rm = T)$estimate
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
  out[i] <- cor.test(sp$dry_gm2[sp$sp_code == filt_vec[i]], sp$timesincelastabsent[sp$sp_code == filt_vec[i]], method = "spearman", na.rm = T)$estimate
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
ggsave("figures/species_correlations.svg")


