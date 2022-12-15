source("code/1_setup.R")
source("code/theme.R")

ua_col = "chocolate1"
si_col = "royalblue1"
kelp_col = "darkgreen"

df <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv", na.strings = "-99999") %>% 
  mutate(id = paste(site, transect, sep = ""), .before = sp_code) %>%
  filter(id != "STEP1")

env <- read.csv("data/intermediary/persistence_metrics.csv")

wide <- df %>% 
  filter(sp_code != "MAPY") %>% #filter out kelp
  filter(group %in% c("epiSI", "ua")) %>%
  select(id, site, transect, sp_code, survey,  dry_gm2) %>%
  pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
  column_to_rownames("id")

mat <- as.matrix(wide[, -c(1:3)])

ord <- metaMDS(mat, dist = "bray", trymax = 100, k = 3) 
ord
stressplot(ord)
plot(ord)

# env_cut <- env %>%
#   mutate(timesince_cut = cut(timesincelastextinct, breaks = c(0,1,5, 10)))

plot_df <- scores(ord, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("id") %>%
  left_join(env) %>%
  mutate(prop_zero_cut = cut(prop_zero, breaks = c(0, 0.4,0.8, 1.1)))


ggplot(plot_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = prop_zero_cut), size = 3, alpha = 0.8) +
  stat_ellipse(aes(color = prop_zero_cut), linetype = 2, size = 1) +
  labs(title = "NMDS")

#----------------
## perMANOVA
#----------------

# Are there differences in the overall (both SI and UA community) benthic community depending on metrics of historic disturbance

env_cut <- env %>%
  mutate(prop_zero_cut = cut(prop_zero, breaks = c(0, 0.4,0.8, 1.1), labels = c("low", "medium", "high")))


benthic_perm <- adonis(mat ~ prop_zero_cut + perturbations + timesincelastabsent, data = env_cut)
benthic_perm

  # seems to be differences by the proportion of time kelp was zero and the number of perturbations


dist_mat <- vegdist(mat,  method = "bray")

beta <- betadisper(dist_mat, env_cut$prop_zero_cut)
temp <- anova(beta)
temp
# the spread is not significant, suggesting that the centroids are in different places




#------------------
## CCA 
#------------------

library(ggvegan)

CCA <- cca(mat ~ perturbations + prop_zero + timesincelastabsent, data = env_cut)
CCA
plot(CCA)

cca_df <- fortify(CCA)

arrows <- filter(cca_df, Score == "biplot") %>%
  mutate(CCA1 = 2*CCA1, 
         CCA2 = 2*CCA2)

cca_df %>%
  filter(Score %in% c("sites", "species")) %>%
  left_join(df %>% distinct(sp_code, group), by = c("Label" = "sp_code")) %>%
  replace_na(list(group = "site")) %>%
ggplot(aes(x = CCA1, y = CCA2))+
  geom_point(aes(color = group), size = 2, alpha = 0.75)+ 
  scale_color_manual(values = c(si_col, "gray", ua_col)) + 
  geom_vline(xintercept = 0, lty = 3)+
  geom_hline(yintercept = 0, lty = 3)+
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
                                               arrow = arrow(length = unit(0.5, "cm")))+
  geom_label(data = arrows, aes(x = CCA1, y = CCA2, label = Label), nudge_x = 0, nudge_y = 0.5)+

  theme_bd()

ggsave("figures/cca_result.png", device = "png", width = 6, height = 4 )
ggsave("figures/cca_result.svg", device = "svg", width = 6, height = 4 )

p1 <- cca_df %>%
  filter(Score %in% c("species")) %>%
  left_join(df %>% distinct(sp_code, group), by = c("Label" = "sp_code")) %>%
  filter(group == "epiSI") %>%
  ggplot(aes(x = forcats::fct_reorder(Label, CCA1), y = CCA1))+
  geom_bar(stat = "identity", color = si_col, fill = si_col)+
  coord_flip()+
  labs(y = "Species code", x= "CCA1")+
  theme_bd()

p2 <- cca_df %>%
  filter(Score %in% c("species")) %>%
  left_join(df %>% distinct(sp_code, group), by = c("Label" = "sp_code")) %>%
  filter(group == "ua") %>%
  ggplot(aes(x = forcats::fct_reorder(Label, CCA1), y = CCA1))+
  geom_bar(stat = "identity", color = ua_col, fill = ua_col)+
  coord_flip()+
  labs(y = "Species code", x= "CCA1")+
  theme_bd()

cowplot::plot_grid(p1, p2)
ggsave("figures/cca1_species_loadings.png", width = 6, height = 12)
ggsave("figures/cca1_species_loadings.svg", width = 6, height = 12)




