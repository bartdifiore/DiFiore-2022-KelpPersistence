source("code/1_setup.R")

df <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv") %>% 
  mutate(id = paste(site, transect, sep = "_"), .before = sp_code)

env <- df %>%
  select(c(id:transect, mean.canopy_previousyear:perturbations))

wide <- df %>% 
  filter(sp_code != "MAPY") %>% #filter out kelp
  filter(group %in% c("epiSI", "ua")) %>%
  select(id, site, transect, sp_code, survey,  wm_gm2) %>%
  pivot_wider(names_from = sp_code, values_from = wm_gm2) %>% 
  column_to_rownames("id")

mat <- as.matrix(wide[, -c(1:3)])

ord <- metaMDS(mat, dist = "bray", trymax = 100, k = 3) 
ord
stressplot(ord)
plot(ord)

plot_df <- scores(ord, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("id") %>%
  left_join(env)

ggplot(plot_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(linetype = 2, size = 1) +
  labs(title = "NMDS")











