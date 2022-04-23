source("code/1_setup.R")

df <- read.csv("data/intermediary/species_withpersistence-pixelscale.csv")

env <- df %>%
  select(c(year:transect, mean.canopy_previousyear:perturbations))

wide <- df %>% 
  select(year, month, site, transect, sp_code, dry_gm2) %>%
  filter(sp_code != "MAPY") %>% #filter out kelp
  pivot_wider(names_from = sp_code, values_from = dry_gm2)

mat <- as.matrix(wide[, -c(1:4)])

ord <- metaMDS(mat, dist = "bray", trymax = 100, k = 3) 
fit <- envfit(ord, df.ord$dist.all, permutations = 1000)
plot(ord, type = "n")
points(ord, "species", pch = 20, cex = 0.5)
points(ord, "sites", pch = 20, cex = 1.2, col = df.ord$dist.all)
ordiellipse(ord, df.ord$patch)
plot(fit, col = "red", lwd = 1.3)


View(mat)






