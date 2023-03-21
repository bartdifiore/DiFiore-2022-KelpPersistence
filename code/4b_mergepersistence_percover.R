source("code/1_setup.R")

persistence <- read.csv("data/intermediary/persistence_metrics.csv")

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
  summarize(percent_cover = mean(percent_cover, na.rm = T)) %>% # average across the different quads on a transect at a site for each species
  group_by(site, transect, group) %>% 
  summarize(percent_cover = sum(percent_cover, na.rm = T)) %>% # add together the percent cover from each species within each group at each transect at each site
  left_join(persistence)

write.csv(df, "data/intermediary/percover_withpersitence.csv", row.names = F)



