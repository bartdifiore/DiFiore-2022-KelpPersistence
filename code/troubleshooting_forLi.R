#LTER data
lt <- read.csv("data/raw/Annual_All_Species_Biomass_at_transect_20210108.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  filter(YEAR == 2018) %>%
  mutate(survey = "core") %>%
  rename_all(tolower)


sp.meta <- lt %>% 
  select(sp_code, scientific_name:coarse_grouping) %>%
  distinct()



# Noncore site data
df <- read.csv("data/raw/BartsBenthic_All_Species_Biomass_at_transect_20220424.csv", stringsAsFactors = F,na.strings ="-99999") %>% 
  mutate(survey = "noncore") %>%
  rename_all(tolower)


sp.bart <- unique(df$sp_code)

sp.lt <- unique(lt$sp_code)


notinbart <- setdiff(sp.lt, sp.bart) # need to zero fill each of these species for each transect
notinlt <- setdiff(sp.bart, sp.lt)

  
dat <- bind_rows(df, lt) %>%
  select(site:sp_code, wm_gm2, survey) %>%
  complete(sp_code, 
           nesting(site, transect, survey), 
           fill = list(wm_gm2 = 0)) %>%
  left_join(sp.meta) %>% 
  mutate(group = case_when(coarse_grouping == "SESSILE INVERT" & !sp_code %in% c("PACA", "CHOV") ~ "epiSI", 
                           coarse_grouping == "GIANT KELP" ~ "kelp",
                           coarse_grouping == "UNDERSTORY ALGAE" ~ "ua", 
                           coarse_grouping == "FISH" ~ "fish", 
                           coarse_grouping == "MOBILE INVERT" ~ "mobile_invert", 
                           coarse_grouping == "SESSILE INVERT" & sp_code %in% c("PACA", "CHOV") ~ "endoSI")) %>%
  select(c(year:common_name, group)) %>%
  left_join(persistence)
  

unique(filter(temp, survey == "noncore")$sp_code)
unique(filter(temp, survey == "core")$sp_code)


  
temp %>% select(site, transect, survey, sp_code, wm_gm2) %>%
  pivot_wider(id_cols = c(site, transect, survey), names_from = sp_code, values_from = wm_gm2) %>%
  arrange(site, transect) %>% View()

  
  
  
  
  
  
  
  
  
  
  
  
  
  d <- tibble(
    site = c("NAPL", "NAPL", "NAPL", "OAKS", "OAKS"),
    transect = c(1, 1, 1, 7, 7),
    sp_code = c("a", "b", "c", "a", "a"),
    value1 = c(1,2,3, 1,1)
  )  

  d %>%
    complete(nesting(site, transect), sp_code, fill = list(value1 = 0))
  

