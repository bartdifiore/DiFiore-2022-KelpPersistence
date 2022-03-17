######################################################
## Incorporate data from LTER and LTE sites
######################################################

#LTER data
lt <- read.csv("data/final_biomass/Annual_All_Species_Biomass_at_transect.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  select("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "PERCENT_COVER", "DENSITY", "DRY_GM2", "SCIENTIFIC_NAME", "COMMON_NAME", "GROUP", "MOBILITY", "GROWTH_MORPH", "COARSE_GROUPING" ) %>%
  filter(YEAR == 2018)

names(lt) <- tolower(c("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "PERCENT_COVER", "DENSITY", "DRY_GM2", "SCIENTIFIC_NAME", "COMMON_NAME", "GROUP", "MOBILITY", "GROWTH_MORPH", "COARSE_GROUPING"))



#BB data
df <- read.csv("data/final_biomass/bartsbenthic_combined_biomass_2018.csv") %>%
  select("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "PERCENT_COVER", "DENSITY", "DRY_GM2", "SCIENTIFIC_NAME", "COMMON_NAME", "GROUP", "MOBILITY", "GROWTH_MORPH", "COARSE_GROUPING" )

names(df) <- tolower(c("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "PERCENT_COVER", "DENSITY", "DRY_GM2", "SCIENTIFIC_NAME", "COMMON_NAME", "GROUP", "MOBILITY", "GROWTH_MORPH", "COARSE_GROUPING"))

df <- bind_rows(lt, df)

write.csv(df, "data/final_biomass/allsites_combined_biomass_2018.csv", row.names = F)
