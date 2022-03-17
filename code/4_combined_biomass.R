##################################
## Combine all biomass estimates
##################################
source("code/1_setup.R")
source("code/2_Write_cleanedcounts.R")
source("code/3a_algae_biomass.R")
source("code/3b_Fish_biomass.R")
source("code/3c_Invert_biomass.R")

# read species data


spec1 <- spp

#################################################
#read transect level biomass data
algae <- read.csv("data/intermediary/algae_Biomass_at_transect.csv",stringsAsFactors = F)

fish <- read.csv("data/intermediary/Fish_Biomass_at_transect.csv",stringsAsFactors = F)

invert <- read.csv("data/intermediary/invert_Biomass_at_transect.csv",stringsAsFactors = F)

###################################################
#combine all the datasets
colname <- c("YEAR","MONTH","DATE","SITE","TRANSECT","SP_CODE","PERCENT_COVER","DENSITY","WM_GM2",
             "DRY_GM2","SFDM","AFDM")

comb <- bind_rows(algae,fish,invert)[,colname]

comb1 <-comb %>%
  left_join(select(spec1,SP_CODE,ADULT),by="SP_CODE") #make sure there is no duplicate here. 

comb2 <- comb1 %>%
  group_by(YEAR,SITE,TRANSECT)%>%
  summarise(MONTH=min(MONTH),DATE=min(DATE)) %>% #minimum date wa the upc date
  ungroup() %>%
  left_join(select(comb1,-MONTH,-DATE),by=c("YEAR","SITE","TRANSECT")) %>%
  group_by(YEAR,MONTH,DATE,SITE,TRANSECT,ADULT) %>%
  summarise_at(c("PERCENT_COVER","DENSITY","WM_GM2","DRY_GM2","SFDM","AFDM"),sum,na.rm=F) %>%
  ungroup() %>%
  rename(SP_CODE=ADULT) %>%
  left_join(spec1,by="SP_CODE") 
# mutate(
#   TROPHIC_LEVEL = case_when(
#     GUILD=='PRODUCER' ~ "0-PRODUCER",
#     GUILD=='PLANKTIVORE' ~ "1B-PLANKTIVORE",
#     GUILD=='GRAZER'|GUILD=='DEPOSITFEEDER' ~ "1A-GRAZER/DEPOSIT FEEDER",
#     GUILD=='MICROCARNIVORE'|GUILD=='MESOCARNIVORE' ~ "2-CARNIVORE'",
#     GUILD=='MACROCARNIVORE' ~ "3-CARNIVORE",
#     TRUE                    ~  "other"
#   ))
#--#######
#coarse grouping
#comb3<-comb2
comb3<-comb2 %>%
  mutate(
    COARSE_GROUPING = case_when(
       GROUP=='FISH' ~ "FISH",
       GROUP=='INVERT'&MOBILITY=='SESSILE' ~ "SESSILE INVERT",
       GROUP=='INVERT'&MOBILITY=='MOBILE' ~ "MOBILE INVERT",
       SP_CODE=='MAPY' ~ "GIANT KELP",
       GROUP=='ALGAE'&!SP_CODE=="MAPY" ~ "UNDERSTORY ALGAE"
  ))

colname1 <- c("YEAR","MONTH","DATE","SITE","TRANSECT","SP_CODE","PERCENT_COVER","DENSITY","WM_GM2",
              "DRY_GM2","SFDM","AFDM","SCIENTIFIC_NAME","COMMON_NAME","TAXON_KINGDOM","TAXON_PHYLUM","TAXON_CLASS",
              "TAXON_ORDER","TAXON_FAMILY","TAXON_GENUS","GROUP","MOBILITY","GROWTH_MORPH", "COARSE_GROUPING")

comb4 <-as.data.frame(comb3[,colname1])

write.csv(comb4,"data/final_biomass/bartsbenthic_combined_biomass_2018.csv",row.names = F,quote=F)



