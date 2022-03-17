

##################################
##calculate fish biomass#########
##################################

#################################################
#read algae biomass relationship
biom <- read.csv("data/SBC_masterfiles/SBCLTER_Fish_Biomass_relationship.csv",stringsAsFactors = F)

biom1 <- biom %>%
  select(SP_CODE,GROUP,X_Variable,Y_Variable,Slope,Intercept,ConvertoXVar,FISH_GROUP) %>%
  mutate(value=1) #use to check the merging process later

###################################################

fish <- comb_bfcf

fish1 <- fish %>%
  select(YEAR,MONTH,DATE,SITE,TRANSECT,QUAD,SIDE,SP_CODE,SIZE,COUNT,AREA) %>%
  left_join(biom1, by="SP_CODE") %>%
  filter(GROUP=="FISH") #in case there is any wired spp. 

fish2 <- fish1 %>%
  mutate(SIZE=if_else(grepl("mm",X_Variable)&!is.na(SIZE),SIZE*10,SIZE),
         SIZE=if_else(grepl("standard",X_Variable)&!is.na(SIZE),SIZE*ConvertoXVar,SIZE),
         WM=(Intercept*SIZE^Slope)*COUNT,
         WM=if_else(COUNT==0,as.numeric(0),WM),
         WM_GM2 = WM/AREA) %>%
  group_by(YEAR,MONTH,DATE,SITE,TRANSECT,QUAD,SIDE,SP_CODE,AREA,FISH_GROUP) %>%  # sum the biomass within a transect/quad/side
  summarise(WM_GM2=round(sum(WM_GM2),2),SIZE=round(weighted.mean(SIZE,COUNT,na.rm=F),2), COUNT=sum(COUNT)) %>%
  ungroup() %>%
  mutate(DRY_GM2=if_else(FISH_GROUP=="Chondrichthyes",WM_GM2*0.285,WM_GM2*0.255),  #separated into chondrichthyes species (softbone species) or bony fish (others). 
         AFDM =if_else(FISH_GROUP=="Chondrichthyes",as.numeric(NA),WM_GM2*0.205),
         DENSITY=COUNT/AREA) %>%
  select(-FISH_GROUP)

# peace <- fish %>%
#   filter(COUNT>0&is.na(SIZE))

#at the Section level, skip this part. At the transect level run this code.( for the cryptic fish only which has four sections)
fish2_1 <-fish2 %>%
  group_by(YEAR,MONTH,DATE,SITE,TRANSECT,SP_CODE) %>%
  summarise_at(c("DENSITY","WM_GM2","DRY_GM2","AFDM"),mean,na.rm=F) %>%
  ungroup()

##########MERGE WITH SPP INFO. 

# read species data

spec1 <- spp %>%
  filter(GROUP=="FISH")%>%
  filter(!SPP_CHANGE=="REMOVE") 


fish3 <- fish2_1 %>%
  left_join(spec1,by="SP_CODE")


colname <- c("YEAR","MONTH","DATE","SITE","TRANSECT","SP_CODE","DENSITY","WM_GM2",
             "DRY_GM2","AFDM","SCIENTIFIC_NAME","COMMON_NAME","TAXON_KINGDOM","TAXON_PHYLUM","TAXON_CLASS",
             "TAXON_ORDER","TAXON_FAMILY","TAXON_GENUS","GROUP","SURVEY","MOBILITY","GROWTH_MORPH")

fish6 <-fish3[,colname]

write.csv(fish6,"data/intermediary/Fish_Biomass_at_transect.csv",row.names = F,quote=F)







