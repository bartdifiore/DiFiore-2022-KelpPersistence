

##################################
##calculate invert biomass#########
##################################

#################################################
#read algae biomass relationship
biom <- read.csv("data/SBC_masterfiles/Invertebrate_Biomass_Relationships.csv",stringsAsFactors = F,na.strings ="-99999")

proxi <-read.csv("data/SBC_masterfiles/Invertebrate_Proxy_List.csv",stringsAsFactors = F,na.strings ="-99999")

biom1 <- biom %>%
  rename(BIOMASS_SPP=SP_CODE)%>%
  right_join(proxi,by="BIOMASS_SPP") %>%
  select(SP_CODE,TYPE,independent_variable,a,b,slope,smearing_estimate,DM_WM,SFDM_WM,AFDM_WM)

###################################################
# read cover data and select invert
cover <- df.upc

cover1<- cover %>%
  select(YEAR,MONTH,DATE,SITE,TRANSECT,QUAD,SIDE,SP_CODE,PERCENT_COVER,GROUP,SURVEY) %>%
  filter(GROUP=="INVERT") %>%
  left_join(biom1,by="SP_CODE")  %>%
  filter(TYPE=="cover"|is.na(TYPE)) %>% # in case there is a new spp that do not have the coefficient. 
  mutate(WM_GM2 = PERCENT_COVER*slope)


######################################################################## 
#READ QUAD SWATH DATA
quad <- comb_qdsw

quad1 <- quad %>%
  select(YEAR,MONTH,DATE,SITE,TRANSECT, QUAD,SIDE,SP_CODE,COUNT,SIZE,AREA,GROUP,SURVEY) %>%
  filter(GROUP=="INVERT") %>%
  mutate(DENSITY=COUNT/AREA) %>%
  left_join(biom1,by="SP_CODE") %>% #note, SABW should be calcuated in cover
  filter(TYPE!="cover") 


# calculating average size to back fill the NA. In the future, can consider using weighted average size. 
quad2 <- quad1 %>%
  group_by(SITE,TRANSECT,SP_CODE) %>% #SIZE AVERAGE OVER TRANSECT
  mutate(SIZE= replace(SIZE, is.na(SIZE)&COUNT>0&SP_CODE!="DIOR", weighted.mean(SIZE,COUNT, na.rm=TRUE))) %>%
  group_by(SITE,SP_CODE) %>% #SIZE AVERAGE OVER SITE
  mutate(SIZE= replace(SIZE, is.na(SIZE)&COUNT>0&SP_CODE!="DIOR", weighted.mean(SIZE,COUNT, na.rm=TRUE))) %>%
  group_by(SP_CODE) %>% # SIZE AVERAGE OVER POPULATION
  mutate(SIZE= replace(SIZE, is.na(SIZE)&COUNT>0&SP_CODE!="DIOR", weighted.mean(SIZE,COUNT, na.rm=TRUE))) %>%
  ungroup() %>%
  #calculate WET MASS
  mutate(SIZE_C=10*SIZE,  #conver size from cm to mm
         WM_GM2 =if_else(TYPE=="density",DENSITY * slope,if_else(TYPE=="length",a*(SIZE_C^b)*smearing_estimate*DENSITY,as.numeric(NA))),
         WM_GM2 =if_else(COUNT==0&!is.na(COUNT),0,WM_GM2))


###################################################
#combine all the datasets

comb <- bind_rows(cover1,quad2) %>%
  select(YEAR,MONTH,DATE,SITE,TRANSECT,QUAD,SIDE,SP_CODE,PERCENT_COVER,DENSITY,SURVEY,WM_GM2,DM_WM,SFDM_WM,AFDM_WM)

########################################################################
#calculate at the section level
# comb1 <- comb %>%
#   mutate(QUAD=if_else(QUAD <=20,20,40)) %>%
#   group_by(YEAR,MONTH,SITE,TRANSECT,QUAD,SIDE,SP_CODE,SURVEY,DM_WM,SFDM_WM,AFDM_WM) %>%
#   summarise(DATE=max(DATE),PERCENT_COVER=mean(PERCENT_COVER),DENSITY=mean(DENSITY),WM_GM2=mean(WM_GM2)) %>% #we want to use the UPC dates. Which happens to be minimum date in the datasets. 
#   ungroup() 
# 
# colname <- c("YEAR","MONTH","DATE","SITE","TRANSECT","TREATMENT","QUAD","SIDE","SP_CODE","PERCENT_COVER","DENSITY","WM_GM2",
#              "DRY_GM2","SFDM","AFDM","GROUP","SCIENTIFIC_NAME","COMMON_NAME","MOBILITY","GROWTH_MORPH")
# 
# ########################################################################

comb1 <- comb %>%
  group_by(YEAR,SITE,TRANSECT,SP_CODE,SURVEY,DM_WM,SFDM_WM,AFDM_WM) %>%
  summarise(MONTH=min(MONTH),DATE=min(DATE),PERCENT_COVER=mean(PERCENT_COVER),DENSITY=mean(DENSITY),WM_GM2=mean(WM_GM2)) %>% #we want to use the UPC dates. Which happens to be minimum date in the datasets. 
  ungroup() 

#calculate dry, sfdm, and afdm.
comb2 <- comb1 %>%
  mutate(DRY_GM2=WM_GM2*DM_WM/100,
         SFDM=WM_GM2*SFDM_WM/100,
         AFDM=WM_GM2*AFDM_WM/100,
         DRY_GM2=if_else(WM_GM2==0&!is.na(WM_GM2),as.numeric(0),DRY_GM2),
         SFDM=if_else(WM_GM2==0&!is.na(WM_GM2),as.numeric(0),SFDM), 
         AFDM=if_else(WM_GM2==0&!is.na(WM_GM2),as.numeric(0),AFDM)
  ) %>%
  left_join(spp,by=c("SP_CODE","SURVEY"))


colname <- c("YEAR","MONTH","DATE","SITE","TRANSECT","SP_CODE","PERCENT_COVER","DENSITY","WM_GM2",
             "DRY_GM2","SFDM","AFDM","SCIENTIFIC_NAME","COMMON_NAME","TAXON_KINGDOM","TAXON_PHYLUM","TAXON_CLASS",
             "TAXON_ORDER","TAXON_FAMILY","TAXON_GENUS","GROUP","SURVEY","MOBILITY","GROWTH_MORPH")

comb6 <-comb2[,colname]


write.csv(comb6,"data/intermediary/invert_Biomass_at_transect.csv",row.names = F,quote=F)

