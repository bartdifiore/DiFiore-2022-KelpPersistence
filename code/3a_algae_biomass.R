

##################################
##calculate algae biomass 
##################################


#################################################
#read algae biomass relationship
biom <- read.csv("data/SBC_masterfiles/Understory_Algal_Biomass_Relationships_v2.csv",stringsAsFactors = F,na.strings ="-99999")

biom1 <- select(biom,SP_CODE,Measurement_type,density_slope, Slope_g_dry_m2)


###################################################
# read cover data and filter out algae
cover <- df.upc

cover1 <- cover %>%
  select(YEAR,MONTH,DATE,SITE,TRANSECT, QUAD,SIDE,SP_CODE,PERCENT_COVER,GROUP,SURVEY) %>%
  filter(GROUP=="ALGAE") %>%
  left_join(biom1,by="SP_CODE") %>%
  filter(Measurement_type=="Percent Cover") %>% #Remove holdfast and a few species "SAMU","SAHO","SMJ","SHJ","BLD" that will be calculated in quad/swath
  mutate(DRY_GM2=PERCENT_COVER*Slope_g_dry_m2,
         DENSITY=NA)
#need a estimate for SAGA
###################################################
#quad and swath data
quad <- comb_qdsw

quad1 <- quad %>%
  select(YEAR,MONTH,DATE,SITE,TRANSECT, QUAD,SIDE,SP_CODE,COUNT,SIZE,AREA,GROUP,SURVEY) %>%
  filter(GROUP=="ALGAE") %>%
  mutate(DENSITY=COUNT/AREA) 

#calculate biomass
quad3<-quad1 %>%
  filter(!((SP_CODE=="SAMU"|SP_CODE=="SAHO")&is.na(COUNT))) %>%  # SAMU and SAHO should get NA filled in the later years. 
  left_join(biom1,by="SP_CODE") %>% #we don't filter measurement type because all species in quad/swath are in the calculation. 
  mutate(PERCENT_COVER=NA,
         DRY_GM2=if_else(!is.na(SIZE),Slope_g_dry_m2 * SIZE * DENSITY,density_slope*DENSITY),# if the data has size, use first equation, otherwise use second equation
         DRY_GM2=if_else(Measurement_type=="Count",Slope_g_dry_m2 * DENSITY,DRY_GM2), #if it is blade count, there is no size info, it is multiply by density only 
         DRY_GM2=if_else(SP_CODE=="MPJ"&SIZE<=33&!is.na(SIZE),(0.37* DENSITY),DRY_GM2),#MPJ size structure
         DRY_GM2=if_else(SP_CODE=="MPJ"&SIZE<=66&SIZE>33&!is.na(SIZE),(1.86* DENSITY),DRY_GM2),
         DRY_GM2=if_else(SP_CODE=="MPJ"&SIZE<=99&SIZE>66&!is.na(SIZE),(4.76* DENSITY),DRY_GM2),
         DRY_GM2=if_else(SP_CODE=="PTCA"&SIZE==0&!is.na(SIZE),(1.37* DENSITY),DRY_GM2), #count >0 and size=0 indicates the species is present and has 0 large blades
         DRY_GM2=if_else(COUNT==0&!is.na(COUNT),0,DRY_GM2)) #for species "SAMU" and "SAHO" that count =0 but no need for biomass-density relationship 



###############################################################################################################
#The commented code below are for the species needed to use density_biomass relationship. CYOS is an exception
#This section of the code only run every few years
# bio_sel <- quad2 %>%
#   filter(!(is.na(SIZE)|is.na(COUNT)|X=="Count")) %>%
#   filter(!(SP_CODE=="SAHO"|SP_CODE=="SAMU")) #THIS SPECIES DOESN'T HAVE ANY COUNT THAT IS GREATER THAN 0 WITHOUT SIZE.
# 
# bio_c <- data.frame() # The linear regression has to force to go throught point 0,0. So remove the intercept in the linear model
# for (i in 1:length(unique(bio_sel$SP_CODE))) {
#   est <- summary(lm(DRY_GM2~DENSITY-1,data=subset(bio_sel,SP_CODE==unique(bio_sel$SP_CODE)[i])))$coefficients["DENSITY","Estimate"]
#   bio_c <-rbind(bio_c,data.frame(SP_CODE=unique(bio_sel$SP_CODE)[i],Estimate=est))
# }
# 
# #for CYOS
# cyos_cover<-cover1%>%
#   filter(SP_CODE=="CYOS")%>%
#   select(YEAR,MONTH,SITE,TRANSECT,QUAD,SIDE,DRY_GM2)
# 
# cyos_quad<-quad1 %>%
#   filter(SP_CODE=="CYOS"&YEAR>2010) %>%
#   select(YEAR,MONTH,SITE,TRANSECT,QUAD,SIDE,DENSITY) %>%
#   full_join(cyos_cover,by=c("YEAR","MONTH","SITE","TRANSECT","QUAD","SIDE")) %>% #note, there are missing data for IVEE in 2011, need to fix the quad/swath data before going further.
#   filter(!is.na(DENSITY))
# summary(lm(DRY_GM2~DENSITY-1,data=cyos_quad))
# #going into the biomass relationship table and edit the numbers under the "density_slope" column
###########################################################################################################


#import giant kelp all year
kelp <- df.kp

kelp_rela <- read.csv("data/SBC_masterfiles/Regression parameters to estimate biomass and production from frond density.csv",stringsAsFactors = F,na=-99999) %>%
  filter(model==1&Dependent_variable=="FSC dry in July") %>% #use the July data to calculate the biomass.
  select(model,slope, intercept)
colnames(kelp_rela) <- toupper(colnames(kelp_rela))

kelp1 <- kelp %>%
  group_by(YEAR,MONTH,DATE,SITE,TRANSECT, QUAD,SIDE,AREA) %>%
  summarise(FRONDS=sum(FRONDS)) %>%
  ungroup()%>%
  mutate(PERCENT_COVER=NA,
         SP_CODE="MAPY", #because there is MPS in the data and we want to treat it as MAPY
         DENSITY=FRONDS/AREA,
         SURVEY="KELP", 
         QUAD = as.numeric(QUAD)) %>%
  mutate(DRY_GM2=(kelp_rela$SLOPE*DENSITY)*1000)  # CALCULATE BIOMASS AND CONVERT from kg m-2 to g m-2. 


###################################################
#combine all the datasets

comb <- bind_rows(cover1,quad3,kelp1) %>%
  select(YEAR,MONTH,DATE,SITE,TRANSECT,QUAD,SIDE,SP_CODE,PERCENT_COVER,DENSITY,DRY_GM2,SURVEY)

########################################################################
#calculate at the section level
#QUAD DATA NEEDED TO RECLASSIFY THE QUAD, and average within a swath

# comb1 <- comb %>%
#   mutate(QUAD=if_else(QUAD <=20,20,40)) %>%
#   group_by(YEAR,MONTH,SITE,TRANSECT,QUAD,SIDE,SP_CODE,SURVEY) %>%
#   summarise(DATE=max(DATE),PERCENT_COVER=mean(PERCENT_COVER),DENSITY=mean(DENSITY),DRY_GM2=mean(DRY_GM2)) %>% # WE want the date to be lined up with UPC date, which happen to match the maximum date in a given transect. 
#   ungroup()

# colname <- c("YEAR","MONTH","DATE","SITE","TRANSECT","TREATMENT","QUAD","SIDE","SP_CODE","PERCENT_COVER","DENSITY","WM_GM2",
#              "DRY_GM2","SFDM","AFDM","GROUP","SCIENTIFIC_NAME","COMMON_NAME","MOBILITY","GROWTH_MORPH")

########################################################################
# calculate at the transect level
comb1 <- comb %>%
  group_by(YEAR,SITE,TRANSECT,SP_CODE,SURVEY) %>%
  summarise(MONTH=min(MONTH),DATE=min(DATE),PERCENT_COVER=mean(PERCENT_COVER),DENSITY=mean(DENSITY),DRY_GM2=mean(DRY_GM2)) %>% # WE want the date to be lined up with UPC date, which happen to match the minimum date in a given transect. 
  ungroup()

# Merge with CHN Ratio table
chn <- read.csv("data/SBC_masterfiles/Understory_Algal_CHN_Ratios.csv",stringsAsFactors = F,na=-99999)

chn1 <-chn %>%
  select(SP_CODE,Dry_wet,afd_Dry)

comb2 <-comb1 %>%
  left_join(chn1,by="SP_CODE") %>%
  mutate(WM_GM2=DRY_GM2*(1/Dry_wet),
         AFDM=DRY_GM2*afd_Dry,
         SFDM=DRY_GM2,
         WM_GM2=if_else(DRY_GM2==0&!is.na(DRY_GM2),as.numeric(0),WM_GM2),
         SFDM=if_else(DRY_GM2==0&!is.na(DRY_GM2),as.numeric(0),SFDM), 
         AFDM=if_else(DRY_GM2==0&!is.na(DRY_GM2),as.numeric(0),AFDM)
  ) %>%
  #mutate(DRY_GM2=if_else(SP_CODE %in% c("AMZO","BO","CAL","CO","LI","EC","UEC"), as.numeric(NA), DRY_GM2)) %>%  # note that the biomass for these spp were decalcified dry mass, so using SFDM is more appropriated. We keep here for now.  
  left_join(spp,by=c("SP_CODE","SURVEY"))


colname <- c("YEAR","MONTH","DATE","SITE","TRANSECT","SP_CODE","PERCENT_COVER","DENSITY","WM_GM2",
             "DRY_GM2","SFDM","AFDM","SCIENTIFIC_NAME","COMMON_NAME","TAXON_KINGDOM","TAXON_PHYLUM","TAXON_CLASS",
             "TAXON_ORDER","TAXON_FAMILY","TAXON_GENUS","GROUP","SURVEY","MOBILITY","GROWTH_MORPH")

comb6 <-comb2[,colname]

write.csv(comb6,"data/intermediary/algae_Biomass_at_transect.csv",row.names = F,quote=F)


###########################################

#checking species frequency. Good
# trial <-comb1 %>%
#   group_by(YEAR,MONTH,SITE,TRANSECT,SP_CODE) %>%
#   summarise(freq=n()) %>%
#   ungroup()
