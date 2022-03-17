###################################################################
## Libraries
###################################################################
library(reshape2)
library(dplyr)
# library(tidyverse)

# library(ggpubr)
# library(rgeos) 
# library(rgdal)
# library(sp)
# library(lme4)
# library(lmerTest)
# library(multcomp)
# library(vegan)


###################################################################
## Build output files from manually cleaned csv's
###################################################################



spp<- read.csv("data/SBC_masterfiles/SBCLTER_Species_List_Master.csv",stringsAsFactors = F,na=c(".","") )%>%
  dplyr::select(-SIZE)%>%
  mutate(START_NA=format(as.Date(START_NA,"%m/%d/%Y"),"%Y-%m-%d"),
         END_NA=format(as.Date(END_NA,"%m/%d/%Y"),"%Y-%m-%d"))

amlw <- data.frame("AMLW", "Patiria minata", "Bat Star", "AMLW", 2.5, 21, "INVERT", "SWATH", "KEEP", "NA", "NA", "NA", "NA", "MESOCARNIVORE", "HARD", "MOBILE", "SOLITARY", 382131, "Species", "Animalia", "Echinodermata", "Asteroidea", "Valvatida", "Asterinidae", "Patiria", "Asterina", "miniata") %>% mutate_if(is.factor, as.character)


names(amlw) <- names(spp)


spp <- spp %>%
  bind_rows(amlw)

############################
## Function to clean all files
###########################

sbc_clean <- function(data, survey){
  if(survey == "bigfish"){
    ################################
    ## Big fish
    ################################
    
    temp1 <- data %>%
      select(-ENTERED_BY,-CHECKED_BY) %>%
      mutate(SP_CODE=trimws(SP_CODE),
             DATE=format(as.Date(DATE,"%m/%d/%Y"),"%Y-%m-%d"),
             SIZE=if_else(COUNT==0&!is.na(COUNT),as.numeric(NA),as.numeric(SIZE)),
             NOTES=gsub(",", ".", NOTES,fixed = TRUE),
             NOTES=if_else(is.na(NOTES),"",NOTES))
    
    nasite <- temp1 %>%  
      filter(is.na(SP_CODE)) %>%
      select(DATE,SITE,TRANSECT) %>%
      mutate(VALUE1="missing")
    
    temp2 <- temp1 %>%
      select(SITE,TRANSECT,DATE,SP_CODE) %>%
      distinct() %>%
      mutate(VALUE=1) %>%
      dcast(DATE+SITE+TRANSECT~SP_CODE,value.var = "VALUE",fill=0) %>%
      melt(id.var=c("SITE","TRANSECT","DATE"),value.name="present",variable.name="SP_CODE") %>%
      mutate(SP_CODE=as.character(SP_CODE)) %>%
      left_join(temp1,by=c("SITE","TRANSECT","DATE","SP_CODE")) %>%
      filter(!(SP_CODE=="NF"|SP_CODE=="NA")) %>% #remove NF after zero fill, NA was generated for the nasite
      left_join(nasite,by=c("DATE","SITE","TRANSECT")) %>%
      mutate(YEAR=as.integer(format(as.Date(DATE,format="%Y-%m-%d"),"%Y")),
             MONTH=as.integer(format(as.Date(DATE,format="%Y-%m-%d"),"%m")),
             NOTES=if_else(is.na(COUNT),"zero filling",NOTES),
             COUNT=if_else(is.na(COUNT)&is.na(VALUE1),as.integer(0),COUNT),
             AREA=80) %>%
      select(-present,-VALUE1)
    
    temp3 <- spp %>%
      filter(SURVEY == "FISH") %>%
      right_join(temp2,by="SP_CODE")
    
    colname<-c("YEAR","MONTH","DATE","SITE", "VIS", "TRANSECT","SP_CODE", "SIZE","COUNT","AREA","NOTES", "SCIENTIFIC_NAME","COMMON_NAME","TAXON_KINGDOM","TAXON_PHYLUM","TAXON_CLASS","TAXON_ORDER","TAXON_FAMILY","TAXON_GENUS","GROUP","SURVEY","MOBILITY","GROWTH_MORPH")
    
    temp3[,colname]
  }
  
  else{
    if(survey == "crypticfish"){
      ####################################
      ## Cryptic fish
      ####################################
      
      temp1 <- data %>%
        select(-ENTERED_BY,-CHECKED_BY) %>%
        mutate(SP_CODE=trimws(SP_CODE),
               DATE=format(as.Date(DATE,"%m/%d/%Y"),"%Y-%m-%d"),
               SECTION=if_else(SECTION=="0-20","20",SECTION),
               SECTION=if_else(SECTION=="21-40","40",SECTION),
               SIDE=if_else(SIDE=="IN","I",SIDE),
               SIDE=if_else(SIDE=="OFF","O",SIDE),
               NOTES=gsub(",", ".", NOTES,fixed = TRUE),
               NOTES=if_else(is.na(NOTES),"",NOTES),
               SIZE=if_else(COUNT==0&!is.na(COUNT),as.numeric(NA),as.numeric(SIZE)),
               AREA=if_else(SIDE=="",40,20)) %>% #CRYPTIC FISH SURVEY THE AREA CHANGED.
        rename(QUAD=SECTION) 
      
      
      # deal with "remove" spp 
      temp2 <- spp %>%
        filter(SURVEY == "CRYPTIC FISH") %>%
        select(SP_CODE,SPP_CHANGE,NA_NOTE) %>%
        right_join(temp1,by="SP_CODE") %>%
        mutate(SP_CODE=if_else(SPP_CHANGE=="REMOVE"&!is.na(SPP_CHANGE),"NF",SP_CODE),
               # the !is.na(spp_change) just in case there is a wrong spp code, we can check it here)
               # this case, there was only one species on the transect, after remove the species, we should fill in NF. 
               SP_CODE=if_else(SPP_CHANGE=="LUMP"&!is.na(SPP_CHANGE),NA_NOTE,SP_CODE)) %>%
        select(-SPP_CHANGE,-NA_NOTE) %>%  #after lump spp, need to take a sum of the count
        group_by(YEAR,MONTH,DATE,SITE,TRANSECT,SIDE,QUAD,SP_CODE,AREA) %>%
        summarise(SIZE=weighted.mean(SIZE,COUNT,na.rm=F), COUNT=sum(COUNT),NOTES=paste(NOTES, collapse="")) %>% #weighted size
        ungroup()
      
      
      #  zero-fill the transect. Check to see if count and sp code is NA. 
      temp3 <-temp2 %>%
        select(SITE,TRANSECT,DATE,QUAD,SIDE,SP_CODE,AREA) %>%
        distinct() %>%
        mutate(VALUE=1) %>%
        dcast(DATE+SITE+TRANSECT+QUAD+SIDE+AREA~SP_CODE,value.var = "VALUE",fill=0) %>%
        melt(id.var=c("SITE","TRANSECT","DATE","QUAD","SIDE","AREA"),value.name="present",variable.name="SP_CODE") %>%
        left_join(select(temp2,-AREA),by=c("SITE","TRANSECT","QUAD","SIDE","DATE","SP_CODE")) %>%
        select(-present) %>%
        filter(!SP_CODE=="NF") %>% #remove NF after zero fill
        mutate(YEAR=as.integer(format(as.Date(DATE,format="%Y-%m-%d"),"%Y")),
               MONTH=as.integer(format(as.Date(DATE,format="%Y-%m-%d"),"%m")),
               NOTES=if_else(is.na(COUNT),"zero filling",NOTES),
               COUNT=if_else(is.na(COUNT),as.integer(0),COUNT))
      
      #assign taxa info
      temp4 <- spp %>%
        filter(SURVEY == "CRYPTIC FISH") %>%
        right_join(temp3,by="SP_CODE") 
      
      colname<-c("YEAR","MONTH","DATE","SITE","TRANSECT","QUAD","SIDE","SP_CODE", "SIZE","COUNT","AREA","NOTES",
                 "SCIENTIFIC_NAME","COMMON_NAME","TAXON_KINGDOM","TAXON_PHYLUM","TAXON_CLASS","TAXON_ORDER",
                 "TAXON_FAMILY","TAXON_GENUS","GROUP","SURVEY","MOBILITY","GROWTH_MORPH")
      
      temp5 <- temp4[,colname]
      
      temp5
      
    }
    else{ if(survey == "swath"){
      
      #################################
      ## Swath
      ################################
      
      temp <- data %>%
        select(-ENTERED_BY,-CHECKED_BY) %>%
        mutate(DATE=format(as.Date(DATE,"%m/%d/%Y"),"%Y-%m-%d"),
               SIDE=trimws(SIDE),
               SP_CODE=trimws(SP_CODE),
               SECTION=SIDE,
               SECTION=if_else(SECTION=="IN 20"|SECTION=="OFF 20","20",SECTION),
               SECTION=if_else(SECTION=="IN 40"|SECTION=="OFF 40","40",SECTION),
               SECTION=as.numeric(SECTION),
               SIDE=if_else(SIDE=="IN 20"|SIDE=="IN 40","I",SIDE),
               SIDE=if_else(SIDE=="OFF 20"|SIDE=="OFF 40","O",SIDE),
               NOTES=gsub(",", ".", NOTES,fixed = TRUE),
               NOTES=if_else(is.na(NOTES),"",NOTES),
               SIZE=if_else(COUNT==0&!is.na(COUNT),as.numeric(NA),as.numeric(SIZE))
        ) %>%
        rename(QUAD=SECTION)
      
      
      # deal with "remove" and "lump" spp 
      temp1 <- spp %>%
        filter(SURVEY == "SWATH") %>%
        select(SP_CODE,SPP_CHANGE,NA_NOTE) %>%
        right_join(temp,by="SP_CODE") %>%
        filter(!(SPP_CHANGE=="REMOVE"&!is.na(SPP_CHANGE))) %>% # the !is.na(spp_change) just in case there is a wrong spp code, we can check it here)
        mutate(SP_CODE=if_else(SPP_CHANGE=="LUMP"&!is.na(SPP_CHANGE),NA_NOTE,SP_CODE)) %>%
        select(-SPP_CHANGE,-NA_NOTE) %>%  #after lump spp, need to take a sum of the count
        group_by(YEAR,MONTH,DATE,SITE,TRANSECT,SIDE,QUAD,SP_CODE) %>%
        summarise(SIZE=weighted.mean(SIZE,COUNT,na.rm=F), COUNT=sum(COUNT),NOTES=paste(NOTES, collapse="")) %>% #weighted size
        ungroup()
      
      #  zero-fill the transect. some spp with NA count such as CYOS (didn't survey)
      
      naspp <- temp1 %>% 
        filter(is.na(COUNT)) %>%
        select(DATE,SITE,TRANSECT,QUAD,SIDE,SP_CODE) %>%
        mutate(VALUE="missing")
      
      
      temp2 <- temp1 %>%
        select(SITE,TRANSECT,DATE,QUAD,SIDE,SP_CODE) %>%
        distinct() %>%
        mutate(VALUE=1) %>%
        dcast(DATE+SITE+TRANSECT+QUAD+SIDE~SP_CODE,value.var = "VALUE",fill=0) %>%
        melt(id.var=c("SITE","TRANSECT","DATE","QUAD","SIDE"),value.name="present",variable.name="SP_CODE") %>%
        mutate(SP_CODE=as.character(SP_CODE)) %>%
        left_join(temp1,by=c("SITE","TRANSECT","QUAD","SIDE","DATE","SP_CODE")) %>%
        filter(!SP_CODE=="NOSP") %>%
        left_join(naspp,by=c("DATE","SITE","TRANSECT","QUAD","SIDE","SP_CODE")) %>%
        mutate(YEAR=as.integer(format(as.Date(DATE,format="%Y-%m-%d"),"%Y")),
               MONTH=as.integer(format(as.Date(DATE,format="%Y-%m-%d"),"%m")),
               NOTES=if_else(is.na(COUNT)&is.na(VALUE),"zero filling",NOTES),
               COUNT=if_else(is.na(COUNT)&is.na(VALUE),as.integer(0),COUNT),
               AREA=20) %>%
        select(-present,-VALUE)
      
      #NA fill some species. 
      temp3 <- spp %>%
        filter(SURVEY == "SWATH") %>%
        filter(SPP_CHANGE=="NA-FILL") %>%
        mutate(START_NA=as.Date(START_NA),END_NA=as.Date(END_NA)) %>%
        select(SP_CODE,SPP_CHANGE,START_NA,END_NA) %>%
        right_join(temp2,by="SP_CODE") %>%
        mutate(COUNT=if_else(SPP_CHANGE=="NA-FILL"&COUNT==0&DATE>START_NA&DATE<END_NA&!is.na(START_NA),as.integer(NA),COUNT)) %>%
        select(-SPP_CHANGE,-START_NA,-END_NA)
      
      
      #assign taxa info
      temp4 <- spp %>%
        filter(SURVEY == "SWATH") %>%
        right_join(temp3,by="SP_CODE")
      
      colname<-c("YEAR","MONTH","DATE","SITE","TRANSECT","QUAD","SIDE","SP_CODE", "SIZE","COUNT","AREA",
                 "SCIENTIFIC_NAME","COMMON_NAME","TAXON_KINGDOM","TAXON_PHYLUM","TAXON_CLASS","TAXON_ORDER",
                 "TAXON_FAMILY","TAXON_GENUS","GROUP","SURVEY","MOBILITY","GROWTH_MORPH")
      
      temp5 <- temp4[,colname]
      temp5
    }
      else{ if(survey == "quad"){
        #################################
        ## Quad
        #################################
        
        temp <- data %>%
          select(-ENTERED_BY,-CHECKED_BY) %>%
          filter(!is.na(SP_CODE)) %>% #REMOVE THE MISSING SHEET SITE. CHECK TO SEE IF THERE IS ANY MISSING SHEET HAPPEN ON THE TRANSECT. 
          mutate(DATE=format(as.Date(DATE,"%m/%d/%Y"),"%Y-%m-%d"),
                 SIDE=trimws(SIDE),
                 SIDE=if_else(SIDE=="ON"&!is.na(SIDE),"I",SIDE),
                 SIDE=if_else(SIDE=="OFF"&!is.na(SIDE),"O",SIDE),
                 SP_CODE=trimws(SP_CODE),
                 NOTES=gsub(",", ".", NOTES,fixed = TRUE),
                 NOTES=if_else(is.na(NOTES),"",NOTES)
          ) 
        
        # deal with "remove" and "lump" spp ; Be cautious about removing a spp and not back fill NA. Check the missing side survey using this dataset
        temp1 <- spp %>%
          filter(SURVEY == "QUAD") %>%
          select(SP_CODE,SPP_CHANGE,NA_NOTE,SURVEY) %>%
          right_join(temp,by="SP_CODE") %>%
          filter(SURVEY=="QUAD") %>% #remove cryptic spp from the list, CHECK ahead of time if there are missing sheet or missing survey for any part of transect. 
          filter(!(SPP_CHANGE=="REMOVE"&!is.na(SPP_CHANGE))) %>% # the !is.na(spp_change) just in case there is a wrong spp code, we can check it here)
          mutate(SP_CODE=if_else(SPP_CHANGE=="LUMP"&!is.na(SPP_CHANGE),NA_NOTE,SP_CODE)) %>%
          select(-SPP_CHANGE,-NA_NOTE) %>%  #after lump spp, need to take a sum of the count
          group_by(YEAR,MONTH,DATE,SITE,TRANSECT,SIDE,QUAD,SP_CODE) %>%
          summarise(SIZE=weighted.mean(SIZE,COUNT,na.rm=F),COUNT=sum(COUNT),NOTES=paste(NOTES, collapse="")) %>%
          ungroup()
        
        
        #assign size for some spp
        temp2 <-temp1 %>%
          mutate(SIZE=if_else(SP_CODE=="MPJ33"&is.na(SIZE),33,SIZE),
                 SIZE=if_else(SP_CODE=="MPJ66"&is.na(SIZE),66,SIZE),
                 SIZE=if_else(SP_CODE=="MPJ99"&is.na(SIZE),99,SIZE),
                 SP_CODE=if_else(SP_CODE=="MPJ33"|SP_CODE=="MPJ66"|SP_CODE=="MPJ99","MPJ",SP_CODE),
                 SIZE=if_else((SP_CODE=="AMS"|SP_CODE=="DIS"|SP_CODE=="PGS"|SP_CODE=="POS"|SP_CODE=="PBS"|SP_CODE=="PHS"|SP_CODE=="SKE"|
                                 SP_CODE=="SPS"|SP_CODE=="SFS"|SP_CODE=="OKS")&is.na(SIZE),2.5,SIZE),
                 SIZE=if_else(SP_CODE=="LIGS"&is.na(SIZE),8.9,SIZE))
        
        temp3 <- temp2 %>% #because we change the MPJ category, we have to sum and count and weighted average the size. 
          group_by(YEAR,MONTH,DATE, SITE, TRANSECT,SIDE,QUAD,SP_CODE) %>%
          summarise(SIZE=weighted.mean(SIZE, COUNT,na.rm=F),COUNT=sum(COUNT),NOTES=paste(NOTES, collapse="")) %>%
          ungroup()
        
        #zero-fill, check to see if there is any record that spp code is NA. If yes, add into nasite piece of the code from the other survey
        naspp <- temp3 %>% 
          filter(is.na(COUNT)) %>%
          select(DATE,SITE,TRANSECT,QUAD,SIDE,SP_CODE) %>%
          mutate(VALUE="missing")
        
        temp4 <-temp3 %>%
          select(SITE,TRANSECT,DATE,QUAD,SIDE,SP_CODE) %>%
          distinct() %>%
          mutate(VALUE=1) %>%
          dcast(DATE+SITE+TRANSECT+QUAD+SIDE~SP_CODE,value.var = "VALUE",fill=0) %>%
          melt(id.var=c("SITE","TRANSECT","DATE","QUAD","SIDE"),value.name="present",variable.name="SP_CODE") %>%
          mutate(SP_CODE=as.character(SP_CODE)) %>%
          left_join(temp3,by=c("SITE","TRANSECT","QUAD","SIDE","DATE","SP_CODE")) %>%
          filter(!SP_CODE=="NOSP") %>%
          left_join(naspp,by=c("DATE","SITE","TRANSECT","QUAD","SIDE","SP_CODE")) %>%
          mutate(YEAR=as.integer(format(as.Date(DATE,format="%Y-%m-%d"),"%Y")),
                 MONTH=as.integer(format(as.Date(DATE,format="%Y-%m-%d"),"%m")),
                 NOTES=if_else(is.na(COUNT)&is.na(VALUE),"zero filling",NOTES),
                 COUNT=if_else(is.na(COUNT)&is.na(VALUE),as.integer(0),COUNT),
                 AREA=1) %>%
          select(-present,-VALUE)
        
        #NA fill some species. 
        temp5 <- spp %>%
          filter(SURVEY == "QUAD") %>%
          filter(SPP_CHANGE=="NA-FILL") %>%
          mutate(START_NA=as.Date(START_NA),END_NA=as.Date(END_NA)) %>%
          select(SP_CODE,SPP_CHANGE,START_NA,END_NA) %>%
          right_join(temp4,by="SP_CODE") %>%
          mutate(COUNT=if_else(SPP_CHANGE=="NA-FILL"&COUNT==0&DATE>START_NA&DATE<END_NA&!is.na(START_NA),as.integer(NA),COUNT)) %>%
          select(-SPP_CHANGE,-START_NA,-END_NA)
        
        #assign taxa info
        temp6 <- spp %>%
          filter(SURVEY == "QUAD") %>%
          right_join(temp5,by="SP_CODE") 
        
        
        colname<-c("YEAR","MONTH","DATE","SITE","TRANSECT","QUAD","SIDE","SP_CODE", "SIZE","COUNT","AREA",
                   "SCIENTIFIC_NAME","COMMON_NAME","TAXON_KINGDOM","TAXON_PHYLUM","TAXON_CLASS","TAXON_ORDER",
                   "TAXON_FAMILY","TAXON_GENUS","GROUP","SURVEY","MOBILITY","GROWTH_MORPH")
        
        temp7 <- temp6[,colname]
        temp7
        
        
      }
        else{ if(survey == "upc"){
          #################################
          ## UPC 
          #################################
          
          #Convert UPC to cover
          #read in UPC data
          
          # assign quad location AND Calculate the total number of observations for each transect, date, quad, side, sp_code
          temp <- data %>%
            filter(!is.na(SP_CODE)) %>%  
            mutate(QUAD=if_else(DISTANCE<=20,20,40), 
                   SP_CODE=trimws(SP_CODE))%>%
            group_by(YEAR,MONTH,DATE,SITE,TRANSECT,QUAD,SIDE,SP_CODE) %>%
            summarise(COUNT=n()) %>%
            ungroup() %>%
            mutate(PERCENT_COVER=(COUNT/20)*100) #one species can't have percent_cover more than 100%
          
          #zero-filling
          temp2 <- temp %>%
            dcast(YEAR+MONTH+DATE+SITE+TRANSECT+QUAD+SIDE~SP_CODE,value.var="PERCENT_COVER",fill=0) %>%
            melt(id.var=c("YEAR","MONTH","DATE","SITE","TRANSECT","QUAD","SIDE"),value.name="PERCENT_COVER",variable.name="SP_CODE") %>%
            filter(!SP_CODE=="NOSP")
          
          temp3 <- spp %>%
            filter(SURVEY == "UPC") %>%
            filter(!SPP_CHANGE=="KEEP") %>%
            mutate(START_NA=as.Date(START_NA),END_NA=as.Date(END_NA)) %>%
            select(SP_CODE,SPP_CHANGE,START_NA,END_NA) %>%
            right_join(temp2,by="SP_CODE") %>%
            mutate(PERCENT_COVER=if_else(SPP_CHANGE=="NA-FILL"&PERCENT_COVER==0&DATE>START_NA&DATE<END_NA&!is.na(START_NA),as.numeric(NA),PERCENT_COVER)) %>%
            select(-SPP_CHANGE,-START_NA,-END_NA)
          
          
          #assign taxa info
          temp4 <- spp %>%
            filter(SURVEY == "UPC") %>%
            right_join(temp3,by="SP_CODE") 
          
          colname<-c("YEAR","MONTH","DATE","SITE","TRANSECT","QUAD","SIDE","SP_CODE","PERCENT_COVER",
                     "SCIENTIFIC_NAME","COMMON_NAME","TAXON_KINGDOM","TAXON_PHYLUM","TAXON_CLASS","TAXON_ORDER",
                     "TAXON_FAMILY","TAXON_GENUS","GROUP","SURVEY","MOBILITY","GROWTH_MORPH")
          
          temp5 <- temp4[,colname]
          temp5
          
          
        }
          
          else{ if(survey == "substrate"){
            #################################
            ## Substrate
            #################################
            
            #substrate percentage calculation
            subs <- data %>%
              filter(POSITION==1) %>% #There might be more then one species were surveyed. so only take a first position for the substrate. 
              mutate(QUAD=if_else(DISTANCE<=20,20,40),
                     SUBSTRTE=trimws(SUBSTRTE))%>%  
              group_by(YEAR,MONTH,DATE,SITE,TRANSECT,QUAD,SIDE,SUBSTRTE) %>%
              summarise(COUNT=n()) %>%
              ungroup() %>%
              rename(SUBSTRATE_TYPE=SUBSTRTE) %>%
              mutate(PERCENT_COVER=(COUNT/20)*100) %>%
              select(-COUNT) %>%
              dcast(YEAR+MONTH+DATE+SITE+TRANSECT+QUAD+SIDE~SUBSTRATE_TYPE,value.var = "PERCENT_COVER",fill=0) %>% #zero filled the substrate
              melt(id.var=c("YEAR","MONTH","DATE","SITE","TRANSECT","QUAD","SIDE"),value.name = "PERCENT_COVER",variable.name="SUBSTRATE_TYPE")
            
            subs
          }
            else{ if(survey == "sand"){
              ###################################
              ## Sand
              ##################################
              # SAND DEPTH PART
              depth <- data %>%
                filter(POSITION==1) %>% #There might be more then one species were surveyed. so only take a first position for the substrate. 
                select(YEAR,MONTH,DATE,SITE,TRANSECT,DISTANCE,SIDE,SAND_DEPTH,OBS_CODE,NOTES)
              
              depth
            }
              else{ if(survey == "kelp"){
                #################################
                ## Kelp 
                #################################
                
                temp <- data %>%
                  select(-ENTERED_BY,-CHECKED_BY) %>%
                  mutate(SP_CODE="MAPY",
                         DATE=format(as.Date(DATE,"%m/%d/%Y"),"%Y-%m-%d"),
                         SECTION=if_else(SECTION=="0-20","20",SECTION),
                         SECTION=if_else(SECTION=="21-40" | SECTION == "20-40" | SECTION == "40-20","40",SECTION),
                         SIDE=if_else(SIDE=="INSHORE","I",SIDE),
                         SIDE=if_else(SIDE=="OFFSHORE","O",SIDE),
                         NOTES=gsub(",", ".", NOTES,fixed = TRUE),
                         NOTES=if_else(is.na(NOTES),"",NOTES),
                         NOTES=if_else(!is.na(MULT_PLANT),paste0(MULT_PLANT,"plants sharing common holdfast;",NOTES),NOTES),
                         AREA=if_else(SIDE=="",40,20)) %>%
                  rename(QUAD=SECTION) %>%
                  select(-MULT_PLANT)
                
                
                #assign taxa info
                temp2 <-spp %>%
                  filter(SURVEY == "KELP") %>%
                  right_join(temp,by="SP_CODE") 
                
                colname<-c("YEAR","MONTH","DATE","SITE","TRANSECT","QUAD","SIDE","SP_CODE","FRONDS","HLD_DIAM","AREA","OBS_CODE","NOTES",
                           "SCIENTIFIC_NAME","COMMON_NAME","TAXON_KINGDOM","TAXON_PHYLUM","TAXON_CLASS","TAXON_ORDER",
                           "TAXON_FAMILY","TAXON_GENUS","GROUP","SURVEY","MOBILITY","GROWTH_MORPH")
                
                temp3 <- temp2[,colname]
                temp3
              }
                else{ if(survey != "kelp" | survey != "sand" | survey != "substrate" | survey != "upc" | survey != "quad" | survey != "crypticfish" | survey != "swath" | survey != "bigfish"){
  print("Survey not recognized. One of kelp, sand, substrate, upc, quad, swath, crypticfish, bigfish")}
              }}}}}}}}}





####################################################
## Function to check the raw data files
#############################################################

checks <- function(data, survey){
  if(survey == "crypticfish"){
    ##################################

    temp <- data %>%
      select(-ENTERED_BY,-CHECKED_BY) %>%
      mutate(SP_CODE=trimws(SP_CODE),
             DATE=format(as.Date(DATE,"%m/%d/%Y"),"%Y-%m-%d"),
             SECTION=if_else(SECTION=="0-20","20",SECTION),
             SECTION=if_else(SECTION=="21-40","40",SECTION),
             SIDE=if_else(SIDE=="IN","I",SIDE),
             SIDE=if_else(SIDE=="OFF","O",SIDE),
             NOTES=gsub(",", ".", NOTES,fixed = TRUE),
             NOTES=if_else(is.na(NOTES),"",NOTES),
             SIZE=if_else(COUNT==0&!is.na(COUNT),as.numeric(NA),as.numeric(SIZE)),
             AREA=if_else(SIDE=="",40,20)) %>% #CRYPTIC FISH SURVEY THE AREA CHANGED.
      rename(QUAD=SECTION) 
  }
  else{ if(survey == "bigfish"){
    ###################################
    temp <- data %>%
      select(-ENTERED_BY,-CHECKED_BY) %>%
      mutate(SP_CODE=trimws(SP_CODE),
             DATE=format(as.Date(DATE,"%m/%d/%Y"),"%Y-%m-%d"),
             SIZE=if_else(COUNT==0&!is.na(COUNT),as.numeric(NA),as.numeric(SIZE)),
             NOTES=gsub(",", ".", NOTES,fixed = TRUE),
             NOTES=if_else(is.na(NOTES),"",NOTES)
      )
  }
    else{ if(survey == "swath"){
      #####################################
      temp <- data %>%
        select(-ENTERED_BY,-CHECKED_BY) %>%
        mutate(DATE=format(as.Date(DATE,"%m/%d/%Y"),"%Y-%m-%d"),
               SIDE=trimws(SIDE),
               SP_CODE=trimws(SP_CODE),
               SECTION=SIDE,
               SECTION=if_else(SECTION=="IN 20"|SECTION=="OFF 20","20",SECTION),
               SECTION=if_else(SECTION=="IN 40"|SECTION=="OFF 40","40",SECTION),
               SECTION=as.numeric(SECTION),
               SIDE=if_else(SIDE=="IN 20"|SIDE=="IN 40","I",SIDE),
               SIDE=if_else(SIDE=="OFF 20"|SIDE=="OFF 40","O",SIDE),
               NOTES=gsub(",", ".", NOTES,fixed = TRUE),
               NOTES=if_else(is.na(NOTES),"",NOTES),
               SIZE=if_else(COUNT==0&!is.na(COUNT),as.numeric(NA),as.numeric(SIZE))
        ) %>%
        rename(QUAD=SECTION)
    }
      else{if(survey == "quad"){
        #####################################
        temp <- data %>%
          select(-ENTERED_BY,-CHECKED_BY) %>%
          filter(!is.na(SP_CODE)) %>% #REMOVE THE MISSING SHEET SITE. CHECK TO SEE IF THERE IS ANY MISSING SHEET HAPPEN ON THE TRANSECT. 
          mutate(DATE=format(as.Date(DATE,"%m/%d/%Y"),"%Y-%m-%d"),
                 SIDE=trimws(SIDE),
                 SIDE=if_else(SIDE=="ON"&!is.na(SIDE),"I",SIDE),
                 SIDE=if_else(SIDE=="OFF"&!is.na(SIDE),"O",SIDE),
                 SP_CODE=trimws(SP_CODE),
                 NOTES=gsub(",", ".", NOTES,fixed = TRUE),
                 NOTES=if_else(is.na(NOTES),"",NOTES)
          ) 
      }
        else{if(survey == "upc"){
         ###################################### 
          temp <- data %>%
            select(-ENTERED_BY,-CHECKED_BY) %>%
            mutate(SP_CODE=trimws(SP_CODE),
                   DATE=format(as.Date(DATE,"%m/%d/%Y"),"%Y-%m-%d"),
                   NOTES=gsub(",", ".", NOTES,fixed = TRUE),
                   NOTES=if_else(is.na(NOTES),"",NOTES),
                   SP_CODE=if_else(is.na(SP_CODE),"NOSP",SP_CODE)) #one species can't have percent_cover more than 100%
        }
          else{if(survey == "kelp"){
            ######################################
            temp <- data %>%
              select(-ENTERED_BY,-CHECKED_BY) %>%
              mutate(SP_CODE="MAPY",
                     DATE=format(as.Date(DATE,"%m/%d/%Y"),"%Y-%m-%d"),
                     SECTION=if_else(SECTION=="0-20","20",SECTION),
                     SECTION=if_else(SECTION=="21-40","40",SECTION),
                     SIDE=if_else(SIDE=="INSHORE","I",SIDE),
                     SIDE=if_else(SIDE=="OFFSHORE","O",SIDE),
                     NOTES=gsub(",", ".", NOTES,fixed = TRUE),
                     NOTES=if_else(is.na(NOTES),"",NOTES),
                     NOTES=if_else(!is.na(MULT_PLANT),paste0(MULT_PLANT,"plants sharing common holdfast;",NOTES),NOTES),
                     AREA=if_else(SIDE=="",40,20)) %>%
              rename(QUAD=SECTION) %>%
              select(-MULT_PLANT)
          }}}}}}
  #CHECK NF and see if it is NF
  c1 <- temp %>%
    filter(SP_CODE=="NF")
  
  # check for non matching spp
  c2 <- temp %>%
    left_join(spp,by="SP_CODE") %>%
    filter(!is.na(SP_CODE)&is.na(SURVEY))
  
  # check fish size
  if(survey == "kelp" | survey == "upc"){
    c3 <- "Not applicable."
    } 
  else{
    c3 <- temp %>%
      left_join(select(spp,SP_CODE,SIZE_MIN_REF,SIZE_MAX_REF),by="SP_CODE") %>%
      filter(!is.na(SP_CODE)&(SIZE>SIZE_MAX_REF|SIZE<SIZE_MIN_REF))
  }
  
  # # check repeated records .
  if(survey == "swath" | survey == "crypticfish" | survey == "quad"){
    c4 <- temp %>%
      group_by(DATE,SITE,TRANSECT,QUAD,SIDE,COUNT,SIZE) %>%
      summarise(freq=n()) %>%
      ungroup() %>%
      mutate(freq=if_else(freq>1,as.integer(2),freq)) %>%
      select(DATE,SITE,TRANSECT,freq) %>%
      distinct()%>%
      group_by(DATE,SITE,TRANSECT) %>%
      summarise(rep=min(freq)) %>%
      ungroup() %>%
      filter(rep==2) %>%
      left_join(temp,by=c("DATE","SITE","TRANSECT"))%>%
      arrange(DATE,SITE,TRANSECT,QUAD,SIDE,COUNT)
    
  } else{if(survey == "bigfish"){
    c4 <- temp %>%
      group_by(DATE,SITE,TRANSECT,SP_CODE,COUNT,SIZE) %>%
      summarise(freq=n()) %>%
      filter(freq>1)
  }
    else{if(survey == "upc"){
      c4 <- temp %>%
        group_by(DATE,SITE,TRANSECT, SIDE, DISTANCE, SP_CODE) %>%
        summarise(freq=n()) %>%
        filter(freq>1)
    }
      else{if(survey == "kelp"){
        c4 <- temp %>%
          group_by(DATE,SITE,TRANSECT, SIDE, QUAD, SP_CODE,FRONDS, HLD_DIAM) %>%
          summarise(freq=n()) %>%
          filter(freq>1)
      }}}}
  # # check if all survey were done in all transect, problems have been fixed directly in the archive sheet. 
  c5 <- temp %>%
    select(SITE,TRANSECT,YEAR) %>%
    distinct() %>%
    mutate(qs=paste0(SITE,TRANSECT),VALUE=1) %>%
    dcast(YEAR~qs,value.var = "VALUE",fill="missing") %>%
    arrange(YEAR)
  
  #check for if there is multiple dates on each transect, good
  c6 <- temp %>%
    select(SITE,TRANSECT,YEAR,MONTH,DATE)%>%
    distinct() %>%
    group_by(SITE,TRANSECT,YEAR,MONTH) %>%
    summarise(freq=n()) %>%
    filter(freq>1)
  
  #UPC specific check 
  
  #check for sand depth
  if(survey == "upc"){
  c7<-temp %>%
      filter(SUBSTRTE=="S") %>%
      filter(is.na(SAND_DEPTH))
  }else{if(survey != 'upc'){
    c7 <- print("Not applicable.")
  }}
  
  
  #check if 40 points on each side
  if(survey == "upc"){
    c8 <- temp %>%
       filter(POSITION==1) %>%
       group_by(YEAR,MONTH,DATE,SITE,TRANSECT,SIDE) %>%
       summarise(freq=n()) %>%
       filter(freq!=40)
  }else{if(survey != 'upc'){
    c8 <- print("Not applicable.")
  }}
  
  
  #check if missing substrate
  if(survey == "upc"){
    c9 <- temp %>%
      filter(POSITION==1) %>%
      filter(is.na(SUBSTRTE))
  }else{if(survey != 'upc'){
    c9 <- print("Not applicable.")
  }}
  
  
  #check if substrate column has any organisms 
  if(survey == "upc"){
    c10 <- temp %>%
       select(-SP_CODE) %>%
       rename(SP_CODE=SUBSTRTE) %>%
       left_join(spp,by="SP_CODE")%>%
       filter(POSITION==1)%>%
       filter(is.na(GROUP)|GROUP!="SUBSTRATE")
  } else{if(survey != 'upc'){
   c10 <- print("Not applicable.")
  }}
  
  # check if all survey were done in all quad-side
  if(survey == "upc"){
    c11 <- temp %>%
      select(SITE,TRANSECT,YEAR,MONTH,SIDE) %>%
      distinct() %>%
      mutate(qs=paste0(SIDE),VALUE=1) %>%
      dcast(SITE+TRANSECT+YEAR+MONTH~qs,value.var = "VALUE",fill="missing") %>%
      arrange(SITE,TRANSECT,YEAR)%>%
      filter(I=="missing"|O=="missing")
  }else{if(survey != 'upc'){
    c11 <- print("Not applicable.")
  }}
  
  
  
  out <- list(c1,c2,c3,c4,c5,c6, c7, c8, c9, c10, c11)
  
  names <- c("Are there species present?",
             "Are there any non-matching species?", 
             "Are species sizes realistic?", 
             "Are there repeated records?", 
             "How many surveys on a transect?", 
             "Are there multiple dates on a single transect?", 
             "If upc: check for sand depths?", 
             "If upc: Are there 40 points?", 
             "If upc: Is any substrate missing?", 
             "If upc: Are organisms and substrate mixed together?", 
             "If upc: Were all surveys done in all quad-sides?")
  
  names(out) <- names
  
  out
  
}


#######################################################
## Misc. functions
#######################################################

flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}



















