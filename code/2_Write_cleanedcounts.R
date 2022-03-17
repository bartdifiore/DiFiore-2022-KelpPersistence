
source("code/1_setup.R")

##########################################
## Big fish
##########################################

raw.bf <- read.csv("data/data_archive/BartsBenthic_Fish_2018.csv", header = T, na=c(".","","-99999"), stringsAsFactors = F)

checks(raw.bf, survey = "bigfish") 
  #manually fixed errors in data frame to build "BartsBenthic_Fish_2018_cleaned.csv"


dat <- read.csv("data/cleaned/BartsBenthic_Fish_2018_cleaned.csv", header=TRUE, na=c(".","","-99999"), stringsAsFactors = F)

checks(dat, survey = "bigfish")
df.bf <- sbc_clean(dat, survey = "bigfish")

##########################################
## Cryptic fish
##########################################

raw.cf <- read.csv("data/data_archive/BartsBenthic_Cryptic_Fish_2018.csv", header = T, na=c(".","","-99999"), stringsAsFactors = F)

checks(raw.cf, survey = "crypticfish") 

dat <- read.csv("data/cleaned/BartsBenthic_Cryptic_Fish_2018_cleaned.csv", header = T, na=c(".","","-99999"), stringsAsFactors = F)

checks(dat, survey = "crypticfish")
df.cf <- sbc_clean(dat, survey = "crypticfish")

##########################################
## Swath
##########################################

raw.sw <- read.csv("data/data_archive/BartsBenthic_Swath_2018.csv", header = T, na=c(".","","-99999"), stringsAsFactors = F)

checks(raw.sw, survey = "swath")


dat <- read.csv("data/cleaned/BartsBenthic_Swath_2018_cleaned.csv", header = T, na=c(".","","-99999"), stringsAsFactors = F) #this cleaned file needs to be recleaned based on a review of the TANK1 data sheet.

checks(dat, survey = "swath")
df.sw <- sbc_clean(dat, survey = "swath")

##########################################
## Quad
##########################################

raw.qd <- read.csv("data/data_archive/BartsBenthic_Quad_2018.csv", header = T, na=c(".","","-99999"), stringsAsFactors = F)

checks(raw.qd, survey = "quad")

dat <- read.csv("data/cleaned/BartsBenthic_Quad_2018_cleaned.csv", header = T, na=c(".","","-99999"), stringsAsFactors = F)

checks(dat, survey = "quad")
df.qd <- sbc_clean(dat, survey = "quad")

##########################################
## UPC
##########################################

raw.upc <- read.csv("data/data_archive/BartsBenthic_UPC_2018.csv", header = T, na=c(".","","-99999"), stringsAsFactors = F)

checks(raw.upc, survey = "upc")

dat <- read.csv("data/cleaned/BartsBenthic_UPC_2018_cleaned.csv", header = T, na=c(".","","-99999"), stringsAsFactors = F)

checks(dat, survey = "upc")
df.upc <- sbc_clean(dat, survey = "upc")
df.sub <- sbc_clean(dat, survey = "substrate")
write.csv(df.sub, "data/intermediary/substrate/BBenthic_substrate.csv", row.names = F)
df.sand <- sbc_clean(dat, survey = "sand")


##########################################
## Kelp
##########################################

raw.kelp <- read.csv("data/data_archive/BartsBenthic_Kelp_2018.csv", header = T, na=c(".","","-99999"), stringsAsFactors = F)

checks(raw.kelp, survey = "kelp")

dat <- read.csv("data/cleaned/BartsBenthic_Kelp_2018_cleaned.csv", header = T, na=c(".","","-99999"), stringsAsFactors = F)

checks(dat, survey = "kelp")
df.kp <- sbc_clean(dat, survey = "kelp")

##########################################################
# merge the quad and swath together
##########################################################


comb_qdsw <- bind_rows(df.qd,df.sw) 


colname<-c("YEAR","MONTH","DATE","SITE","TRANSECT","QUAD","SIDE","SP_CODE","SIZE","COUNT","AREA",
           "SCIENTIFIC_NAME","COMMON_NAME","TAXON_KINGDOM","TAXON_PHYLUM","TAXON_CLASS","TAXON_ORDER",
           "TAXON_FAMILY","TAXON_GENUS","GROUP","SURVEY","MOBILITY","GROWTH_MORPH")


comb_qdsw <- comb_qdsw[,colname]

##########################################################
# merge the annual big fish and the cryptic fish together#
##########################################################

cf <- df.cf %>%
  mutate(VIS=NA)

bf <- df.bf %>%
  mutate(QUAD=NA,SIDE=NA)

data <- rbind(bf,cf)

colname<-c("YEAR","MONTH","DATE","SITE","TRANSECT","QUAD","SIDE","VIS","SP_CODE","SIZE","COUNT","AREA",
           "SCIENTIFIC_NAME","COMMON_NAME","TAXON_KINGDOM","TAXON_PHYLUM","TAXON_CLASS","TAXON_ORDER",
           "TAXON_FAMILY","TAXON_GENUS","GROUP","SURVEY","MOBILITY","GROWTH_MORPH")

comb_bfcf <- data[,colname]




