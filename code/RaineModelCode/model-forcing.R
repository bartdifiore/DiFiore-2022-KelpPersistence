# README: code for using the theoretical model to run timeseries at each transect with kelp abundance forced 
# to be equal to the observed quarterly means (with biomass scaled to be between 0-1)

# plotting colors
ua_col <-  "chocolate1"
si_col <-  "royalblue1"
kelp_col <-  "darkgreen"

# packages
library("tidyverse")
library("zoo")#for working with dates in quarter format
library("lubridate")

# ------------------
# data formatting
# ------------------

# Import the canopy biomass data
cnpydt <- read.csv("data/intermediary/kelpcanopytimeseries_wide.csv")

# Rescale the kelp biomass to the max value observed in each transect
scaledcnpy <- NaN*cnpydt
scaledcnpy[ ,1] <- cnpydt[, 1]#first column (date) stays the same

for(i in 2:ncol(cnpydt)){
  scaledcnpy[, i] <- cnpydt[, i]/max(cnpydt[, i], na.rm = TRUE)#AQUE3 had one Na value
}

#Format the dates
#Start date = Jan 1, 2008, end date = Oct. 1, 2018.

#convert quarter to days: R stores dates numerically as days since 1/1/1970, so to convert day to a number 
#representing the sampling day, convert the year.quarter values to yearqtr format that can be converted to 
#date format, convert this to numeric (days since 1/1/1970), then subtract the number of days between 
#1/1/1970 and 1/1/2008 to make 1/1/2008 correspond to day 1
daysbtw <- length(seq(from = as.Date("1970-01-01"), to = as.Date("2008-01-01"), by = 'day'))-2
day.1 <- as.numeric(as.Date(as.yearqtr(scaledcnpy$year.quarter)))-daysbtw

#Get the timepoints in the model
#get the list of time points in days, with day 1 being the first day of observation/model simulation
tot_days <- length(seq(from = as.Date("2008-01-01"), to = as.Date("2018-10-01"), by = 'day')) #3927 days

#get the model timepoints: want to keep the interval in the for loop similar to what it was in the simulations 
#for Detmer et al. 2021 paper, which were ~ 5475 days long with a length.out of 82125, meaning each step in the 
#for loop is ~1/15 of a day
timepoints.1 <- seq(from=1, to=tot_days, length.out=tot_days*82125/5475)

#the exact days of the starts of the quarters aren't necessarily in timepoints.1 
#(e.g., timepoints.1 might contain day 92.111 instead of 92), so use unique() to make sure the exact days 
#are in the model timepoints
timepoints.2 <- unique(sort(c(timepoints.1, day.1)))

#Want to make column for days since start (day 1=1), date, month, year. 
#First make the dates, then the days since start, then the month, year

#add model timepoint, date, month, and year columns to the dataset (will be filled in later)
tmpt <- NaN*day.1
date <- NaN*day.1
month <- NaN*day.1
year <- NaN*day.1

#add these columns and the time in days (day.1) to the dataset
scaledcnpy2 <- cbind(tmpt, day.1,date, month, year, scaledcnpy)

#now want to expand the data so the number of rows is equal to the number of timepoints in 
#the model simulations, where the observed values at all dates in a quarter correspond to the observed 
#mean biomass for that quarter. 
cnpydt_full <- data.frame(matrix(ncol=ncol(scaledcnpy2), nrow=length(timepoints.2)))
colnames(cnpydt_full) <- colnames(scaledcnpy2)

cnpydt_full[, 1] <- timepoints.2#make the model timepoints the first column

match_rows <- which(cnpydt_full$tmpt %in% scaledcnpy2$day.1)#these are the rows for which the day corresponds 
#to the year.quarter day

#fill in these rows with the corresponding row from the scaledcnpy dataframe
for(i in 1:length(scaledcnpy2$day.1)){
  
  row.i <- match_rows[i]
  
  cnpydt_full[row.i, 2:ncol(cnpydt_full)] <- scaledcnpy2[i, 2:ncol(scaledcnpy2)]#don't need to fill in 
  #the first column because this is already timepoints.2
  
}

#fill in the date, month, and year columns. floor(timepoints.2) rounds each value of timepoints.2 down 
#to the nearest integer (day starting on 1/1/2008)
cnpydt_full1 <- cnpydt_full %>% 
  mutate(date=as.Date(floor(timepoints.2)+daysbtw), 
         month=month(as.Date(floor(timepoints.2)+daysbtw)), 
         year=year(as.Date(floor(timepoints.2)+daysbtw)), 
         quarter = if_else(as.numeric(month) %in% c(1, 2, 3), 1, 
                           if_else(as.numeric(month) %in% c(4, 5, 6), 2, 
                                   if_else(as.numeric(month) %in% c(7, 8, 9), 3, 4)))) %>% 
  mutate(year.quarter = paste(year, quarter, sep = "."))


#fill in all the Nas with the quarterly kelp averages
cnpydt_full2 <- cnpydt_full1 %>% fill(names(cnpydt_full)) %>% select(-quarter)#fill in all the Nas


# ------------------
# model initialization
# ------------------

# make the holding data frames for model simulations

#add empty "group" column to the full data frame
group <- NaN*timepoints.2
cnpydt_full3 <- cbind(group, cnpydt_full2)

kelp.runs <- data.frame(matrix(ncol=ncol(cnpydt_full3), nrow=length(timepoints.2)))
colnames(kelp.runs) <- colnames(cnpydt_full3)
kelp.runs[ ,2:7] <- cnpydt_full3[ ,2:7]
kelp.runs[ ,1] <- rep("kelp", length(timepoints.2))
kelp.runs[1, 8:ncol(cnpydt_full3)] <- cnpydt_full3[1, 8:ncol(cnpydt_full3)]#set initial kelp abundance to the observed mean for the first quarter

invert.runs <- data.frame(matrix(ncol=ncol(cnpydt_full3), nrow=length(timepoints.2)))
colnames(invert.runs) <- colnames(cnpydt_full3)
invert.runs[ ,2:7] <- cnpydt_full3[ ,2:7]
invert.runs[ ,1] <- rep("invert", length(timepoints.2))

algae.runs <- data.frame(matrix(ncol=ncol(cnpydt_full3), nrow=length(timepoints.2)))
colnames(algae.runs) <- colnames(cnpydt_full3)
algae.runs[ ,2:7] <- cnpydt_full3[ ,2:7]
algae.runs[ ,1] <- rep("algae", length(timepoints.2))

#number of sites (individual transects) to iterate over
n_sites <- length(8:ncol(cnpydt_full3))

# get initial invert and understory macroalgal percent cover at each site by running the model to equilibrium 
#with kelp help constant at the initial observed quarterly mean

#model parameters 
L_surface<-1000 #surface irradiance (mol/m^2/s)
K_A<-1 #giant kelp frond density carrying capacity (prop. max. fronds/m^2)
k_l<- -log(0.1)/K_A #giant kelp frond extinction coeff. (m^2/fronds)
sigma_ext<-0.0001#supply rate of giant kelp gametophytes (via dispersal of spores) from external populations (ind./m^2/d)
sigma_A<-0.01#rate of gametophyte production (via spore production) by local giant kelp sporophytes (ind./fronds/d)
sigma_M<-0.001#rate of increase in fractional cover of understory macroalgae from external populations (1/d)
sigma_I<-0.0005#rate of increase in fractional cover of sessile inverts from external populations (1/d)
r_G<-0.00000001#recruitment rate of giant kelp gametophytes to juvenile sporophytes (m^2*s/d/mol)
r_J<-0.00001#maturation rate of juvenile sporophytes to adult sporophyte fronds (m^2*s/mol/d); note it is assumed that an individual juvenile sporophyte has 1 frond
g_A<-0.00009#growth rate of adult fronds (m^4*s/mol/d/fronds)
g_M<-0.00006#growth rate of understory macroalgae (m^2*s/mol/d)
g_I<-0.008#growth rate of sessile invertebrates (1/d)
m_G<-1#mortality rate of giant kelp gametophytes (m^2/ind./d)
m_J<-1#mortality rate of juvenile giant kelp sporophytes (m^2/ind./d)
s_M<-0.0099#senescence rate of understory macroalgae (1/d); DEFAULT= 0.009, CHANGED to 0.0099 
s_I<-0.0022#senescence rate of sessile invertebrates (1/d); DEFAULT= 0.002, CHANGED to 0.0022 
S_T<-1#total substrate space available for understory macroalgal and invertebrate growth
alpha<-0.8#competition coefficient of sessile inverts on understory macroalgae
beta<-1.25#competition coefficient of understory macroalgae on sessile inverts


tset.1.1 <- timepoints.2#vector of time steps to iterate over

# run the model for each site
for(j in 1:n_sites){
  
  # initial conditions
  A0<-kelp.runs[1, j+7]#initial kelp abundance at each site
  M0<-0.5
  I0 <- 0.5
  
  #creating holding vectors for the outputs of the for loop; setting the first element in each holding vector
  M.simu.1 <- NaN*tset.1.1; M.simu.1[1]<-M0 #understory macroalgae
  I.simu.1 <- NaN*tset.1.1; I.simu.1[1] <- I0 #sessile invertebrates
  
  #simulate the model using a for loop
  for(i in 2:length(tset.1.1)){
    dt <- tset.1.1[i]-tset.1.1[i-1]#time step
    A <- A0
    M <- M.simu.1[i-1] 
    I <- I.simu.1[i-1]
    
    #calculate the change in benthic macroalgae and sessile inverts
    dM <- (g_M*M*L_surface*exp(-k_l*A)*((S_T-M-alpha*I)/S_T)-s_M*M+sigma_M)*dt#equation for change in understory macroalgae cover
    dI <- (g_I*I*((S_T-I-beta*M)/S_T)-s_I*I+sigma_I)*dt#equation for change in sessile invert cover
    
    M.simu.1[i] <- M + dM
    I.simu.1[i] <- I + dI
    
    
  }
  
  invert.runs[1, j+7] <- I.simu.1[length(tset.1.1)]#I0 for the jth site
  algae.runs[1, j+7] <- M.simu.1[length(tset.1.1)]#M0 for the jth site
  
}


# ------------------
# model forcing simulations
# ------------------

# holding data frames
invert.runs2 <- data.frame(matrix(ncol=ncol(cnpydt_full3), nrow=length(timepoints.2)))
colnames(invert.runs2) <- colnames(cnpydt_full3)
invert.runs2[ ,2:7] <- cnpydt_full3[ ,2:7]
invert.runs2[ ,1] <- rep("invert", length(timepoints.2))
invert.runs2[1, 8:ncol(cnpydt_full3)] <- invert.runs[1, 8:ncol(cnpydt_full3)]#set initial conditions

algae.runs2 <- data.frame(matrix(ncol=ncol(cnpydt_full3), nrow=length(timepoints.2)))
colnames(algae.runs2) <- colnames(cnpydt_full3)
algae.runs2[ ,2:7] <- cnpydt_full3[ ,2:7]
algae.runs2[ ,1] <- rep("algae", length(timepoints.2))
algae.runs2[1, 8:ncol(cnpydt_full3)] <- algae.runs[1, 8:ncol(cnpydt_full3)]#set initial conditions

# simulate the model (takes ~30 seconds or so)

tset.1.1 <- timepoints.2#vector of time steps to iterate over

for(j in 1:n_sites){#for each transect at each site
  
  #set initial conditions
  A0<-cnpydt_full3[1, j+7]#initial kelp abundance
  
  M0<-algae.runs2[1, j+7]#initial benthic macroalgal abundance
  I0<-invert.runs2[1, j+7]#initial invert abundance
  
  #creating holding vectors for the outputs of the for loop; setting the first element in each holding vector
  
  A.simu.1 <- NaN*tset.1.1; A.simu.1[1]<-A0 #adult giant kelp fronds
  M.simu.1 <- NaN*tset.1.1; M.simu.1[1]<-M0 #understory macroalgae
  I.simu.1 <- NaN*tset.1.1; I.simu.1[1] <- I0 #sessile invertebrates
  
  
  #simulate the model using a for loop
  for(i in 2:length(tset.1.1)){
    dt <- tset.1.1[i]-tset.1.1[i-1]#time step
    A <- A.simu.1[i-1]
    M <- M.simu.1[i-1] 
    I <- I.simu.1[i-1]
    
    A.simu.1[i] <-cnpydt_full3[i, j+7]#set A at time t equal to observed mean
    
    #calculate the change in all the other parts of the model (benthic macroalgae, sessile inverts)
    dM <- (g_M*M*L_surface*exp(-k_l*A)*((S_T-M-alpha*I)/S_T)-s_M*M+sigma_M)*dt#equation for change in understory macroalgae cover
    dI <- (g_I*I*((S_T-I-beta*M)/S_T)-s_I*I+sigma_I)*dt#equation for change in sessile invert cover
    M.simu.1[i] <- M + dM
    I.simu.1[i] <- I + dI
    
  }
  
  invert.runs2[ , j+7] <- I.simu.1#invert time series for the jth site/transect
  algae.runs2[ , j+7] <- M.simu.1#algae for the jth site/transect
  
}

# get the quarterly mean values (to match the kelp data)
invert.runs.qmeans <- invert.runs2 %>% group_by(group, year.quarter) %>% summarize_at(colnames(invert.runs2[8:75]), list(~mean(.)))

algae.runs.qmeans <- algae.runs2 %>% group_by(group, year.quarter) %>% summarize_at(colnames(algae.runs2[8:75]), list(~mean(.)))

kelp.runs.qmeans <- cnpydt_full3 %>% mutate(group = "kelp")%>% group_by(group, year.quarter) %>% summarize_at(colnames(cnpydt_full3[8:75]), list(~mean(.)))

# check the kelp values are the same as the original data
#plot(x = scaledcnpy$year.quarter, y = scaledcnpy$AQUE2, type = "l", col = "black", lwd=2)
#lines(x = scaledcnpy$year.quarter, y = kelp.runs.qmeans$AQUE2, type = "o", col = "palegreen3")

# bind all the model runs and kelp together
FullModelRunsAllSites <- rbind(kelp.runs.qmeans, invert.runs.qmeans, algae.runs.qmeans) 

# ------------------
# export the data
# ------------------

# export the data
write.csv(FullModelRunsAllSites, "data/Large_files/Raine_modeloutput/FullModelRunsAllSites.csv")

# ------------------
# check plots
# ------------------

# code copied from "5b_percover-analysis.R"

dist <- read.csv("data/intermediary/persistence_metrics.csv")

rd_ts  <- FullModelRunsAllSites %>%
  pivot_longer(cols = ABUR1:WOOD3, names_to = "id", values_to = "percent_cover") %>%
  group_by(id) %>%
  separate(id, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])", remove = F) %>%
  mutate(transect = as.numeric(transect)) %>%
  left_join(dist) %>%
  separate(year.quarter, into = c("year", "quarter"), sep = "[.]") %>%
  mutate(quarter = case_when(quarter == 1 ~ ".0", 
                             quarter == 2 ~ ".25", 
                             quarter == 3 ~ ".5", 
                             quarter == 4 ~ ".75")) %>%
  mutate(year.quarter = as.numeric(paste(year, quarter, sep = "")))


rd_predictions <- rd_ts %>% 
  filter(year.quarter == 2018.5, group != "kelp")

p1 <- ggplot(rd_predictions, aes(x = timesincelastabsent, y = percent_cover))+
  geom_point(aes(color = group))+
  geom_smooth(method = "lm", aes(color = group), se = F)+
  scale_color_manual(values = c(ua_col, si_col))+
  labs(x = "", y = "Theoretical prediction\nrelative abundance", color = "")+
  #theme_bd()+
  theme(legend.position = "top")

p4 <- ggplot(rd_predictions, aes(x = prop_zero, y = percent_cover))+
  geom_point(aes(color = group))+
  geom_smooth(method = "lm", aes(color = group), se = F)+
  scale_color_manual(values = c(ua_col, si_col))+
  labs(x = "", y = "Theoretical prediction\nrelative abundance", color = "")+
  #theme_bd()+
  theme(legend.position = "top")

p1

p4








