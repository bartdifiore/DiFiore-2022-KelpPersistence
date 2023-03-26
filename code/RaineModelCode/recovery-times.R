# README: code to compute the time it takes the benthic community to recover to within 10% of their pre-disturbance 
#values as a function of number of disturbances

library("tidyverse")

# parameters
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

# frond senescence
a <- 0.05 # amplitude of synchronized giant kelp frond senescence (1/d)
b <- 0.0125 # baseline rate of giant kelp frond senescence (1/d)
c <- 0.015#decay rate of synchrony in giant kelp frond senescence
d_syn <- -70 # displacement of the peak in synchrony is 70 days
period <- 110 # 110 day period between peaks in senescence
delay <- 150 # time to delay initiation of senescence is 150 days

#--------------------------------------
# get the undisturbed reference state (initial conditions)
#--------------------------------------

#model timepoints
tset.1 <- seq(from=1, to=4000, length.out=4000*82125/5475)

psi.A.mild<-.95#proportion of adult giant kelp fronds remaining after a mild storm
psi.A.sev<-0#proportion of adult giant kelp fronds remaining after a severe storm (change this value to change intensity of severe storms)
psi.b<-1#proportion of J, M, and I remaining after a severe storm (here no scraping)

#mild storms
t_mild_storm_set <- NaN # no mild storms
#t_mild_storm_set<-c(2:60, 335:425, 700:790, 1065:1155, 1430:1521, 1796:1886, 2161:2251, 2526:2616, 2891:2982, 3257:3347, 3622:3712)#time points at which mild storms occur (in this case, every day between Dec. 1 and Mar. 1 of each year)

t_severe_storm_set.1<-NaN # no severe storms

#tset.2 <- unique(sort(c(t_severe_storm_set.1, sim_tpts.1)))#new vector of time steps after any dates with storms that weren't already in sim_tpts.1 have been incorporated

t_all_storm_set.1<-NaN#all time points at which a storm occurs


#creating holding vectors for the outputs of the for loop; setting the first element in each holding vector
a_t <- NaN*tset.1;a_t[1]<-0#senescence rate of synchronized fronds
t_sen_set<-NaN*tset.1;t_sen_set[1]<-0#t_sen is a separate time vector; a_t is a function of t_sen
sA_set.1 <- NaN*tset.1; sA_set.1[1]<- 0#overall rate of frond senescence
f_syn <- NaN*tset.1 #proportion of the fronds that are synchronized 
G.simu.1 <- NaN*tset.1; G.simu.1[1]<-0.08 #giant kelp gametophytes
J.simu.1 <- NaN*tset.1; J.simu.1[1]<-0.00007 #juvenile giant kelp sporophytes
A.simu.1 <- NaN*tset.1; A.simu.1[1]<-0.5 #adult giant kelp fronds (start at 0.5)
M.simu.1 <- NaN*tset.1; M.simu.1[1]<-0.5 #understory macroalgae (start at 0.5)
I.simu.1 <- NaN*tset.1; I.simu.1[1] <- 0.5 #sessile invertebrates (start at 0.5)
psi.simu.1<-NaN*tset.1; psi.simu.1[1:2]<-1 #psi_A (proportion of giant kelp fronds remaining after a storm)


#simulate the model using a for loop
for(i in 2:length(tset.1)){
  dt <- tset.1[i]-tset.1[i-1]
  t_sen_set[i]<-t_sen_set[i-1]+dt
  G <- G.simu.1[i-1]
  J <- J.simu.1[i-1]
  A <- A.simu.1[i-1]
  M <- M.simu.1[i-1] 
  I <- I.simu.1[i-1]
  
  if(tset.1[i] %in% t_all_storm_set.1){
    t_sen_set[i]<-0
  }#any time a storm occurs, t_sen gets set to 0 (re-starting the synchronized senescence function a_t)
  
  if (t_sen_set[i] < delay){
    a_t[i] <- 0#the rate of synchronized senescence is 0 for the first 150 days after the most recent storm (delay=150)
  } else{
    a_t[i] <- a*(sin(2*pi/period*(t_sen_set[i]+d_syn))+1)#after delay days since the most recent storm, synchronized senescence is initiated
  }
  
  
  if (t_sen_set[i] < delay){
    f_syn[i] <- 1-psi.simu.1[i-1]#during the delay period following a storm, the proportion of synchronized fronds is equal to the prop. that got removed (1-psi_A) and are growing back
  } else{
    f_syn[i] <- (1-psi.simu.1[i-1])*exp(-c*(t_sen_set[i]-delay))#after the delay period, the proportion of synchronized fronds decays exponentially
  }
  
  
  sA <- a_t[i]*f_syn[i]+b*(1-f_syn[i]) #overall senescence rate is a composite of the senescence of synchronous fronds (f_syn) and asynchronous fronds (1-f_syn) 
  
  if(tset.1[i] %in% t_mild_storm_set){ #when a weak storm occurs
    psi.simu.1[i]<-psi.A.mild #psi_A=psi.A.mild 
    A.simu.1[i]<-psi.A.mild*A #A=psi.A.mild*A
    J.simu.1[i]<-1*J#J, M, and I are not affected by mild storms
    M.simu.1[i]<-1*M
    I.simu.1[i]<-1*I
    
  } else if (tset.1[i] %in% t_severe_storm_set.1){ #when a severe storm occurs
    psi.simu.1[i]<-psi.A.sev #psi_A=psi.A.sev
    A.simu.1[i]<-psi.A.sev*A #A=psi.A.sev*A
    J.simu.1[i]<-psi.b*J#the proportion of benthic organisms (J, M, and I) remaining after a storm is determined by psi.b
    M.simu.1[i]<-psi.b*M
    I.simu.1[i]<-psi.b*I
    
  } else {
    psi.simu.1[i]<-psi.simu.1[i-1]
    dJ <- (r_G*L_surface*exp(-k_l*A)*G-r_J*L_surface*exp(-k_l*A)*J-m_J*J^2)*dt #equation for change in juvenile giant kelp sporophyte population
    dA <- (r_J*L_surface*exp(-k_l*A)*J+g_A*A*L_surface*((K_A-A)/K_A)-sA*A)*dt#equation for change in adult giant kelp frond density
    dM <- (g_M*M*L_surface*exp(-k_l*A)*((S_T-M-alpha*I)/S_T)-s_M*M+sigma_M)*dt#equation for change in understory macroalgae cover
    dI <- (g_I*I*((S_T-I-beta*M)/S_T)-s_I*I+sigma_I)*dt#equation for change in sessile invert cover
    J.simu.1[i] <- J + dJ
    A.simu.1[i] <- A + dA
    M.simu.1[i] <- M + dM
    I.simu.1[i] <- I + dI}
  
  dG <- (sigma_ext+sigma_A*A-r_G*L_surface*exp(-k_l*A)*G-m_G*G^2)*dt#equation for change in giant kelp gametophyte population
  G.simu.1[i] <- G + dG
  sA_set.1[i] <- sA 
}


# plot results to check it equilibrated
#plot(x = tset.1, y = A.simu.1, col = "darkgreen", type = "l", lwd = 1.5, xlab = "time (years)", ylab = "Relative abundance", xaxt = "n", ylim = c(0,1))
#axis(side = 1, at = seq(from = 0, to = 4000, by = 365), labels = seq(from = 0, to = 10, by = 1))
#lines(x = tset.1, y = I.simu.1, col = "royalblue1", type = "l", lwd = 1.5)
#lines(x = tset.1, y = M.simu.1, col = "chocolate1", type = "l", lwd = 1.5)

#store the final values (the pre-disturbance state)
Gref <- G.simu.1[length(tset.1)]
Jref <- J.simu.1[length(tset.1)]
Aref <- A.simu.1[length(tset.1)]
Iref <- I.simu.1[length(tset.1)]
Mref <- M.simu.1[length(tset.1)]


# get the time steps for the simulations
#days in simulation: go to 50 years
recov_days <- length(seq(from = 0, to = 50*365, by = 1)) #18251 days
#model timepoints
recov_tpts.1 <- seq(from=1, to=recov_days, length.out=recov_days*82125/5475)

# get the disturbance sets
# number of disturbances to iterate over (from 1 to 25)
ndist <- c(1:25)
# full disturbance set (days of each disturbance: one per year)
t_severe_storm_set_full <- seq(from = 365, to = 25*365, by = 365)
#mild storms
t_mild_storm_set <- NaN # no mild storms


#--------------------------------------
# get recovery times for each number of disturbances
#--------------------------------------


# make a function that takes the scraping parameter as an input and returns the time it takes the benthic community to 
#recover (i.e., both algae and inverts are within 10% of their predisturbed values) for each number of disturbances

recov_fun <- function(psi.b.val){# psi.b = fraction of pre-disturbance algae and inverts remaining after disturbance
  
  # holding vector of recovery times
  recov_times.f <- NaN*ndist
  
  # holding list for the timeseries (to check what they look like)
  #tsims <- vector(mode='list', length=length(ndist))
  #Isims <- vector(mode='list', length=length(ndist))
  #Msims <- vector(mode='list', length=length(ndist))
  #Asims <- vector(mode='list', length=length(ndist))
  
  for(j in 1:length(ndist)){
    
    # number of disturbances
    t_severe_storm_set.j <- t_severe_storm_set_full[1:j]
    
    tset.j <- unique(sort(c(t_severe_storm_set.j, recov_tpts.1)))#new vector of time steps after any dates with storms that weren't already in sim_tpts.1 have been incorporated
    
    t_all_storm_set.j<-t_severe_storm_set.j#all time points at which a storm occurs
    
    # holding vectors
    #creating holding vectors for the outputs of the for loop; setting the first element in each holding vector
    a_t <- NaN*tset.j;a_t[1]<-0#senescence rate of synchronized fronds
    t_sen_set<-NaN*tset.j;t_sen_set[1]<-0#t_sen is a separate time vector; a_t is a function of t_sen
    sA_set.1 <- NaN*tset.j; sA_set.1[1]<- 0#overall rate of frond senescence
    f_syn <- NaN*tset.j #proportion of the fronds that are synchronized 
    G.simu.j <- NaN*tset.j; G.simu.j[1]<-Gref #giant kelp gametophytes
    J.simu.j <- NaN*tset.j; J.simu.j[1]<-Jref #juvenile giant kelp sporophytes
    A.simu.j <- NaN*tset.j; A.simu.j[1]<-Aref #adult giant kelp fronds (choose A0, I0, M0 so model starts near dynamic equilibrium)
    M.simu.j <- NaN*tset.j; M.simu.j[1]<-Mref #understory macroalgae 
    I.simu.j <- NaN*tset.j; I.simu.j[1] <- Iref #sessile invertebrates
    psi.simu.j<-NaN*tset.j; psi.simu.j[1:2]<-1 #psi_A (proportion of giant kelp fronds remaining after a storm)
    
    
    #simulate the model using a for loop
    for(i in 2:length(tset.j)){
      dt <- tset.j[i]-tset.j[i-1]
      t_sen_set[i]<-t_sen_set[i-1]+dt
      G <- G.simu.j[i-1]
      J <- J.simu.j[i-1]
      A <- A.simu.j[i-1]
      M <- M.simu.j[i-1] 
      I <- I.simu.j[i-1]
      
      if(tset.j[i] %in% t_all_storm_set.j){
        t_sen_set[i]<-0
      }#any time a storm occurs, t_sen gets set to 0 (re-starting the synchronized senescence function a_t)
      
      if (t_sen_set[i] < delay){
        a_t[i] <- 0#the rate of synchronized senescence is 0 for the first 150 days after the most recent storm (delay=150)
      } else{
        a_t[i] <- a*(sin(2*pi/period*(t_sen_set[i]+d_syn))+1)#after delay days since the most recent storm, synchronized senescence is initiated
      }
      
      
      if (t_sen_set[i] < delay){
        f_syn[i] <- 1-psi.simu.j[i-1]#during the delay period following a storm, the proportion of synchronized fronds is equal to the prop. that got removed (1-psi_A) and are growing back
      } else{
        f_syn[i] <- (1-psi.simu.j[i-1])*exp(-c*(t_sen_set[i]-delay))#after the delay period, the proportion of synchronized fronds decays exponentially
      }
      
      
      sA <- a_t[i]*f_syn[i]+b*(1-f_syn[i]) #overall senescence rate is a composite of the senescence of synchronous fronds (f_syn) and asynchronous fronds (1-f_syn) 
      
      if(tset.j[i] %in% t_mild_storm_set){ #when a weak storm occurs
        psi.simu.j[i]<-psi.A.mild #psi_A=psi.A.mild 
        A.simu.j[i]<-psi.A.mild*A #A=psi.A.mild*A
        J.simu.j[i]<-1*J#J, M, and I are not affected by mild storms
        M.simu.j[i]<-1*M
        I.simu.j[i]<-1*I
        
      } else if (tset.j[i] %in% t_severe_storm_set.j){ #when a severe storm occurs
        psi.simu.j[i]<-psi.A.sev #psi_A=psi.A.sev
        A.simu.j[i]<-psi.A.sev*A #A=psi.A.sev*A
        J.simu.j[i]<-psi.b.val*J#the proportion of benthic organisms (J, M, and I) remaining after a storm is determined by psi.b.val
        M.simu.j[i]<-psi.b.val*M
        I.simu.j[i]<-psi.b.val*I
        
      } else {
        psi.simu.j[i]<-psi.simu.j[i-1]
        dJ <- (r_G*L_surface*exp(-k_l*A)*G-r_J*L_surface*exp(-k_l*A)*J-m_J*J^2)*dt #equation for change in juvenile giant kelp sporophyte population
        dA <- (r_J*L_surface*exp(-k_l*A)*J+g_A*A*L_surface*((K_A-A)/K_A)-sA*A)*dt#equation for change in adult giant kelp frond density
        dM <- (g_M*M*L_surface*exp(-k_l*A)*((S_T-M-alpha*I)/S_T)-s_M*M+sigma_M)*dt#equation for change in understory macroalgae cover
        dI <- (g_I*I*((S_T-I-beta*M)/S_T)-s_I*I+sigma_I)*dt#equation for change in sessile invert cover
        J.simu.j[i] <- J + dJ
        A.simu.j[i] <- A + dA
        M.simu.j[i] <- M + dM
        I.simu.j[i] <- I + dI}
      
      dG <- (sigma_ext+sigma_A*A-r_G*L_surface*exp(-k_l*A)*G-m_G*G^2)*dt#equation for change in giant kelp gametophyte population
      G.simu.j[i] <- G + dG
      sA_set.1[i] <- sA 
      
      # test values
      testI <- abs((I.simu.j[i]-Iref)/Iref)
      testM <- abs((M.simu.j[i]-Mref)/Mref)
      testtime <- tset.j[i]
      
      if(testI<0.1 & testM<0.1 & testtime > t_severe_storm_set_full[j] + 80){ # if the system has recovered (the benthic community is within 10% of starting conditions, and it has been at least 80d since the last disturbance to avoid initial post-disturbance abundances near starting conditions), break the loop
        
        break
        
      }
    }
    
    # record the ith timepoint (the point when the benthic community recovered)
    recov_times.f[j] <- tset.j[i] - 365*ndist[j] # days since last disturbance
    
    # save the timeseries up to point of recovery
    #tsims[[j]] <- tset.j[which(tset.j<=tset.j[i])]
    #Isims[[j]] <- I.simu.j[which(tset.j<=tset.j[i])]
    #Msims[[j]] <- M.simu.j[which(tset.j<=tset.j[i])]
    #Asims[[j]] <- A.simu.j[which(tset.j<=tset.j[i])]
    
    
  }
  
  return(recov_times.f)
  
}


# run this function for different scraping scenarios (takes around a min)
# no scraping
rtimes1 <- recov_fun(1)

# scraping = 0.25 (75% remain)
rtimes75 <- recov_fun(0.75)

# scraping = 0.5 (50% remain)
rtimes50 <- recov_fun(0.5)

# scraping = 0.75 (25% remain)
rtimes25 <- recov_fun(0.25)


# turn into data frame (also add the 0,0 cases -> when there is no disturbance, recovery time = 0d))
rtimesdf <- tibble(
  c(0,ndist),
  c(0,rtimes1),
  c(0,rtimes75),
  c(0,rtimes50),
  c(0,rtimes25)
)

# write this to a csv file
write.csv(rtimesdf, "data/intermediary/ModelRecovTimes.csv")

# plot results (note here added point at 0,0 -> when there is no disturbance, recovery time = 0d)
plot(x= rtimesdf$ndist, y = rtimesdf$rtimes1, type = "l", xlab = "number of disturbances", ylab = "recovery time (d)", ylim = c(0, 1100))
lines(x= rtimesdf$ndist, y = rtimesdf$rtimes75, type = "l", lty = 2)
lines(x= rtimesdf$ndist, y = rtimesdf$rtimes50, type = "l", lty = 3)
lines(x= rtimesdf$ndist, y = rtimesdf$rtimes25, type = "l", lty = 4)
legend("bottomright", legend = c("0", "0.25","0.5", "0.75"), title = "Fraction benthos \ndisturbed", lty = c(1, 2, 3, 4), bty = "n")






