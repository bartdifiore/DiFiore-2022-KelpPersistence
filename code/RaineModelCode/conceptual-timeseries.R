# README: code to make the time series of theoretical high and low kelp variability sites

# plotting colors
ua_col <-  "chocolate1"
si_col <-  "royalblue1"
kelp_col <-  "darkgreen"

# ------------------
# define parameters and get initial conditions
# ------------------

# model parameters (same as for model forcing simulations- no scraping)
# same as for model forcing simulations
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

# timesteps: simulate for 10 years
#days in simulation
sim_days <- length(seq(from = as.Date("2008-01-01"), to = as.Date("2018-10-01"), by = 'day')) #3927 days

#model timepoints
sim_tpts.1 <- seq(from=1, to=sim_days, length.out=sim_days*82125/5475)

# initial conditions: will start kelp at a lowish value (25% of carrying capacity) in the timeseries, so first get the invert and algae abundances 
#that correspond to this by running the model to equilibrium with kelp held at this value
Astart <- 0.25

tset.1 <- sim_tpts.1#vector of time steps to iterate over

A0<-Astart#initial kelp abundance 

M0<-0.5
I0 <- 0.5

#creating holding vectors for the outputs of the for loop; setting the first element in each holding vector
M.simu.1 <- NaN*tset.1; M.simu.1[1]<-M0 #understory macroalgae
I.simu.1 <- NaN*tset.1; I.simu.1[1] <- I0 #sessile invertebrates

#simulate the model using a for loop
for(i in 2:length(tset.1)){
  dt <- tset.1[i]-tset.1[i-1]#time step
  A <- A0
  M <- M.simu.1[i-1] 
  I <- I.simu.1[i-1]
  
  #calculate the change in benthic macroalgae and sessile inverts
  dM <- (g_M*M*L_surface*exp(-k_l*A)*((S_T-M-alpha*I)/S_T)-s_M*M+sigma_M)*dt#equation for change in understory macroalgae cover
  dI <- (g_I*I*((S_T-I-beta*M)/S_T)-s_I*I+sigma_I)*dt#equation for change in sessile invert cover
  
  M.simu.1[i] <- M + dM
  I.simu.1[i] <- I + dI
  
  
}

Istart <- I.simu.1[length(tset.1)]#starting value for inverts
Mstart <- M.simu.1[length(tset.1)]#starting value for understory macroalgae


# ------------------
# low variability site
# ------------------

psi.A.mild<-.95#proportion of adult giant kelp fronds remaining after a mild storm
psi.A.sev<-0#proportion of adult giant kelp fronds remaining after a severe storm (change this value to change intensity of severe storms)
psi.b<-1#proportion of J, M, and I remaining after a severe storm (here no scouring)

#mild storms
t_mild_storm_set <- NaN # no mild storms

t_severe_storm_set.1<-c(730) # time points at which severe storms occur

tset.2 <- unique(sort(c(t_severe_storm_set.1, sim_tpts.1)))#new vector of time steps after any dates with storms that weren't already in sim_tpts.1 have been incorporated

t_all_storm_set.1<-t_severe_storm_set.1#all time points at which a storm occurs


#creating holding vectors for the outputs of the for loop; setting the first element in each holding vector
a_t <- NaN*tset.2;a_t[1]<-0#senescence rate of synchronized fronds
t_sen_set<-NaN*tset.2;t_sen_set[1]<-0#t_sen is a separate time vector; a_t is a function of t_sen
sA_set.1 <- NaN*tset.2; sA_set.1[1]<- 0#overall rate of frond senescence
f_syn <- NaN*tset.2 #proportion of the fronds that are synchronized 
G.simu.2 <- NaN*tset.2; G.simu.2[1]<-G0 #giant kelp gametophytes
J.simu.2 <- NaN*tset.2; J.simu.2[1]<-J0 #juvenile giant kelp sporophytes
A.simu.2 <- NaN*tset.2; A.simu.2[1]<-Astart #adult giant kelp fronds (choose A0, I0, M0 so model starts near dynamic equilibrium)
M.simu.2 <- NaN*tset.2; M.simu.2[1]<-Mstart #understory macroalgae 
I.simu.2 <- NaN*tset.2; I.simu.2[1] <- Istart #sessile invertebrates
psi.simu.2<-NaN*tset.2; psi.simu.2[1:2]<-1 #psi_A (proportion of giant kelp fronds remaining after a storm)


#simulate the model using a for loop
for(i in 2:length(tset.2)){
  dt <- tset.2[i]-tset.2[i-1]
  t_sen_set[i]<-t_sen_set[i-1]+dt
  G <- G.simu.2[i-1]
  J <- J.simu.2[i-1]
  A <- A.simu.2[i-1]
  M <- M.simu.2[i-1] 
  I <- I.simu.2[i-1]
  
  if(tset.2[i] %in% t_all_storm_set.1){
    t_sen_set[i]<-0
  }#any time a storm occurs, t_sen gets set to 0 (re-starting the synchronized senescence function a_t)
  
  if (t_sen_set[i] < delay){
    a_t[i] <- 0#the rate of synchronized senescence is 0 for the first 150 days after the most recent storm (delay=150)
  } else{
    a_t[i] <- a*(sin(2*pi/period*(t_sen_set[i]+d_syn))+1)#after delay days since the most recent storm, synchronized senescence is initiated
  }
  
  
  if (t_sen_set[i] < delay){
    f_syn[i] <- 1-psi.simu.2[i-1]#during the delay period following a storm, the proportion of synchronized fronds is equal to the prop. that got removed (1-psi_A) and are growing back
  } else{
    f_syn[i] <- (1-psi.simu.2[i-1])*exp(-c*(t_sen_set[i]-delay))#after the delay period, the proportion of synchronized fronds decays exponentially
  }
  
  
  sA <- a_t[i]*f_syn[i]+b*(1-f_syn[i]) #overall senescence rate is a composite of the senescence of synchronous fronds (f_syn) and asynchronous fronds (1-f_syn) 
  
  if(tset.2[i] %in% t_mild_storm_set){ #when a weak storm occurs
    psi.simu.2[i]<-psi.A.mild #psi_A=psi.A.mild 
    A.simu.2[i]<-psi.A.mild*A #A=psi.A.mild*A
    J.simu.2[i]<-1*J#J, M, and I are not affected by mild storms
    M.simu.2[i]<-1*M
    I.simu.2[i]<-1*I
    
  } else if (tset.2[i] %in% t_severe_storm_set.1){ #when a severe storm occurs
    psi.simu.2[i]<-psi.A.sev #psi_A=psi.A.sev
    A.simu.2[i]<-psi.A.sev*A #A=psi.A.sev*A
    J.simu.2[i]<-psi.b*J#the proportion of benthic organisms (J, M, and I) remaining after a storm is determined by psi.b
    M.simu.2[i]<-psi.b*M
    I.simu.2[i]<-psi.b*I
    
  } else {
    psi.simu.2[i]<-psi.simu.2[i-1]
    dJ <- (r_G*L_surface*exp(-k_l*A)*G-r_J*L_surface*exp(-k_l*A)*J-m_J*J^2)*dt #equation for change in juvenile giant kelp sporophyte population
    dA <- (r_J*L_surface*exp(-k_l*A)*J+g_A*A*L_surface*((K_A-A)/K_A)-sA*A)*dt#equation for change in adult giant kelp frond density
    dM <- (g_M*M*L_surface*exp(-k_l*A)*((S_T-M-alpha*I)/S_T)-s_M*M+sigma_M)*dt#equation for change in understory macroalgae cover
    dI <- (g_I*I*((S_T-I-beta*M)/S_T)-s_I*I+sigma_I)*dt#equation for change in sessile invert cover
    J.simu.2[i] <- J + dJ
    A.simu.2[i] <- A + dA
    M.simu.2[i] <- M + dM
    I.simu.2[i] <- I + dI}
  
  dG <- (sigma_ext+sigma_A*A-r_G*L_surface*exp(-k_l*A)*G-m_G*G^2)*dt#equation for change in giant kelp gametophyte population
  G.simu.2[i] <- G + dG
  sA_set.1[i] <- sA 
}

# plot results
plot(x = tset.2, y = A.simu.2, col = kelp_col, type = "l", lwd = 1.5, xlab = "time (years)", ylab = "Relative abundance", xaxt = "n")
axis(side = 1, at = seq(from = 0, to = 4000, by = 365), labels = seq(from = 0, to = 10, by = 1))
lines(x = tset.2, y = I.simu.2, col = si_col, type = "l", lwd = 1.5)
lines(x = tset.2, y = M.simu.2, col = ua_col, type = "l", lwd = 1.5)
abline(h=A.simu.2[length(A.simu.2)], lty = 2)# add line for kelp abundance at the end of the simulation
abline(v = max(tset.2), lty = 2) # add line at the final timepoint

# ------------------
# high variability site
# ------------------

psi.A.mild<-.95#proportion of adult giant kelp fronds remaining after a mild storm
psi.A.sev<-0#proportion of adult giant kelp fronds remaining after a severe storm (change this value to change intensity of severe storms)
psi.b<-1#proportion of J, M, and I remaining after a severe storm (here no scouring)

#mild storms
t_mild_storm_set <- NaN # no mild storms

t_severe_storm_set.2<-c(365, 730, 1095, 2190, 2555, 2700, 2920, 3285, 3400, 3760) -5 # time points at which severe storms occur

tset.3 <- unique(sort(c(t_severe_storm_set.2, sim_tpts.1)))#new vector of time steps after any dates with storms that weren't already in sim_tpts.1 have been incorporated

t_all_storm_set.2<-t_severe_storm_set.2#all time points at which a storm occurs


#creating holding vectors for the outputs of the for loop; setting the first element in each holding vector
a_t <- NaN*tset.3;a_t[1]<-0#senescence rate of synchronized fronds
t_sen_set<-NaN*tset.3;t_sen_set[1]<-0#t_sen is a separate time vector; a_t is a function of t_sen
sA_set.2 <- NaN*tset.3; sA_set.2[1]<- 0#overall rate of frond senescence
f_syn <- NaN*tset.3 #proportion of the fronds that are synchronized 
G.simu.3 <- NaN*tset.3; G.simu.3[1]<-G0 #giant kelp gametophytes
J.simu.3 <- NaN*tset.3; J.simu.3[1]<-J0 #juvenile giant kelp sporophytes
A.simu.3 <- NaN*tset.3; A.simu.3[1]<-Astart #adult giant kelp fronds (choose A0, I0, M0 so model starts near dynamic equilibrium)
M.simu.3 <- NaN*tset.3; M.simu.3[1]<-Mstart #understory macroalgae 
I.simu.3 <- NaN*tset.3; I.simu.3[1] <- Istart #sessile invertebrates
psi.simu.3<-NaN*tset.3; psi.simu.3[1:2]<-1 #psi_A (proportion of giant kelp fronds remaining after a storm)


#simulate the model using a for loop
for(i in 2:length(tset.3)){
  dt <- tset.3[i]-tset.3[i-1]
  t_sen_set[i]<-t_sen_set[i-1]+dt
  G <- G.simu.3[i-1]
  J <- J.simu.3[i-1]
  A <- A.simu.3[i-1]
  M <- M.simu.3[i-1] 
  I <- I.simu.3[i-1]
  
  if(tset.3[i] %in% t_all_storm_set.2){
    t_sen_set[i]<-0
  }#any time a storm occurs, t_sen gets set to 0 (re-starting the synchronized senescence function a_t)
  
  if (t_sen_set[i] < delay){
    a_t[i] <- 0#the rate of synchronized senescence is 0 for the first 150 days after the most recent storm (delay=150)
  } else{
    a_t[i] <- a*(sin(2*pi/period*(t_sen_set[i]+d_syn))+1)#after delay days since the most recent storm, synchronized senescence is initiated
  }
  
  
  if (t_sen_set[i] < delay){
    f_syn[i] <- 1-psi.simu.3[i-1]#during the delay period following a storm, the proportion of synchronized fronds is equal to the prop. that got removed (1-psi_A) and are growing back
  } else{
    f_syn[i] <- (1-psi.simu.3[i-1])*exp(-c*(t_sen_set[i]-delay))#after the delay period, the proportion of synchronized fronds decays exponentially
  }
  
  
  sA <- a_t[i]*f_syn[i]+b*(1-f_syn[i]) #overall senescence rate is a composite of the senescence of synchronous fronds (f_syn) and asynchronous fronds (1-f_syn) 
  
  if(tset.3[i] %in% t_mild_storm_set){ #when a weak storm occurs
    psi.simu.3[i]<-psi.A.mild #psi_A=psi.A.mild 
    A.simu.3[i]<-psi.A.mild*A #A=psi.A.mild*A
    J.simu.3[i]<-1*J#J, M, and I are not affected by mild storms
    M.simu.3[i]<-1*M
    I.simu.3[i]<-1*I
    
  } else if (tset.3[i] %in% t_severe_storm_set.2){ #when a severe storm occurs
    psi.simu.3[i]<-psi.A.sev #psi_A=psi.A.sev
    A.simu.3[i]<-psi.A.sev*A #A=psi.A.sev*A
    J.simu.3[i]<-psi.b*J#the proportion of benthic organisms (J, M, and I) remaining after a storm is determined by psi.b
    M.simu.3[i]<-psi.b*M
    I.simu.3[i]<-psi.b*I
    
  } else {
    psi.simu.3[i]<-psi.simu.3[i-1]
    dJ <- (r_G*L_surface*exp(-k_l*A)*G-r_J*L_surface*exp(-k_l*A)*J-m_J*J^2)*dt #equation for change in juvenile giant kelp sporophyte population
    dA <- (r_J*L_surface*exp(-k_l*A)*J+g_A*A*L_surface*((K_A-A)/K_A)-sA*A)*dt#equation for change in adult giant kelp frond density
    dM <- (g_M*M*L_surface*exp(-k_l*A)*((S_T-M-alpha*I)/S_T)-s_M*M+sigma_M)*dt#equation for change in understory macroalgae cover
    dI <- (g_I*I*((S_T-I-beta*M)/S_T)-s_I*I+sigma_I)*dt#equation for change in sessile invert cover
    J.simu.3[i] <- J + dJ
    A.simu.3[i] <- A + dA
    M.simu.3[i] <- M + dM
    I.simu.3[i] <- I + dI}
  
  dG <- (sigma_ext+sigma_A*A-r_G*L_surface*exp(-k_l*A)*G-m_G*G^2)*dt#equation for change in giant kelp gametophyte population
  G.simu.3[i] <- G + dG
  sA_set.2[i] <- sA 
}


# plot results
plot(x = tset.3, y = A.simu.3, col = kelp_col, type = "l", lwd = 1.5, xlab = "time (years)", ylab = "Relative abundance", xaxt = "n")
axis(side = 1, at = seq(from = 0, to = 4000, by = 365), labels = seq(from = 0, to = 10, by = 1))
lines(x = tset.3, y = I.simu.3, col = si_col, type = "l", lwd = 1.5)
lines(x = tset.3, y = M.simu.3, col = ua_col, type = "l", lwd = 1.5)
abline(h=A.simu.2[length(A.simu.2)], lty = 2)# add line for kelp abundance at the end of the low variability simulation
abline(v = max(tset.3), lty = 2) # add line at final time point

















