
# Alternative ways of estimating variability in time series. Based on The consecutive disparity index, D: a measure of temporal variability in ecological studies, Ecosphere. And code found online at https://github.com/T-Engel/CValternatives/blob/master/R/Functions.R

PV <- function (Z){
  n = length(Z)
  pairs = combn(Z,2)
  min_z = apply(pairs,2, min)
  max_z = apply(pairs,2, max)
  z = 1- (min_z/max_z)
  PV=2*sum(z)/(n*(n-1))
  return(PV)
}

# According to Heath 2006 Oikos, the function PV should be adapted to accomodate zeros in the dataset. Specifically, if zi = zj then z should equal zero, if zi != zj then z is 1-min(zi,zj)/max(zi,zj)

PV <- function (z){
  tryCatch(
    {
  Z = as.vector(z)
  n = length(Z)
  pairs = combn(Z,2)
  min_z = apply(pairs,2, min)
  max_z = apply(pairs,2, max)
  z = vector()
  for(i in 1:length(min_z)){
    if(min_z[i] == max_z[i]){
      z[i] = 0
    }else(
      z[i] = 1- (min_z[i]/max_z[i])
    )
  }
  PV=2*sum(z)/(n*(n-1))
  return(PV)
  },error = function(cond) {
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  })
}


D  <- function(p){
  P = as.vector(p)
  n = length(P)
  flixbus <- NA
  for(i in (1: (n-1))){
    flixbus[i]=abs(log(P[i+1]/P[i]))
  }
  D=sum(flixbus)/(n-1)
  return(D)
  
}

# Adapted function to include k offset based on < 1% of the mean
Dc <- function(p){
  P = as.vector(p)
  n = length(P)
  flixbus <- NA
  k <- mean(P, na.rm = T)*0.01 # according to FERNANDEZ-MARTINEZ et al. 2018 who proposed this metric, a constant can be added to both numerator and denominator. However, in order to compare between time series and not bias the variance estimate at low values, the constant should be some small fraction of the mean value (< 1% of the mean). They claim that changing the values by < 1% will change D < 1% which should be acceptible for ecological time series (see page 4).
  for(i in (1: (n-1))){
    flixbus[i]=abs(log((P[i+1]+k)/(P[i]+k)))
  }
  D=sum(flixbus)/(n-1)
  return(D)
  
}























