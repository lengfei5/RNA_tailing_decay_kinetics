########################################################
########################################################
# Section :
#  this is simulation function by Annamaria
########################################################
########################################################
library(deSolve)
library(minpack.lm)
#library(viridis) #allows you to use viridis, magma, ... color schemes

reactionrates=function(t, c, Parms){
  
  #define r as a vector of zeroes, has the length of the number of required differential equations
  r=rep(0, length(c))
  
  #in the following we replace the zeroes in r with differential equations
  #differential equation(DE) of substrate (U0) consumption:
  r[1]= - Parms[1]*c["U0"]
  #DEs of first till n-1th tailing events:
  for(i in 2:length(r)){
    r[i]=Parms[i-1]*c[paste("U", i-2, sep="")]-Parms[i]*c[paste("U" ,i-1, sep="")]
  }
  
  #DE of last tailing event
  #r[length(r)-1]= Parms[[ncol(Data)-2]]*c[paste("U", ncol(Data)-3, sep="")]
  
  return(list(r))
  #the resulting DEs should be:
  #  r[1]=-k1*c["U0"]
  #  r[2]=k1*c["U0"]-k2*c["U1"]
  #  r[3]=k2*c["U1"]-k3*c["U2"]
  #  r[4]=k3*c["U2"]-k4*c["U3"]
  #  r[5]=k4*c["U3"]-k5*c["U4"]
  #  r[6]=k5*c["U4"]-k6*c["U5"]
  #  r[7]=k6*c["U5"] - K7*c["U6"]
}

##########################################
# define the error function for fitting using optim
##########################################
#Define difference as function to calculate the difference between model and data for least square regression:
f2min=function(pars.init, times, dat, cinit){
  #cinit=ICs
  #t=Data[,1]
  #solve ODE as before (out)
  model=ode(y=cinit, times=times, func=reactionrates, parms=pars.init)
  #calculate difference between model and data:
  error2 <- sum((model[, -1]-dat)^2)
  return(error2)
}


