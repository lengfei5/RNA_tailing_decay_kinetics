#The following script is designed to model uridylation tailing events in a chain-reaction.

#Basically, new units of "U" are added to a Substrate (U0) in the first reaction:
#U0 + U -> U1       [reaction rate k1]

#The following reactions then add further U to the resulting products:
#U1 + U -> U2       [reaction rate k2]
#U2 + U -> U3       [reaction rate k3]
#U3...              [reaction rate k...]

#For the analysis here we assume that the reaction is pseudo-first order 
#(U is in excess and the enyzme catalyzing the reaction is constant over time)
#and that addition of U is irreversible

#REQUIRES is a .csv file with headers (!). The first column should contain time values,
#the second column substrate (U0) concentrations, and the  following
#columns concentrations of "tailed" products (U1, U2, U3, ...)

#The script was based on the code of the following post:
#https://notesofdabbler.wordpress.com/2013/06/30/learning-r-parameter-fitting-for-models-involving-differential-equations/

# (January 2019)


library(deSolve)
library(minpack.lm)
#library(viridis) #allows you to use viridis, magma, ... color schemes

#read your .csv file into variable A and save as matrix called Data;
A <- read.csv(file="/Users/annamaria.sgromo/Desktop/R/R_file_TAAG.csv", header=TRUE, sep=";")
Data <- as.matrix(A)
print(Data)
print(Data[,2:ncol(Data)])

#define Parms as list of rates k; here we define k1 to k... = 1
Parms <- as.list(rep(1, ncol(Data)-2))

reactionrates=function(Time, c, Parms){
  #define r as a vector of zeroes, has the length of the number of required differential equations
  r=rep(0, ncol(Data)-1)
  
  #in the following we replace the zeroes in r with differential equations
  #differential equation(DE) of substrate (U0) consumption:
  r[1]=-Parms[[1]]*c["U0"]
  
  #DEs of first till n-1th tailing events:
  for(i in 2:(ncol(Data)-2)){
    r[i]=Parms[[i-1]]*c[paste("U", i-2, sep="")]-Parms[[i]]*c[paste("U" ,i-1, sep="")]
  }

  #DE of last tailing event
  r[ncol(Data)-1]=Parms[[ncol(Data)-2]]*c[paste("U", ncol(Data)-3, sep="")]
  
  return(list(r))
  #the resulting DEs should be:
  #  r[1]=-k1*c["U0"]
  #  r[2]=k1*c["U0"]-k2*c["U1"]
  #  r[3]=k2*c["U1"]-k3*c["U2"]
  #  r[4]=k3*c["U2"]-k4*c["U3"]
  #  r[5]=k4*c["U3"]-k5*c["U4"]
  #  r[6]=k5*c["U4"]-k6*c["U5"]
  #  r[7]=k6*c["U5"]
}

#define initial concentrations as list specifying the initial concentration of each species
#here we define U0=1; all others (U1, U2, U3, ...) are 0.
Initial_concentrations <-list()
Name <- "U0"
Initial_concentrations[[Name]] <- 1
for(i in 1:(ncol(Data)-2)){
  Name <- paste("U", i, sep="")
  Initial_concentrations[[Name]] <- 0
}

#to put Initial concentrations into ODE solver they need to be in vector format. ICs now defined as initial concentrations in vector format:
ICs <- unlist(Initial_concentrations)

#solving ODEs to get predicted concentrations (c in reactionrates) with guessed Parms; out  is the solution to the ODEs (matrix format)
cinit=ICs
t=Data[,1]
out=ode(y=cinit, times=t, func=reactionrates, parms=Parms)
#head(out)
print(out)

#The following is a solution of the ODEs with high time resolution for plotting of the model with guessed parameters
#note that the solution to the ODEs is almost the same as above for out (different times here)
modeltimes <- seq(from=0, to=20, by=0.1)
plotmodel=ode(y=cinit, times=modeltimes, func=reactionrates, parms=Parms)
matplot(modeltimes, plotmodel[,2:ncol(plotmodel)], type="l", lty="solid", lwd=2, pch=19, col=rainbow(ncol(Data)-1))
legend("right", c(names(Initial_concentrations)), fill=rainbow(ncol(Data)-1))

#Define difference as function to calculate the difference between model and data for least square regression:
Difference=function(Parms){
  cinit=ICs
  t=Data[,1]
  #solve ODE as before (out)
  model=ode(y=cinit, times=t, func=reactionrates, parms=Parms)
  #calculate difference between model and data:
  Diff <- matrix(model-Data)
  return(Diff)
}

#Nonlinear Least Square (NLS) regression using Levenberg-Marquardt (LM) algorithm to fit model to data; define Ks as list of optimized rates:
fit=nls.lm(par=Parms, fn=Difference)
summary(fit)
Ks <- as.list(coef(fit))
print(Ks)

#calculate fitted model by solving ODE with Ks, see also function "reactionrates":
plotfit=function(Time, c, Parms){
  f=rep(0, ncol(Data)-1)
  
  f[1]=-Ks[[1]]*c["U0"]
  
  for(i in 2:(ncol(Data)-2)){
    f[i]=Ks[[i-1]]*c[paste("U", i-2, sep="")]-Ks[[i]]*c[paste("U" ,i-1, sep="")]
  }
  
  f[ncol(Data)-1]=Ks[[ncol(Data)-2]]*c[paste("U", ncol(Data)-3, sep="")]
  
  return(list(f))
}

cinit=ICs
t=Data[,1]

#plot the Data and the fit
#here it would actually be good to INCREASE THE t RESOLUTION TO MAKE PLOT LESS "EDGY" !!!
optimized=ode(y=cinit, times=t, func=plotfit, parms=Parms)
matplot(Data[,1], Data[,2:ncol(Data)], pch=19, col=rainbow(ncol(Data)-1))
matlines(Data[,1], optimized[,2:ncol(optimized)], lty="solid", pch=19, col=rainbow(ncol(Data)-1))

#to plot each single events
df <- as.data.frame(Data)
library(tidyverse)
library(readxl)
library(cowplot)
df %>% gather(tail, value, Substrate, starts_with("Tail"), factor_key = TRUE) %>% ggplot() + geom_point(aes(x = Time, y = value)) + facet_wrap(~tail)
