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

source("solution_kinetic_model_for_tailing.R")
##########################################
# import adn prepare data table
##########################################
A <- read.csv(file="R_file_TAAG.csv", header=TRUE, sep=";")
B <- read.csv(file="Final_Tailing_RNAseq.csv", header=TRUE, sep=";")

Data <- data.frame(B, stringsAsFactors = FALSE)
print(Data)
#print(Data[,2:ncol(Data)])
head(Data)

##########################################
# specify parameters and initia conditions
##########################################
times = Data$Time;
Parms <- rep(1, ncol(Data)-1) # simply just a vector

cinit = c(1, rep(0, ncol(Data)-2))
names(cinit) = paste0("U", c(0:(length(cinit)-1)))  

## test the ODE solution
out=ode(y=cinit, times=times, func=reactionrates, parms=Parms)
print(out)
par(mfrow = c(1, 1))
modeltimes <- seq(from=0, to=20, by=0.1)
plotmodel=ode(y=cinit, times=modeltimes, func=reactionrates, parms=Parms)
matplot(modeltimes, plotmodel[,2:ncol(plotmodel)], type="l", lty="solid", lwd=2, pch=19, col=rainbow(ncol(Data)-1))
legend("right", c(names(cinit)), fill=rainbow(ncol(Data)-1))

##########################################
# initialize argument for optim function
# then run optim for fitting the data to the model
# at the end get estimated parameters
# due to large number of parameters, the fitting step can be quite slow
# one solution is to chose good starting points of parameters, which is tested below 
##########################################
#pars.init = rep(1, length(Parms));
K0.warmup= abs(lm(log(Data[, 2]) ~ times - 1 )$coefficients) # fit the log(U0) with linear regression and get the first guess for K0
par.init = rep(K0.warmup, length(Parms))
lower.bounds = rep(0, length(Parms)) # lower boundaries for parameters
upper.bounds = rep(K0.warmup*10, length(Parms))

fit = optim(pars.init, f2min, times = times, dat = Data[, -1], cinit = cinit,
            method = 'L-BFGS-B', lower = lower.bounds, upper = upper.bounds,
            control = list(maxit=500,trace=0), hessian = FALSE)

#summary()
Ks <- fit$par
print(Ks)

#plot the Data and the fit
#here it would actually be good to INCREASE THE t RESOLUTION TO MAKE PLOT LESS "EDGY" !!!
optimized=ode(y=cinit, times=modeltimes, func=reactionrates, parms=Ks)

par(mfrow = c(3, 4))
cols = rainbow((ncol(Data)))

for(n in 2:ncol(Data)){
  plot(times, Data[,n], pch=19, col=cols[n], main = colnames(Data)[n], ylim = range(c(Data[,n], optimized[, n])))
  points(modeltimes, optimized[, n], lty="solid", pch=19, col=cols[n], type = "l")
}


#to plot each single events
df <- as.data.frame(Data)
library(tidyverse)
library(readxl)
library(cowplot)
df %>% gather(tail, value, Substrate, starts_with("Tail"), factor_key = TRUE) %>% ggplot() + geom_point(aes(x = Time, y = value)) + facet_wrap(~tail)
