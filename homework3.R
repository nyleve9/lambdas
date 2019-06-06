graphics.off()# closes fig window
rm(list=ls())  # clear all variables

#-----INPUT PARAMETERS----------------------
# list of years that data were collected
yrs=c(1972,
1973,
1974,
1975,
1976,
1977,
1978,
1979,
1980,
1981,
1982,
1983,
1984,
1985,
1986,
1987)


#list of population sizes
Ns=c(978,
2339,
553,
306,
71,
227,
395,
1168,
274,
165,
203,
101,
36,
118,
300,
536
)

l=length(Ns)
l
lambdas<-((Ns[2:l])/(Ns[1:(l-1)]))
lambdas
Zlams=cbind(yrs,lambdas)
Zlams
plot(Zlams,type="o")

hist((log(lambdas)))

lams2=cbind(Ns,lambdas)
plot(lams2,type="p")

#---- END OF INPUTS ---------------------------

# first, estimate lambdas and the right quantities to do the regression:
numNs = length(Ns)
lams = (Ns[2:numNs]/Ns[1:(numNs-1)])^(1/(yrs[2:numNs]-yrs[1:(numNs-1)]  ))

xx=sqrt(yrs[2:numNs]-yrs[1:(numNs-1)])
yy=log(Ns[2:numNs]/Ns[1:(numNs-1)])/xx

regressout=lm(yy~0+xx) # this conducts the linear regression that is described in the text, and in lecture using the fairybell example
sumregress=summary(regressout)
anovaregress=anova(regressout)


#mu<-mean(log(lambdas))
#mu
#sigsq<-var(log(lambdas))
#sigsq
#geomean<-exp(mu)
#geomean

mu=sumregress$coefficients[1] # estimate of mu from the regression

sig2=anovaregress$Mean[2]	# estimate of sigma^2 from the regression

mu
sig2
sumregress
anovaregress

#begin extcdf function from homework2:

#mu<-mean(log(lambdas))
#mu
#sig2<-var(log(lambdas))
#sig2
#geomean<-exp(mu)
#geomean
mu
sig2
d<-log(Ns[16])-log(100)
d
tmax<-100
#making an extinction CDF plot
extcdf=function(mu,sig2,d,tmax) {
  t=1:tmax
  G=pnorm((-d-mu*t)/sqrt(sig2*t))+
    exp(-2*mu*d/sig2)*pnorm((-d+mu*t)/sqrt(sig2*t))
  return(G)
}
sample=extcdf(mu,sig2,d,tmax)
plot(sample,ylab="Cumulative probability of extinction",xlab="Time (Years)")

#-------------------------------------------------------
#SimpleGrowthChooser.r
# This program allow you to simulate a 
# very simple growth process, with you providing a discrete
# list of lambda values and probabilities that each will occur

graphics.off()# closes fig window
#rm(list=ls())  # clear all variables

# -----INPUT PARAMETERS----------------------
#lambdas = c(0.4,1.5,0.3,4.8,0.5,4.1,1.6,1.1,0.9,0.8,0.7,3.4,1.2,0.1,2.3,4.4,0.2,1.3,2.5,0.8) # list of lambdas that can occur
lams=lambdas
lams

maxyr = 100 # the number of years to simulate

startN = 536 # starting population size

Reps = 1000 # number of replicate runs to make  	

Nqe = 100 # quasi-extinction threshold

maxcap = 1000
#---- END OF INPUTS ---------------------------

Ns = matrix(0,maxyr,Reps) # set up a matrix for the pop sizes to fill


for (jj in 1:Reps) {
  Ns[1,jj]=startN # initialize with starting pop size
  for (ii in 1:(maxyr-1)) {
    lam_t=sample(lams,1) # choose a random lambda value
    
    
    Ns[(ii+1),jj]=Ns[ii,jj]*lam_t # this grows the population one year
    
    #enforcing the quasi-extinction threshold:
    if (Ns[(ii+1),jj] <= Nqe) {Ns[(ii+1),jj]=0} 
    #enforcing the maximum cap on population numbers:
    if (Ns[(ii+1),jj] >= maxcap)  {Ns[(ii+1),jj]=maxcap}
  } # end of ii loop
} # end of jj loop


#perparing a matrix to use in estimating extinction risk:
Ns2=Ns 
Ns2[Ns2 >0] = 1 #set all values of Ns that are greater than one to 
alive = apply(Ns2,1,sum)
dead = 100*(1-alive/Reps)


allyrs = c(1:maxyr)
windows()

plot(allyrs,dead,type = "l",xlab="Year",ylab="Extinction CDF")

windows()
matplot(allyrs,Ns, type = "l",xlab="Year",ylab="Nt",main="Arithmetic scale")

windows()
matplot(allyrs,log(Ns), type = "l",xlab="Year",ylab="log(Nt)",main="Log scale")
