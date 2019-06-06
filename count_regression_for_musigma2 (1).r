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
lams=cbind(yrs,lambdas)
lams
plot(lams,type="o")

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
d<-log(Ns[16])-log(120)
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


