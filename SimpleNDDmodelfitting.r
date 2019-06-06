#SimpleNDDmodelfitting.r
# This program takes as inputs count data and returns the estimates of exonential, ricker, and theta-logistic 
# model fits to the data. NOTE: for this 
# program to work, there can't be any missing data in the time series

# this program is set up with the data set for the Bay Checkerspot butterfly, and also for the Sand Lizard 
# of Zion Natl. Park 

graphics.off()# closes fig window
rm(list=ls())  # clear all variables

#-----INPUT PARAMETERS----------------------
# requires that simpleext.functions.R be in the same home directory

# list of years that data were collected

# data for Bay Checkspot Butterfly
yrs=1960:1986
Ns=c(90,
175,
40,
45,
175,
200,
425,
425,
800,
256,
713,
198,
1819,
575,
567,
1819,
7227,
852,
216,
244,
267,
1753,
999,
1788,
143,
79,
94)

#data for the Sand Lizard
yrs=1969:1979
Ns=c(166,
201,
180,
96,
92,
141,
114,
117,
175,
106,
144)

nn=length(yrs)

windows()
plot(yrs,Ns, type='b')
windows()
plot(yrs[2:nn], log(Ns[2:nn]/Ns[1:(nn-1)]))
windows()
plot(Ns[1:(nn-1)], log(Ns[2:nn]/Ns[1:(nn-1)]))

# note that in the final line of code, for the legend, you need to adjust the first two numbers to fit the 
# legend into the space of the figure. these are the x,y coordinates for the figure you are making

#---- END OF INPUTS ---------------------------
 numNs = length(Ns)
# first, estimate lambdas and the right quantities to do the regression:
lams = (Ns[2:numNs]/Ns[1:(numNs-1)])^(1/(yrs[2:numNs]-yrs[1:(numNs-1)]  ))
#xx=sqrt(yrs[2:numNs]-yrs[1:(numNs-1)])
#yy=log(Ns[2:9]/Ns[1:(numNs-1)])/xx

outs = matrix(0,3,5) # matrix to hold the outputs for model comparison and simulation
			   # this will hold: r,K,theta,residual variance, parameter number

# ----------------fit exponential model-------------------------
xx=matrix(1,length(lams),1)
regressout=lm((log(lams)~xx))
sumregress=summary(regressout)
anovaregress=anova(regressout)

outs[1,1]=sumregress$coefficients[1] # estimate of r from the regression
outs[1,4]=anovaregress$Sum/length(lams)	# estimate of residual var from the regression
outs[1,5]=2 #number of parameters

sumregress
outs[1,1]
outs[1,4]
# ----------------fit ricker model-------------------------
xx= Ns[-length(Ns)]
regressout=lm(log(lams)~xx)
sumregress=summary(regressout)
anovaregress=anova(regressout)


outs[2,1]=sumregress$coefficients[1,1]# estimate of r from the regression
outs[2,2]=-sumregress$coefficients[1,1]/sumregress$coefficients[2,1] # estimate of K from the regression
outs[2,4]=anovaregress$Sum[2]/length(lams) # res. variance 
outs[2,5]=3 #number of parameters

#regressout=nls(log(lams)~r*(1-(xx/K)), start=list(r=-0.01,K=mean(Ns))  ) 
#, algorithm = "port",upper=c(r=10,K=10*max(Ns)),lower=c(r=-10,K=1)) # these add-ons can help with some 
# convergence problems


#outs[2,1]=sumregress$coefficients[1,1] # estimate of r from the regression
#outs[2,2]=sumregress$coefficients[2,1] # estimate of K from the regression

#outs[2,4]=sum( (sumregress$residuals)^2)/length(lams)

# ----------------fit theta logistic model using the R function nlm ("nonlinear minimization")-------------------------


 # First, define the "sum of squares" function to minimize
SStlog=function(p) { sum( (log(lams) - (p[1]*(1-(xx/p[2])^p[3])) )^2 ) }


out2=nlm(SStlog,c(outs[2,1],outs[2,2],1))	# perform the nonlinear minimization, using the Ricker values for r and K 
					# and theta=1 as the starting parameter values
btlog=out2$estimate		# retrieve the parameter estimates
Vt=out2$minimum/length(lams); 

outs[3,1:3]=btlog
outs[3,4]=Vt
outs[3,5]=4



# estimate and plot the predicted loglambdas vs Nstart values for each model:
Ex=matrix(outs[1,1],length(xx),1)
Rick=outs[2,1]*(1-(xx/outs[2,2]))
ThetaLog=outs[3,1]*(1-(xx/outs[3,2])^outs[3,3])



Y=cbind(log(lams),Ex,Rick,ThetaLog)


matplot(xx,Y,
	xlab="N(t)",ylab="log{N(t+1)/N(t)}",
	type=c("p","p","p","p"),pch=21,col=c("black","blue","green","red"),lty=1,lwd=2)

legend(50000,0.2,c("Data","DI model","Ricker","Theta Logistic"),pch=c(21,-1,-1,-1),
	col=c("black","blue","green","red"),lty=c("blank","solid","solid","solid"),lwd=2)

# now, get the AIC values for each model to compare them:

q = length(lams) # number of data points
LLs=-.5*q*(log(2*pi*outs[,4])+1)	# use Vd from the DI model to compute the log likelihood
AICs=-2*LLs + (2*outs[,5]*q/(q-outs[,5]-1))

modresults=matrix( c(outs),3,5,
dimnames=list(c("exp", "ricker", "theta"),
               c("r", "K", "theta", "res.variance", "num parameters")))
modresults
print("AIC values for the Exponential, Ricker, and Theta-Logistic Models")
print(AICs)

