rm(list = ls(all.names = TRUE))
library(survsim)
library(survival)
library(foreign)

#http://stats.stackexchange.com/questions/109237/how-to-generate-survival-data-with-time-dependent-covariates-using-r

#First set of simulated, time-invariant, censorship-free survival times
x<-list(c("bern", 0.5), c("bern", 0.5), c("bern", 0.5), c("bern", 0.5))
b<-list(-0.7, -0.7, -0.1, -0.1)

firstData<-simple.surv.sim(n=500, foltime=365, dist.ev="weibull", anc.ev<-1, beta0.ev=3, anc.cens=1, beta0.cens=500, beta=b, x=x)
firstData

ds<-firstData
colnames(ds)<-c("nid", "status", "start", "stop", "Z", "X1", "X2", "X3", "X4")

survobject<-Surv(time=ds$stop, ds$status==1)

coxphTest<-coxph(survobject~X1 + X2 + X3 + X4, data=ds)

weibullTest<-survreg(survobject ~ X1 + X2+ X3 + X4, data=ds)

str(weibullTest)
weibullTest$loglik[2]


#Functional, but apparently cannot accomodate interactions, or time-dependent covariates. So, different tactic:
#Bender et al, 2005
#Control Data
#X5 and X6 are interactions between 1/3 and 2/4
X1<-rbinom(5000, 1, 0.5)
X2<-rbinom(5000, 1, 0.5)
X3<-rbinom(5000, 1, 0.5)
X4<-rbinom(5000, 1, 0.5)
X5<-X1*X3
X6<-X2*X4
U<-runif(5000)
scale<-1
negUlog<--1*log(U)
betaVector<-c(0.7, 0.7, 0.1, 0.1, 0.7, 0.7)
hf<-scale*(exp(X1*betaVector[1]+X2*betaVector[2]+X3*betaVector[3]+X4*betaVector[4]+X5*betaVector[5]+X6*betaVector[6]))
T<-negUlog/hf
T<-T*100
status<-1
dsNew<-cbind(X1, X2, X3, X4, X5, X6, scale, U, negUlog, hf, T, status)
dsNew<-as.data.frame(dsNew)


survobject<-Surv(time=dsNew$T, event=dsNew$status)
coxphTest<-coxph(survobject~X1 + X2 + X3 + X4 + X5 + X6, data=dsNew)

weibullTest<-survreg(survobject~X1 + X2 + X3 + X4 + X5 + X6, data=dsNew, dist="weibull")

#Now, the time-dependent effects


