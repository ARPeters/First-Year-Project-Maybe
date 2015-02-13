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

