rm(list = ls(all.names = TRUE))
library(survsim)
library(survival)
library(foreign)
library(PermAlgo)

n=5000
m=365

xmat<-matrix(nrow=n*m, ncol=13)

xmat[,1]<-rep(round(rbinom(n,1, 0.5)), each=m)
xmat[,2]<-rep(round(rbinom(n,1, 0.5)), each=m)
xmat[,3]<-rep(round(rbinom(n,1, 0.5)), each=m)
xmat[,4]<-rep(round(rbinom(n,1, 0.5)), each=m)


#5th column is going to be column number/time interval as if we had applied counting process to a data set
xmat[,5]<-rep.int(1:m, times=n)

#6th and 7th columns are going to represent the log of time * Strong1 and Strong2, respectively.
xmat[,6]<-log(xmat[,5])*xmat[,1]
xmat[,7]<-log(xmat[,5])*xmat[,2]

#8th and 9th columns are going to represent time squared * Strong1 and Strong2, respectively
xmat[,8]<-xmat[,5]*xmat[,5]*xmat[,1]
xmat[,9]<-xmat[,5]*xmat[,5]*xmat[,2]

#10th and 11th columns are going to represent interactions between Strong1*Weak1 and Strong2*Weak2
xmat[,10]<-xmat[,1]*xmat[,3]
xmat[,11]<-xmat[,2]*xmat[,4]

#12 and 13th columns are going to represent interactions between heaviside function and strong1/strong2 respectively
for(i in 1:length(xmat[,5])){
  xmat[i,12]<-ifelse(xmat[i,5]<=182, 0, 1)  
  xmat[i,13]<-ifelse(xmat[i,5]<=182, 0, 1)
}

dsMaster<-as.data.frame(xmat)

colnames(dsMaster)<-c("Strong1", "Strong2", "Weak1", "Weak2", "intForFt", "logtStrong1", "logtStrong2", "t2Strong1", "t2Strong2", "Strong1Weak1", "Strong2Weak2", "Strong1H", "Strong2H")

#Now permute it!
dsLogt<-as.matrix(dsMaster[,c(1:4, 6, 7)])
eventTimesMaybe<-runif(n, 1, m)

dataLogt<-permalgorithm(n, m, Xmat=dsLogt, XmatNames=c("Strong1", "Strong2", "Weak1", "Weak2", "logtStrong1", "logtStrong2"), eventRandom=eventTimesMaybe, betas=c(0.7, 0.7, 0.1, 0.1, 0.7, 0.7))

attach(dataLogt)
survobjectlogt<-Surv(time=Start, time2=Stop, Event==1)

test<-coxph(survobjectlogt ~ Strong1 + Strong2 + Weak1 + Weak2 + logtStrong1 + logtStrong2, data=dataLogt, ties="breslow")
test
detach(dataLogt)

#Strong and Weak variables, plust two tsquared elements
#dst2<-as.matrix(dsMaster[,c(1:4, 8, 9)])
#eventTimesMaybe<-runif(n, 1, m)

#datat2<-permalgorithm(n, m, Xmat=dst2, XmatNames=c("Strong1", "Strong2", "Weak1", "Weak2", "t2Strong1", "t2Strong2"), eventRandom=eventTimesMaybe, betas=c(0.7, 0.7, 0.1, 0.1, 0.7, 0.7))

#attach(datat2)
#survobjectt2<-Surv(time=Start, time2=Stop, Event==1)

#testt2<-coxph(survobjectt2 ~ Strong1 + Strong2 + Weak1 + Weak2 + Strong1Weak1 + Strong2Weak2, data=datat2, ties="breslow")
#testt2
#detach(datat2)


#Strong and Weak variables, plust two interaction variables
dsControl<-as.matrix(dsMaster[,c(1:4, 10, 11)])
eventTimesMaybe<-runif(n, 1, m)

dataControl<-permalgorithm(n, m, Xmat=dsControl, XmatNames=c("Strong1", "Strong2", "Weak1", "Weak2", "Strong1Weak1", "Strong2Weak2"), eventRandom=eventTimesMaybe, betas=c(0.7, 0.7, 0.1, 0.1, 0.7, 0.7))

attach(dataControl)
survobjectControl<-Surv(time=Start, time2=Stop, Event==1)

testControl<-coxph(survobjectControl ~ Strong1 + Strong2 + Weak1 + Weak2 + Strong1Weak1 + Strong2Weak2, data=dataControl, ties="breslow")
testControl
detach(dataControl)

#Strong and Weak variables, plust two heaviside function variables
dsH<-as.matrix(dsMaster[,c(1:4, 12, 13)])
eventTimesMaybe<-runif(n, 1, m)

dataH<-permalgorithm(n, m, Xmat=dsH, XmatNames=c("Strong1", "Strong2", "Weak1", "Weak2", "Strong1H", "Strong2H"), eventRandom=eventTimesMaybe, betas=c(0.7, 0.7, 0.1, 0.1, 0.7, 0.7))

attach(dataH)
survobjectH<-Surv(time=Start, time2=Stop, Event==1)

testH<-coxph(survobjectH ~ Strong1 + Strong2 + Weak1 + Weak2 + Strong1H + Strong2H, data=dataH, ties="breslow")
testH
head(dataH)
Model1<-lm(dataH$Stop~Strong1 + Strong2 + Weak1 + Weak2 + Strong1H + Strong2H, data=dataH)

detach(dataH)
