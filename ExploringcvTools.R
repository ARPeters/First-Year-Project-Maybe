rm(list = ls(all.names = TRUE))

#Exploring cross validation
#install.packages("cvTools")
library(cvTools)
library(PermAlgo)
library(survival)
library(ROCR)
betas<-c(0.7, 0.7, 0.1, 0.1, 0, 0, 0, 0, 0, 0, 0.7, 0.7)


n=1000
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
for(j in 1:length(xmat[,5])){
  xmat[j,12]<-ifelse(xmat[j,5]<=182, 0, xmat[j,1])  
  xmat[j,13]<-ifelse(xmat[j,5]<=182, 0, xmat[j,2])
}

dsMaster<-as.data.frame(xmat)

colnames(dsMaster)<-c("Strong1", "Strong2", "Weak1", "Weak2", "intForFt", "logtStrong1", "logtStrong2", "t2Strong1", "t2Strong2", "Strong1Weak1", "Strong2Weak2", "Strong1H", "Strong2H")

#Strong and Weak variables, plust two heaviside function variables
dsH<-as.matrix(dsMaster[,c(1:4,6:13)])
eventTimesMaybe<-runif(n, 1, m)

head(dsH)
dataH<-permalgorithm(n, m, Xmat=dsH, XmatNames=c("Strong1", "Strong2", "Weak1", "Weak2", "logtStrong1", "logtStrong2", "t2Strong1", "t2Strong2", "Strong1Weak1", "Strong2Weak2", "Strong1H", "Strong2H"), eventRandom=eventTimesMaybe, betas=betas)

attach(dataH)
survobjectH<-Surv(time=Start, time2=Stop, Event==1)

testH<-coxph(survobjectH ~ Strong1 + Strong2 + Weak1 + Weak2 + Strong1H + Strong2H, data=dataH, ties="breslow")

testControl<-coxph(survobjectH ~ Strong1 + Strong2 + Weak1 + Weak2 + Strong1Weak1 + Strong2Weak2, data=dataH)

testLog<-coxph(survobjectH ~ Strong1 + Strong2 + Weak1 + Weak2 + logtStrong1 + logtStrong2, data=dataH)

testControlActual<-coxph(survobjectH ~ Strong1 + Strong2 + Weak1 + Weak2, data=dataH)

#Setting up the cvFolds
cvF<-cvFolds(n=length(dataH$Event), K=10)
cvTestH<-cvTool(call=testH, y=dataH$Event, folds=cvF)
cvTestControl<-cvTool(call=testControl, y=dataH$Event, folds=cvF)
cvTestLog<-cvTool(call=testLog, y=dataH$Event, folds=cvF)
cvTestControlActual<-cvTool(call=testControlActual, y=dataH$Event, folds=cvF)

detach(dataH)

cvTestH
cvTestControl
cvTestLog
cvTestControlActual

##Testing c statistic

predTestH<-prediction(testH$linear.predictors, dataH$Event)
perfTestH<-performance(predTestH, measure="auc")
as.numeric(perfTestH@y.values)

