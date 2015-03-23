rm(list = ls(all.names = TRUE))
library(cvTools)
library(survsim)
library(survival)
library(foreign)
library(PermAlgo)
library(ROCR)

#betas<-c(0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0, 0, ftw[l], ftw[l], 0, 0, 0, 0)

#Beta vector for f(t) = log(t)
#betas<-c(0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0.7, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

#Beta vector for f(t) = t
betas<-c(0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 0.0, 0.0)

#Beta vector for f(t) = interaction
#betas<-c(0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0)

#Beta vector for f(t) = heaviside
#betas<-c(0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 0.7)

#Betas for playing around with ordering of variables
#betas<-c(0.7, 0.7, 0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

#Creating a table of AICs and BIC values
reps<-5

fitTable<-data.frame(matrix(ncol=20, nrow=reps, ))
colnames(fitTable)<-c("AICC",  "AICH", "AICI", "AICLog", "AICT", "BICC","BICH", "BICI","BICLog", "BICT" , "cvC", "cvH", "cvI", "cvLog", "cvT", "AUCC", "AUCH", "AUCI", "AUCLog", "AUCT")

#For this weight of specific time function, create this many sets of data  
for(i in 1:reps){
  
  n=5000
  m=24
  mhalf<-m/2
  
  xmat<-matrix(nrow=n*m, ncol=15)
  
  xmat[,1]<-rep(round(rbinom(n,1, 0.5)), each=m)
  xmat[,2]<-rep(round(rbinom(n,1, 0.5)), each=m)
  xmat[,3]<-rep(round(rbinom(n,1, 0.5)), each=m)
  xmat[,4]<-rep(round(rbinom(n,1, 0.5)), each=m)
  xmat[,5]<-rep(round(rbinom(n,1, 0.5)), each=m)
  xmat[,6]<-rep(round(rbinom(n,1, 0.5)), each=m)
  
  #7th column is going to be column number/time interval as if we had applied counting process to a data set
  xmat[,7]<-rep.int(1:m, times=n)
  
  #8th and 9th columns are going to represent the log of time * Strong1 and Strong2, respectively.
  xmat[,8]<-log(xmat[,7])*xmat[,1]
  xmat[,9]<-log(xmat[,7])*xmat[,2]
  
  #10th and 11th columns are going to represent time * Strong1 and Strong2, respectively
  xmat[,10]<-(xmat[,7]/m)*xmat[,1]
  xmat[,11]<-(xmat[,7]/m)*xmat[,2]
  
  #12th and 13th columns are going to represent interactions between Strong1*Weak1 and Strong2*Weak2
  xmat[,12]<-xmat[,1]*xmat[,4]
  xmat[,13]<-xmat[,2]*xmat[,5]
  
  #14th and 15th columns are going to represent interactions between heaviside function and strong1/strong2 respectively
  for(j in 1:length(xmat[,7])){
    xmat[j,14]<-ifelse(xmat[j,7]<=mhalf, 0, xmat[j,1])  
    xmat[j,15]<-ifelse(xmat[j,7]<=mhalf, 0, xmat[j,2])
  }
  
  dsMaster<-as.data.frame(xmat)
  
  colnames(dsMaster)<-c("Strong1", "Strong2", "Strong3", "Weak1", "Weak2", "Weak3", "intForFt", "logtStrong1", "logtStrong2", "tStrong1", "tStrong2", "Strong1Weak1", "Strong2Weak2", "Strong1H", "Strong2H")
  
  #dsMasterNewT<-dsMaster[,c(10,11,1:6,7,8:9,12:13,14:15)]
  
  #Strong and Weak variables, plust time function variables
  dsT<-as.matrix(dsMaster[,c(1:6, 8:15)])
  #dsT<-as.matrix(dsMasterNewT[,c(1:8, 10:15)])
  
  eventTimesMaybe<-as.integer(runif(n, 1, m))
  
  dataT<-permalgorithm(n, m, Xmat=dsT, XmatNames=c("Strong1", "Strong2", "Strong3", "Weak1", "Weak2", "Weak3", "logtStrong1", "logtStrong2", "tStrong1", "tStrong2", "Strong1Weak1", "Strong2Weak2", "Strong1H", "Strong2H"), betas=betas, eventRandom=eventTimesMaybe)
  #dataT<-permalgorithm(n, m, Xmat=dsT, XmatNames=c("tStrong1", "tStrong2", "Strong1", "Strong2", "Strong3", "Weak1", "Weak2", "Weak3", "logtStrong1", "logtStrong2", "Strong1Weak1", "Strong2Weak2", "Strong1H", "Strong2H"), betas=betas, groupByD=TRUE)
  
  attach(dataT)
  survobjectT<-Surv(time=Start, time2=Stop, Event==1)
  
  #Creating the cox models for each set of variables
  
  testControl<-coxph(survobjectT ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3, data=dataT, ties="breslow")
  
  testH<-coxph(survobjectT ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3 + Strong1H + Strong2H, data=dataT, ties="breslow")
  
  testInt<-coxph(survobjectT ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3 + Strong1Weak1 + Strong2Weak2, data=dataT, ties="breslow")
  
  testLog<-coxph(survobjectT ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3 + logtStrong1 + logtStrong2, data=dataT, ties="breslow")
  
  testT<-coxph(survobjectT ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3 + tStrong1 + tStrong2, data=dataT, ties="breslow")
  #testT<-coxph(survobjectT ~  tStrong1 + tStrong2 + Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3, data=dataT, ties="breslow")
  print(testT)
  #Getting cv error for each model
  cvF<-cvFolds(n=length(dataT$Event), K=10)
  cvTestH<-cvTool(call=testH, y=dataT$Event, folds=cvF)
  cvTestI<-cvTool(call=testInt, y=dataT$Event, folds=cvF)
  cvTestLog<-cvTool(call=testLog, y=dataT$Event, folds=cvF)
  cvTestC<-cvTool(call=testControl, y=dataT$Event, folds=cvF)
  cvTestT<-cvTool(call=testT, y=dataT$Event, folds=cvF)
  
  predTestControl<-prediction(testControl$linear.predictors, dataT$Event)
  perfTestControl<-performance(predTestControl, measure="auc")
  AUCC<-as.numeric(perfTestControl@y.values)
  
  predTestH<-prediction(testH$linear.predictors, dataT$Event)
  perfTestH<-performance(predTestH, measure="auc")
  AUCH<-as.numeric(perfTestH@y.values)
  
  predTestI<-prediction(testInt$linear.predictors, dataT$Event)
  perfTestI<-performance(predTestI, measure="auc")
  AUCI<-as.numeric(perfTestI@y.values)
  
  predTestLog<-prediction(testLog$linear.predictors, dataT$Event)
  perfTestLog<-performance(predTestLog, measure="auc")
  AUCLog<-as.numeric(perfTestLog@y.values)
  
  predTestT<-prediction(testT$linear.predictors, dataT$Event)
  perfTestT<-performance(predTestT, measure="auc")
  AUCT<-as.numeric(perfTestT@y.values)
  detach(dataT)
}
