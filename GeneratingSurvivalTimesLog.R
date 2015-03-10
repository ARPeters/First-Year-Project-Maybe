rm(list = ls(all.names = TRUE))
library(survsim)
library(survival)
library(foreign)
library(PermAlgo)


#ftw<-c(0, 0.01, 0.35, 0.7)
ftw<-c(1:70)/100
#ftw<-c(0:20)/5

AICPropCorrect<-vector(length=length(ftw))
AICPropWC<-vector(length=length(ftw))
AICPropWH<-vector(length=length(ftw))


BICPropCorrect<-vector(length=length(ftw))
BICPropWC<-vector(length=length(ftw))
BICPropWH<-vector(length=length(ftw))


#For: weight of a specific time function (heaviside, the last one)
for(l in 1:length(ftw)){
  #Declaring Betas: first, we are looking at the heaviside function, so the last two time-dependent 
  #variables are given betas, all else are given 0 weight. 
  betas<-c(0.7, 0.7, 0.1, 0.1, ftw[l], ftw[l], 0, 0, 0, 0, 0, 0)
  
  #Creating a table of AICs and BIC values
  reps<-10
  
  fitTable<-data.frame(matrix(ncol=6, nrow=reps, ))
  colnames(fitTable)<-c("AICH", "BICH", "AICC","BICC", "AICLog","BICLog")
  
  #For this weight of specific time function, create this many sets of data  
  for(i in 1:reps){
    
    n=500
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
    dsLog<-as.matrix(dsMaster[,c(1:4,6:13)])
    eventTimesMaybe<-runif(n, 1, m)
    
    dataLog<-permalgorithm(n, m, Xmat=dsLog, XmatNames=c("Strong1", "Strong2", "Weak1", "Weak2", "logtStrong1", "logtStrong2", "t2Strong1", "t2Strong2", "Strong1Weak1", "Strong2Weak2", "Strong1H", "Strong2H"), eventRandom=eventTimesMaybe, betas=betas)
    
    attach(dataLog)
    survobjectLog<-Surv(time=Start, time2=Stop, Event==1)
    
    testH<-coxph(survobjectLog ~ Strong1 + Strong2 + Weak1 + Weak2 + Strong1H + Strong2H, data=dataLog, ties="breslow")
    
    testControl<-coxph(survobjectLog ~ Strong1 + Strong2 + Weak1 + Weak2 + Strong1Weak1 + Strong2Weak2, data=dataLog, ties="breslow")
    
    testLog<-coxph(survobjectLog ~ Strong1 + Strong2 + Weak1 + Weak2 + logtStrong1 + logtStrong2, data=dataLog, ties="breslow")
    
    detach(dataLog)
    
    AICH<-(-2*testH$loglik[2])+(2*(length(testH$coefficients)))
    BICH<-(-2*testH$loglik[2])+(log(n)*(length(testH$coefficients)))
    
    AICC<-(-2*testControl$loglik[2])+(2*(length(testControl$coefficients)))
    BICC<-(-2*testControl$loglik[2])+(log(n)*(length(testControl$coefficients)))
    
    AICLog<-(-2*testLog$loglik[2])+(2*(length(testLog$coefficients)))
    BICLog<-(-2*testLog$loglik[2])+(log(n)*(length(testLog$coefficients)))
    
    
    fitTable[i,1]<-AICH
    fitTable[i,2]<-BICH
    fitTable[i,3]<-AICC
    fitTable[i,4]<-BICC
    fitTable[i,5]<-AICLog
    fitTable[i,6]<-BICLog
    
    AICPropTable<-fitTable[,c(1,3,5)]
    BICPropTable<-fitTable[,c(2,4,6)]
    
    
    #AICPropTable$WAICH<-c(0)
    #AICPropTable$WAICH[i]<-ifelse(AICPropTable[i,1]==min(AICPropTable[i,1:3]),1,0)
    
    #AICPropTable$WAICC<-c(0)
    #AICPropTable$WAICC[i]<-ifelse(AICPropTable[i,2]==min(AICPropTable[i,1:3]),1,0)
    
    #AICPropTable$CAICLog<-c(0)
    #AICPropTable$CAICLog[i]<-ifelse(AICPropTable[i,3]==min(AICPropTable[i,1:3]),1,0)
    
    #BICPropTable$WBICH<-c(0)
    #BICPropTable$WBICH[i]<-ifelse(BICPropTable[i,1]==min(BICPropTable[i,1:3]),1,0)
    
    #BICPropTable$WBICC<-c(0)
    #BICPropTable$WBICC[i]<-ifelse(BICPropTable[i,2]==min(BICPropTable[i,1:3]),1,0)
    
    #BICPropTable$CBICLog<-c(0)
    #BICPropTable$CBICLog[i]<-ifelse(BICPropTable[i,3]==min(BICPropTable[i,1:3]),1,0)
    
  }
  
  for(m in 1:length(AICPropTable[,1])){
    AICPropTable$WAICH[m]<-ifelse(AICPropTable[m,1]==min(AICPropTable[m,1:3]),1,0)
    AICPropTable$WAICC[m]<-ifelse(AICPropTable[m,2]==min(AICPropTable[m,1:3]),1,0)
    AICPropTable$CAICLog[m]<-ifelse(AICPropTable[m,3]==min(AICPropTable[m,1:3]),1,0)
    BICPropTable$WBICH[m]<-ifelse(BICPropTable[m,1]==min(BICPropTable[m,1:3]),1,0)
    BICPropTable$WBICC[m]<-ifelse(BICPropTable[m,2]==min(BICPropTable[m,1:3]),1,0)
    BICPropTable$CBICLog[m]<-ifelse(BICPropTable[m,3]==min(BICPropTable[m,1:3]),1,0)
  }
  
  print(AICPropTable)
  print(BICPropTable)
  
  #AICPropTable<-fitTable[,c(1,3,5)]
  #  for(k in 1:reps){ 
  #  AICPropTable$CorrectN<-ifelse(AICPropTable[k,1]==min(AICPropTable[k,1:3]),1,0)
  #  print(AICPropTable)
  #}
  
  
  AICPropWH[l]<-as.numeric(sum(AICPropTable$WAICH)/reps)
  AICPropWC[l]<-sum(AICPropTable$WAICC)/reps
  AICPropCorrect[l]<-sum(AICPropTable$CAICLog)/reps
  
  BICPropWH[l]<-sum(BICPropTable$WBICH)/reps
  BICPropWC[l]<-sum(BICPropTable$WBICC)/reps
  BICPropCorrect[l]<-sum(BICPropTable$CBICLog)/reps
  
  print(c("Proportion of times AIC selected heaviside model across weights."))
  print(AICPropWH)
  
  print(c("Proportion of times AIC selected Control model across weights"))
  print(AICPropWC)
  
  
  print(c("Proportion of times AIC selected (correct) Log model across weights"))
  print(AICPropCorrect)
  
  print(c("Proportion of times BIC selected heaviside model across weights."))
  print(BICPropWH)
  
  print(c("Proportion of times BIC selected Control model across weights"))
  print(BICPropWC)
  
  print(c("Proportion of times BIC selected (correct) Log model across weights"))
  print(BICPropCorrect)  
}

