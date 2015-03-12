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
AICPropWI<-vector(length=length(ftw))

BICPropCorrect<-vector(length=length(ftw))
BICPropWC<-vector(length=length(ftw))
BICPropWH<-vector(length=length(ftw))
BICPropWI<-vector(length=length(ftw))

#For: weight of a specific time function (heaviside, the last one)
for(l in 1:length(ftw)){
  #Declaring Betas: first, we are looking at the heaviside function, so the last two time-dependent 
  #variables are given betas, all else are given 0 weight. 
  betas<-c(0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0, 0, 0, 0, ftw[l], ftw[l], 0, 0)
  
  #Creating a table of AICs and BIC values
  reps<-10
  
  fitTable<-data.frame(matrix(ncol=8, nrow=reps, ))
  colnames(fitTable)<-c("AICC", "BICC", "AICH", "BICH", "AICI","BICI", "AICLog","BICLog")
  
  #For this weight of specific time function, create this many sets of data  
  for(i in 1:reps){
    
    n=500
    m=365
    
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
    
    #10th and 11th columns are going to represent time squared * Strong1 and Strong2, respectively
    xmat[,10]<-xmat[,7]*xmat[,7]*xmat[,1]
    xmat[,11]<-xmat[,7]*xmat[,7]*xmat[,2]
    
    #12th and 13th columns are going to represent interactions between Strong1*Weak1 and Strong2*Weak2
    xmat[,12]<-xmat[,1]*xmat[,4]
    xmat[,13]<-xmat[,2]*xmat[,5]
    
    #14th and 15th columns are going to represent interactions between heaviside function and strong1/strong2 respectively
    for(j in 1:length(xmat[,7])){
      xmat[j,14]<-ifelse(xmat[j,7]<=182, 0, xmat[j,1])  
      xmat[j,15]<-ifelse(xmat[j,7]<=182, 0, xmat[j,2])
    }
    
    dsMaster<-as.data.frame(xmat)
    head(dsMaster)
    colnames(dsMaster)<-c("Strong1", "Strong2", "Strong3", "Weak1", "Weak2", "Weak3", "intForFt", "logtStrong1", "logtStrong2", "t2Strong1", "t2Strong2", "Strong1Weak1", "Strong2Weak2", "Strong1H", "Strong2H")
    
    #Strong and Weak variables, plust two heaviside function variables
    dsLog<-as.matrix(dsMaster[,c(1:6, 8:15)])
    eventTimesMaybe<-runif(n, 1, m)
    
    dataLog<-permalgorithm(n, m, Xmat=dsLog, XmatNames=c("Strong1", "Strong2", "Strong3", "Weak1", "Weak2", "Weak3", "logtStrong1", "logtStrong2", "t2Strong1", "t2Strong2", "Strong1Weak1", "Strong2Weak2", "Strong1H", "Strong2H"), eventRandom=eventTimesMaybe, betas=betas)

    attach(dataLog)
    survobjectLog<-Surv(time=Start, time2=Stop, Event==1)
    
    testControl<-coxph(survobjectLog ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3, data=dataLog, ties="breslow")
    
    testH<-coxph(survobjectLog ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3 + Strong1H + Strong2H, data=dataLog, ties="breslow")
    
    testInt<-coxph(survobjectLog ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3 + Strong1Weak1 + Strong2Weak2, data=dataLog, ties="breslow")
    
    testLog<-coxph(survobjectLog ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3 + logtStrong1 + logtStrong2, data=dataLog, ties="breslow")
    
    detach(dataLog)
    
    AICC<-(-2*testControl$loglik[2])+(2*(length(testControl$coefficients)))
    BICC<-(-2*testControl$loglik[2])+(log(n)*(length(testControl$coefficients)))
    
    AICH<-(-2*testH$loglik[2])+(2*(length(testH$coefficients)))
    BICH<-(-2*testH$loglik[2])+(log(n)*(length(testH$coefficients)))
    
    AICI<-(-2*testControl$loglik[2])+(2*(length(testInt$coefficients)))
    BICI<-(-2*testControl$loglik[2])+(log(n)*(length(testInt$coefficients)))
    
    AICLog<-(-2*testLog$loglik[2])+(2*(length(testLog$coefficients)))
    BICLog<-(-2*testLog$loglik[2])+(log(n)*(length(testLog$coefficients)))
    
    fitTable[i,1]<-AICC
    fitTable[i,2]<-BICC
    fitTable[i,3]<-AICH
    fitTable[i,4]<-BICH
    fitTable[i,5]<-AICI
    fitTable[i,6]<-BICI
    fitTable[i,7]<-AICLog
    fitTable[i,8]<-BICLog
    
    AICPropTable<-fitTable[,c(1,3,5,7)]
    BICPropTable<-fitTable[,c(2,4,6,8)]
  }
  
  for(m in 1:length(AICPropTable[,1])){
    
    AICPropTable$WAICC[m]<-ifelse(AICPropTable[m,1]==min(AICPropTable[m,1:4]),1,0)
    AICPropTable$WAICH[m]<-ifelse(AICPropTable[m,2]==min(AICPropTable[m,1:4]),1,0)
    AICPropTable$CAICI[m]<-ifelse(AICPropTable[m,3]==min(AICPropTable[m,1:4]),1,0)
    AICPropTable$WAICLog[m]<-ifelse(AICPropTable[m,4]==min(AICPropTable[m,1:4]),1,0)
    BICPropTable$WBICC[m]<-ifelse(BICPropTable[m,1]==min(BICPropTable[m,1:4]),1,0)
    BICPropTable$WBICH[m]<-ifelse(BICPropTable[m,2]==min(BICPropTable[m,1:4]),1,0)
    BICPropTable$CBICI[m]<-ifelse(BICPropTable[m,3]==min(BICPropTable[m,1:4]),1,0)
    BICPropTable$WBICLog[m]<-ifelse(BICPropTable[m,4]==min(BICPropTable[m,1:4]),1,0)
  }
  
  print(AICPropTable)
  print(BICPropTable)
  
  AICPropWC[l]<-as.numeric(sum(AICPropTable$WAICC)/reps)
  AICPropWH[l]<-as.numeric(sum(AICPropTable$WAICH)/reps)
  AICPropCI[l]<-sum(AICPropTable$CAICI)/reps
  AICPropWL[l]<-sum(AICPropTable$WAICLog)/reps
  
  BICPropWC[l]<-sum(BICPropTable$WBICC)/reps
  BICPropWH[l]<-sum(BICPropTable$WBICH)/reps
  BICPropCorrect[l]<-sum(BICPropTable$CBICI)/reps
  BICPropWL[l]<-sum(BICPropTable$WBICLog)/reps
  
  print(c("Proportion of times AIC selected no-time model across weights."))
  print(AICPropWC)
  
  print(c("Proportion of times AIC selected heaviside model across weights."))
  print(AICPropWH)
  
  print(c("Proportion of times AIC selected (correct) Interaction model across weights"))
  print(AICPropCorrect)
  
  print(c("Proportion of times AIC selected Log model across weights"))
  print(AICPropWL)
  
  print(c("Proportion of times BIC selected no-time model across weights."))
  print(BICPropWC)
  
  print(c("Proportion of times BIC selected heaviside model across weights."))
  print(BICPropWH)
  
  print(c("Proportion of times BIC selected (correct) Interaction model across weights"))
  print(BICPropCorrect)
  
  print(c("Proportion of times BIC selected Log model across weights"))
  print(BICPropWL)  
}

GraphVectorAIC<-cbind(ftw, AICPropCorrect)
plot(GraphVectorAIC)

GraphVectorBIC<-cbind(ftw, BICPropCorrect)
plot(GraphVectorBIC)
