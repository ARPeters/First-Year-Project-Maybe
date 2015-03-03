rm(list = ls(all.names = TRUE))
library(survsim)
library(survival)
library(foreign)
library(PermAlgo)


ftw<-c(0.01,0.35, 0.7)

AICPropCorrect<-vector(length=length(ftw))
BICPropCorrect<-vector(length=length(ftw))

#For: weight of a specific time function (heaviside, the last one)
for(l in 1:length(ftw)){
  #Declaring Betas: first, we are looking at the heaviside function, so the last two time-dependent 
  #variables are given betas, all else are given 0 weight. 
  betas<-c(0.7, 0.7, 0.1, 0.1, 0, 0, 0, 0, 0, 0, ftw[l], ftw[l])

  #Creating a table of AICs and BIC values
  reps<-2
  
  fitTable<-data.frame(matrix(ncol=6, nrow=reps, ))
  colnames(fitTable)<-c("AICH", "BICH", "AICC","BICC", "AICLog","BICLog")
  
  #For this weight of specific time function, create this many sets of data  
  for(i in 1:reps){
    
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
    for(j in 1:length(xmat[,5])){
      xmat[j,12]<-ifelse(xmat[j,5]<=182, 0, xmat[j,1])  
      xmat[j,13]<-ifelse(xmat[j,5]<=182, 0, xmat[j,2])
    }
    
    dsMaster<-as.data.frame(xmat)
    
    colnames(dsMaster)<-c("Strong1", "Strong2", "Weak1", "Weak2", "intForFt", "logtStrong1", "logtStrong2", "t2Strong1", "t2Strong2", "Strong1Weak1", "Strong2Weak2", "Strong1H", "Strong2H")
    
    #Strong and Weak variables, plust two heaviside function variables
    dsH<-as.matrix(dsMaster[,c(1:4,6:13)])
    eventTimesMaybe<-runif(n, 1, m)
    
    dataH<-permalgorithm(n, m, Xmat=dsH, XmatNames=c("Strong1", "Strong2", "Weak1", "Weak2", "logtStrong1", "logtStrong2", "t2Strong1", "t2Strong2", "Strong1Weak1", "Strong2Weak2", "Strong1H", "Strong2H"), eventRandom=eventTimesMaybe, betas=betas)
    
    attach(dataH)
    survobjectH<-Surv(time=Start, time2=Stop, Event==1)
    
    testH<-coxph(survobjectH ~ Strong1 + Strong2 + Weak1 + Weak2 + Strong1H + Strong2H, data=dataH, ties="breslow")
    
    testControl<-coxph(survobjectH ~ Strong1 + Strong2 + Weak1 + Weak2 + Strong1Weak1 + Strong2Weak2, data=dataH)
    
    testLog<-coxph(survobjectH ~ Strong1 + Strong2 + Weak1 + Weak2 + logtStrong1 + logtStrong2, data=dataH)
    
    detach(dataH)
    
    AICH<-(-2*testH$loglik[2])+2*(length(testH$coefficients))
    BICH<-(-2*testH$loglik[2])+log(n)*(length(testH$coefficients))
    
    AICC<-(-2*testControl$loglik[2])+2*(length(testControl$coefficients))
    BICC<-(-2*testControl$loglik[2])+log(n)*(length(testControl$coefficients))
    
    AICLog<-(-2*testLog$loglik[2])+2*(length(testLog$coefficients))
    BICLog<-(-2*testLog$loglik[2])+log(n)*(length(testH$coefficients))
    
    
    fitTable[i,1]<-AICH
    fitTable[i,2]<-BICH
    fitTable[i,3]<-AICC
    fitTable[i,4]<-BICC
    fitTable[i,5]<-AICLog
    fitTable[i,6]<-BICLog
    
    AICPropTable<-fitTable[,c(1,3,5)]
    BICPropTable<-fitTable[,c(2,4,6)]
    
    
    AICPropTable$Correct<-ifelse(AICPropTable[i,1]==min(AICPropTable[i,1:3]),1,0)
    BICPropTable$Correct<-ifelse(BICPropTable[i,1]==min(BICPropTable[i,1:3]),1,0)
  }
  
  #AICPropTable<-fitTable[,c(1,3,5)]
  #  for(k in 1:reps){ 
  #  AICPropTable$CorrectN<-ifelse(AICPropTable[k,1]==min(AICPropTable[k,1:3]),1,0)
  #  print(AICPropTable)
  #}
  
  AICPropCorrect[l]<-sum(AICPropTable$Correct)/reps
  BICPropCorrect[l]<-sum(BICPropTable$Correct)/reps
  print(AICPropCorrect)
  print(BICPropCorrect)
}

print(AICPropCorrect)
print(BICPropCorrect)
