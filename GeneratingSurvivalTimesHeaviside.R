rm(list = ls(all.names = TRUE))
library(cvTools)
library(survsim)
library(survival)
library(foreign)
library(PermAlgo)
library(ROCR)

#ftw<-c(0.7)
#ftw<-c(0, 0.01, 0.35, 0.7)
ftw<-c(1:100)/100

cvPropWL<-vector(length=length(ftw))
cvPropWC<-vector(length=length(ftw))
cvPropCH<-vector(length=length(ftw))
cvPropWI<-vector(length=length(ftw))
cvPropWT<-vector(length=length(ftw))

AUCPropWL<-vector(length=length(ftw))
AUCPropWC<-vector(length=length(ftw))
AUCPropCH<-vector(length=length(ftw))
AUCPropWI<-vector(length=length(ftw))
AUCPropWT<-vector(length=length(ftw))

AICPropWL<-vector(length=length(ftw))
AICPropWC<-vector(length=length(ftw))
AICPropCH<-vector(length=length(ftw))
AICPropWI<-vector(length=length(ftw))
AICPropWT<-vector(length=length(ftw))

BICPropWL<-vector(length=length(ftw))
BICPropWC<-vector(length=length(ftw))
BICPropCH<-vector(length=length(ftw))
BICPropWI<-vector(length=length(ftw))
BICPropWT<-vector(length=length(ftw))

#For: weight of a specific time function (heaviside, the last one)
for(l in 1:length(ftw)){
  #Declaring Betas: first, we are looking at the heaviside function, so the last two time-dependent 
  #variables are given betas, all else are given 0 weight. 
  betas<-c(0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0, 0, 0, 0, 0, 0, ftw[l], ftw[l])
  
  #Creating a table of AICs and BIC values
  reps<-30
  
  fitTable<-data.frame(matrix(ncol=20, nrow=reps, ))
  colnames(fitTable)<-c("AICC",  "AICH", "AICI", "AICLog", "AICT", "BICC","BICH", "BICI","BICLog", "BICT" , "cvC", "cvH", "cvI", "cvLog", "cvT", "AUCC", "AUCH", "AUCI", "AUCLog", "AUCT")
  
  #For this weight of specific time function, create this many sets of data  
  for(i in 1:reps){
    
    n=1000
    m=24
    mhalf=m/2
    
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
    xmat[,10]<-(xmat[,7])*xmat[,1]
    xmat[,11]<-(xmat[,7])*xmat[,2]
    
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
    
    #Strong and Weak variables, plust two heaviside function variables
    dsH<-as.matrix(dsMaster[,c(1:6, 8:15)])
    eventTimesMaybe<-runif(n, 1, m)
    
    dataH<-permalgorithm(n, m, Xmat=dsH, XmatNames=c("Strong1", "Strong2", "Strong3", "Weak1", "Weak2", "Weak3", "logtStrong1", "logtStrong2", "tStrong1", "tStrong2", "Strong1Weak1", "Strong2Weak2", "Strong1H", "Strong2H"), eventRandom=eventTimesMaybe, betas=betas)
    
    attach(dataH)
    survobjectH<-Surv(time=Start, time2=Stop, Event==1)
    
    #Creating the cox models for each set of variables
    
    testControl<-coxph(survobjectH ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3, data=dataH, ties="breslow")
    
    testH<-coxph(survobjectH ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3 + Strong1H + Strong2H, data=dataH, ties="breslow")
    
    testInt<-coxph(survobjectH ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3 + Strong1Weak1 + Strong2Weak2, data=dataH, ties="breslow")
    
    testLog<-coxph(survobjectH ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3 + logtStrong1 + logtStrong2, data=dataH, ties="breslow")
    
    testT<-coxph(survobjectH ~ Strong1 + Strong2 + Strong3 + Weak1 + Weak2 + Weak3 + tStrong1 + tStrong2, data=dataH, ties="breslow")
    
    #Getting cv error for each model
    cvF<-cvFolds(n=length(dataH$Event), K=10)
    cvTestH<-cvTool(call=testH, y=dataH$Event, folds=cvF)
    cvTestI<-cvTool(call=testInt, y=dataH$Event, folds=cvF)
    cvTestLog<-cvTool(call=testLog, y=dataH$Event, folds=cvF)
    cvTestC<-cvTool(call=testControl, y=dataH$Event, folds=cvF)
    cvTestT<-cvTool(call=testT, y=dataH$Event, folds=cvF)
    
    predTestControl<-prediction(testControl$linear.predictors, dataH$Event)
    perfTestControl<-performance(predTestControl, measure="auc")
    AUCC<-as.numeric(perfTestControl@y.values)
    
    predTestH<-prediction(testH$linear.predictors, dataH$Event)
    perfTestH<-performance(predTestH, measure="auc")
    AUCH<-as.numeric(perfTestH@y.values)
    
    predTestI<-prediction(testInt$linear.predictors, dataH$Event)
    perfTestI<-performance(predTestI, measure="auc")
    AUCI<-as.numeric(perfTestI@y.values)
    
    predTestLog<-prediction(testLog$linear.predictors, dataH$Event)
    perfTestLog<-performance(predTestLog, measure="auc")
    AUCLog<-as.numeric(perfTestLog@y.values)
    
    predTestT<-prediction(testT$linear.predictors, dataH$Event)
    perfTestT<-performance(predTestT, measure="auc")
    AUCT<-as.numeric(perfTestT@y.values)
    
    detach(dataH)
    
    #Getting AIC and BIC for each model
    AICC<-(-2*testControl$loglik[2])+(2*(length(testControl$coefficients)))
    BICC<-(-2*testControl$loglik[2])+(log(n)*(length(testControl$coefficients)))
    
    AICH<-(-2*testH$loglik[2])+(2*(length(testH$coefficients)))
    BICH<-(-2*testH$loglik[2])+(log(n)*(length(testH$coefficients)))
    
    AICI<-(-2*testInt$loglik[2])+(2*(length(testInt$coefficients)))
    BICI<-(-2*testInt$loglik[2])+(log(n)*(length(testInt$coefficients)))
    
    AICLog<-(-2*testLog$loglik[2])+(2*(length(testLog$coefficients)))
    BICLog<-(-2*testLog$loglik[2])+(log(n)*(length(testLog$coefficients)))
    
    AICT<-(-2*testT$loglik[2])+(2*(length(testT$coefficients)))
    BICT<-(-2*testT$loglik[2])+(log(n)*(length(testT$coefficients)))
    
    fitTable[i,1]<-AICC
    fitTable[i,2]<-AICH
    fitTable[i,3]<-AICI
    fitTable[i,4]<-AICLog
    fitTable[i,5]<-AICT
    fitTable[i,6]<-BICC
    fitTable[i,7]<-BICH
    fitTable[i,8]<-BICI
    fitTable[i,9]<-BICLog
    fitTable[i,10]<-BICT
    fitTable[i,11]<-cvTestC
    fitTable[i,12]<-cvTestH
    fitTable[i,13]<-cvTestI
    fitTable[i,14]<-cvTestLog
    fitTable[i,15]<-cvTestT
    fitTable[i,16]<-AUCC
    fitTable[i,17]<-AUCH
    fitTable[i,18]<-AUCI
    fitTable[i,19]<-AUCLog
    fitTable[i,20]<-AUCT
    
    AICPropTable<-fitTable[,c(1:5)]
    BICPropTable<-fitTable[,c(6:10)]
    cvPropTable<-fitTable[,c(11:15)]
    AUCPropTable<-fitTable[,c(16:20)]
  }
  
  for(m in 1:length(AICPropTable[,1])){
    
    AICPropTable$WAICC[m]<-ifelse(AICPropTable[m,1]==min(AICPropTable[m,1:5]),1,0)
    AICPropTable$CAICH[m]<-ifelse(AICPropTable[m,2]==min(AICPropTable[m,1:5]),1,0)
    AICPropTable$WAICI[m]<-ifelse(AICPropTable[m,3]==min(AICPropTable[m,1:5]),1,0)
    AICPropTable$WAICLog[m]<-ifelse(AICPropTable[m,4]==min(AICPropTable[m,1:5]),1,0)
    AICPropTable$WAICT[m]<-ifelse(AICPropTable[m,5]==min(AICPropTable[m,1:5]),1,0)
    
    BICPropTable$WBICC[m]<-ifelse(BICPropTable[m,1]==min(BICPropTable[m,1:5]),1,0)
    BICPropTable$CBICH[m]<-ifelse(BICPropTable[m,2]==min(BICPropTable[m,1:5]),1,0)
    BICPropTable$WBICI[m]<-ifelse(BICPropTable[m,3]==min(BICPropTable[m,1:5]),1,0)
    BICPropTable$WBICLog[m]<-ifelse(BICPropTable[m,4]==min(BICPropTable[m,1:5]),1,0)
    BICPropTable$WBICT[m]<-ifelse(BICPropTable[m,5]==min(BICPropTable[m,1:5]),1,0)
    
    cvPropTable$WcvC[m]<-ifelse(cvPropTable[m,1]==min(cvPropTable[m,1:5]),1,0)
    cvPropTable$CcvH[m]<-ifelse(cvPropTable[m,2]==min(cvPropTable[m,1:5]),1,0)
    cvPropTable$WcvI[m]<-ifelse(cvPropTable[m,3]==min(cvPropTable[m,1:5]),1,0)
    cvPropTable$WcvLog[m]<-ifelse(cvPropTable[m,4]==min(cvPropTable[m,1:5]),1,0)
    cvPropTable$WcvT[m]<-ifelse(cvPropTable[m,5]==min(cvPropTable[m,1:5]),1,0)
    
    AUCPropTable$WAUCC[m]<-ifelse(AUCPropTable[m,1]==max(AUCPropTable[m,1:5]),1,0)
    AUCPropTable$CAUCH[m]<-ifelse(AUCPropTable[m,2]==max(AUCPropTable[m,1:5]),1,0)
    AUCPropTable$WAUCI[m]<-ifelse(AUCPropTable[m,3]==max(AUCPropTable[m,1:5]),1,0)
    AUCPropTable$WAUCLog[m]<-ifelse(AUCPropTable[m,4]==max(AUCPropTable[m,1:5]),1,0)
    AUCPropTable$WAUCT[m]<-ifelse(AUCPropTable[m,5]==max(AUCPropTable[m,1:5]),1,0)
    
  }
  
  print(AICPropTable)
  print(BICPropTable)
  print(cvPropTable)
  print(AUCPropTable)
  
  AICPropWC[l]<-sum(AICPropTable$WAICC)/reps
  AICPropCH[l]<-sum(AICPropTable$CAICH)/reps
  AICPropWI[l]<-sum(AICPropTable$WAICI)/reps
  AICPropWL[l]<-sum(AICPropTable$WAICLog)/reps
  AICPropWT[l]<-sum(AICPropTable$WAICT)/reps
  
  BICPropWC[l]<-sum(BICPropTable$WBICC)/reps
  BICPropCH[l]<-sum(BICPropTable$CBICH)/reps
  BICPropWI[l]<-sum(BICPropTable$WBICI)/reps
  BICPropWL[l]<-sum(BICPropTable$WBICLog)/reps
  BICPropWT[l]<-sum(BICPropTable$WBICT)/reps
  
  cvPropWC[l]<-sum(cvPropTable$WcvC)/reps
  cvPropCH[l]<-sum(cvPropTable$CcvH)/reps
  cvPropWI[l]<-sum(cvPropTable$WcvI)/reps
  cvPropWL[l]<-sum(cvPropTable$WcvLog)/reps
  cvPropWT[l]<-sum(cvPropTable$WcvT)/reps
  
  AUCPropWC[l]<-sum(AUCPropTable$WAUCC)/reps
  AUCPropCH[l]<-sum(AUCPropTable$CAUCH)/reps
  AUCPropWI[l]<-sum(AUCPropTable$WAUCI)/reps
  AUCPropWL[l]<-sum(AUCPropTable$WAUCLog)/reps  
  AUCPropWT[l]<-sum(AUCPropTable$WAUCT)/reps
  
  print(testH)
  
  print(c("Proportion of times AIC selected no-time model across weights."))
  print(AICPropWC)
  
  print(c("Proportion of times AIC selected (correct) heaviside model across weights."))
  print(AICPropCH)
  
  print(c("Proportion of times AIC selected Interaction model across weights"))
  print(AICPropWI)
  
  print(c("Proportion of times AIC selected Log model across weights"))
  print(AICPropWL)
  
  print(c("Proportion of times AIC selected time model across weights"))
  print(AICPropWT)
  
  print(c("Proportion of times BIC selected no-time model across weights."))
  print(BICPropWC)
  
  print(c("Proportion of times BIC selected (correct) heaviside model across weights."))
  print(BICPropCH)
  
  print(c("Proportion of times BIC selected Interaction model across weights"))
  print(BICPropWI)
  
  print(c("Proportion of times BIC selected Log model across weights"))
  print(BICPropWL)  
  
  print(c("Proportion of times BIC selected time model across weights"))
  print(BICPropWT)  
  
  print(c("Proportion of times cv selected no-time model across weights."))
  print(cvPropWC)
  
  print(c("Proportion of times cv selected (correct) heaviside model across weights."))
  print(cvPropCH)
  
  print(c("Proportion of times cv selected Interaction model across weights"))
  print(cvPropWI)
  
  print(c("Proportion of times cv selected Log model across weights"))
  print(cvPropWL)  
  
  print(c("Proportion of times cv selected time model across weights"))
  print(cvPropWT)
  
  print(c("Proportion of times AUC selected no-time model across weights."))
  print(AUCPropWC)
  
  print(c("Proportion of times AUC selected (correct) heaviside model across weights."))
  print(AUCPropCH)
  
  print(c("Proportion of times AUC selected Interaction model across weights"))
  print(AUCPropWI)
  
  print(c("Proportion of times AUC selected Log model across weights"))
  print(AUCPropWL)  
  
  print(c("Proportion of times AUC selected time model across weights"))
  print(AUCPropWT)  
}

GraphVectorAIC<-cbind(ftw, AICPropCI)
plot(GraphVectorAIC)

GraphVectorBIC<-cbind(ftw, BICPropCI)
plot(GraphVectorBIC)

GraphVectorcv<-cbind(ftw, cvPropCorrect)
plot(GraphVectorcv)

GraphVectorAUC<-cbind(ftw, AUCPropCorrect)
plot(GraphVectorAUC)

#Run these when sim is completed. 
#SimInt500N320<-rbind(ftw, AICPropWC, AICPropWH, AICPropCI, AICPropWL, BICPropWC, BICPropWH, BICPropCI, BICPropWL, cvPropWC, cvPropWH, cvPropCI, cvPropWL, AUCPropWC, AUCPropWH, AUCPropCI, AUCPropWL)
#write.csv(SimInt500N320, file="SimInt500N320.csv", na=".")
