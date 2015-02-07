#Clearing global environment
rm(list = ls(all = TRUE)) 

library(survival)
library(ROCR)
library(foreign)

install.packages("survivalROC.r")

library(Epi)
#getting practice data
dsVets <- read.dta("http://web1.sph.emory.edu/dkleinb/allDatasets/surv2datasets/vets.dta")
names(dsVets)<-c("tx", "Large", "Adeno", "Small", "Squamous", "survt", "perf", "DisDur", "age", "priortx", "status")

head(dsVets)

################################################
#Some Survival Analysis
#Fitting full COX PH Model
coxphPractice1<-coxph(Surv(dsVets$survt, dsVets$status==1)~tx+Large+Adeno+Small+perf+DisDur+age+priortx, ties="breslow", data=dsVets)
coxphPractice1

#Fitting a Poisson regression model
glm_pr_Practice1<-glm(status~tx+Large+Adeno+Small+perf+DisDur+age+priortx +offset(log(survt)), family="poisson", data=dsVets)
glm_pr_Practice1

#Testing PH Assumption for each predictor
cox.zph(coxphPractice1)

dsVetsEvents<-dsVets[dsVets$status==1,]
coxphPractice1Events<-coxph(Surv(dsVetsEvents$survt, dsVetsEvents$status==1)~tx+Large+Adeno+Small+perf+DisDur+age+priortx, 
                            ties="breslow", data=dsVetsEvents)

cox.zph(coxphPractice1Events, transform="rank")

#Quick Comparison
PoissonPredictorCoefficients<-glm_pr_Practice1$coefficients[2:9]
coefficientComparisonTable<-cbind(as.numeric(coxphPractice1$coefficients), PoissonPredictorCoefficients)
colnames(coefficientComparisonTable)<-c("COX coefficients", "Poisson coefficients")
coefficientComparisonTable

#So, ROC of poisson model; prediction and performance object
predPoisson<-prediction(glm_pr_Practice1$linear.predictors, dsVets$status)
perfPoisson<-performance(predPoisson, measure="tpr", x.measure="fpr")

#Ditto, coxPH model
predCoxph<-prediction(coxphPractice1$linear.predictors, dsVets$status)
perfCoxph<-performance(predCoxph, measure="tpr", x.measure="fpr")

#Graphs!
plot(perfPoisson, col=c("red"))
par(new=TRUE)
plot(perfCoxph)
abline(a=0,b=1)

#But...
fairCompcheck<-cbind(glm_pr_Practice1$linear.predictors, coxphPractice1$linear.predictors)

#Counting process of coxph model
dsAddicts<- read.dta("http://web1.sph.emory.edu/dkleinb/allDatasets/surv2datasets/addicts.dta")

names(dsAddicts)<-c("Subject", "Clinic", "Status", "SurvivalTime", "Prison", "Dose")

#Test Question 1
#So, this gives the regular cox model
Ch4Test1 <- coxph(Surv(SurvivalTime, Status) ~ Clinic + Prison + Dose, ties="breslow", data=dsAddicts)
summary(Ch4Test1)

#This uses point process model
Ch4Test1Poisson <- glm(Status ~ Clinic + Prison + Dose + offset(log(SurvivalTime)), (family=poisson), data=dsAddicts)
summary(Ch4Test1Poisson)

eventTimes <- unique(dsAddicts$SurvivalTime[dsAddicts$Status==1])
eventTimes <- eventTimes[order(eventTimes)]

require(plyr)

createPTable <- function(d){  
  dNew <- d
  for(i in 1:length(eventTimes)){
    dNew[i,] <- d
    dNew[i,"r"] <- i
    dNew[i,"tr"] <- eventTimes[i]
    dNew[i,"dir"] <- ifelse(i==1,min(d$SurvivalTime,eventTimes[i]),min(d$SurvivalTime,eventTimes[i]) - eventTimes[i-1])
    dNew[i,"yir"] <- 0
    if(d$SurvivalTime <= eventTimes[i]) {
      if(d$Status %in% 1) dNew[i,"yir"] <- 1
      break      
    }      
  }    
  return(dNew)
}

ptProcessDat <- ddply(.data=dsAddicts,.variables=.(Subject),.fun = createPTable)
ptProcessDat[ptProcessDat$Subject %in% c(166,111),]
dim(ptProcessDat)
colnames(dsAddicts)
head(ptProcessDat)

#The counting process version of a Cox regression model
Ch4Test1PoissonNew <- glm(yir ~ I(as.factor(r)) + Clinic + Prison + Dose + offset(I(log(dir))), family=poisson(link = "log"), data=ptProcessDat)
summary(Ch4Test1PoissonNew)

#compare to the coxph output
summary(Ch4Test1Poisson)
summary(Ch4Test1)

CompTable2<-cbind(Ch4Test1PoissonNew$coefficients[141:143], Ch4Test1$coefficients)
colnames(CompTable2)<-c("Poisson coefficients", "Cox coefficients")

#So, this table provides the two different estimates of the regression coefficients. So close, comparing them might be pointless.
CompTable2

#Compare anyway!
#So, ROC of poisson model; prediction and performance object
predPoisson2<-prediction(Ch4Test1PoissonNew$linear.predictors, ptProcessDat$yir)
perfPoisson2<-performance(predPoisson2, measure="tpr", x.measure="fpr")

#Ditto, coxPH model
predCoxph2<-prediction(Ch4Test1$linear.predictors, dsAddicts$Status)
perfCoxph2<-performance(predCoxph2, measure="tpr", x.measure="fpr")

#Graphs!
plot(perfPoisson2, col=c("red"))
par(new=TRUE)
plot(perfCoxph2)
abline(a=0,b=1)


#So...
#Calculate area under curve at each value of r?
# Expand ptProcessDat to allow for time-variant predictor variables; look at graphs of AUC across time to compare
# predictive efficiency across time. 

#Code for calculating AUC; a performance object (made from a prediction object) with "measure" set to "auc"
perfPoisson2AUC<-performance(predPoisson2, measure="auc")
perfPoisson2AUC@y.values

perfCoxph2AUC<-performance(predCoxph2, measure="auc")
perfCoxph2AUC@y.values

#First, a regular cox ph model
#Graph the change in AUC across time

aucList<-length(max(ptProcessDat$r))-1
accList<-length(max(ptProcessDat$r))-1
fprList<-length(max(ptProcessDat$r))-1
tprList<-length(max(ptProcessDat$r))-1
fnrList<-length(max(ptProcessDat$r))-1
tnrList<-length(max(ptProcessDat$r))-1


for(i in 2:max(ptProcessDat$r)){
  
  dsSubset<-ptProcessDat[ptProcessDat$r<=i,]
  glmSubset<-glm(yir ~ I(as.factor(r)) + Clinic + Prison + Dose + offset(I(log(dir))), family=poisson(link = "log"), data=dsSubset)
  predSubset<-prediction(glmSubset$linear.predictors, dsSubset$yir)
  perfSubset<-performance(predSubset, measure="auc")
  
  aucList[i-1]<-as.numeric(perfSubset@y.values)
  
  }

#Playing with the code that will go in the for loop
#dsSubset<-ptProcessDat[ptProcessDat$r<=2,]
#glmSubset<-glm(yir ~ I(as.factor(r)) + Clinic + Prison + Dose + offset(I(log(dir))), family=poisson(link = "log"), data=dsSubset)
#predSubset<-prediction(glmSubset$linear.predictors, dsSubset$yir)
#perfSubset<-performance(predSubset, measure="auc")
#perfSubset@y.values
#aucList[1]<-perfSubset@y.values

AucIntervalPairs<-cbind(c(1:length(aucList)),aucList)

plot(AucIntervalPairs, xlab=c("Intervals, r, defined by unique event times"), ylab=c("AUC"))

#Okay, so...schoenfield residuls of clinic show a significant correlation with time
survAddicts<-Surv(dsAddicts$survt, dsAddicts$status==1)
scurve<-survfit(survAddicts~1)
summary(scurve)
plot(scurve, conf.int=TRUE)
coxphTest1<-coxph(Surv(dsAddicts$survt, dsAddicts$status==1)~clinic+prison+dose, ties="breslow", data=dsAddicts)
coxphTest1
cox.zph(coxphTest1, transform="rank")

#Let's have an interaction between Clinic and a heaviside function that goes from 0 to 1 at r=90
#head(ptProcessDat)

HFR90<-as.integer()
length(HFR90)<-length(ptProcessDat$Subject)

for(i in 1:length(ptProcessDat$Subject)){
  HFR90[i]<-ifelse(ptProcessDat[i,8]<428, 0, 1)
}

dsAddictsEx<-cbind(ptProcessDat, HFR90)

#Extended Cox Model
exPoissonHFR90<-glm(yir ~ I(as.factor(r)) + Clinic + Clinic*HFR90 + Prison + Dose + offset(I(log(dir))), family=poisson(link = "log"), data=dsAddictsEx)

#Now make the graph of AUC across time

aucListHFR90<-length(max(dsAddictsEx$r))-1

for(i in 2:max(dsAddictsEx$r)){
  
  dsSubset<-dsAddictsEx[dsAddictsEx$r<=i,]
  glmSubset<-glm(yir ~ I(as.factor(r)) + Clinic + Clinic*HFR90 + Prison + Dose + offset(I(log(dir))), family=poisson(link = "log"), data=dsSubset)
  predSubset<-prediction(glmSubset$linear.predictors, dsSubset$yir)
  perfSubset<-performance(predSubset, measure="auc")

  aucListHFR90[i-1]<-as.numeric(perfSubset@y.values)
  
}


#Graph AUC values across event-time intervals (r), for ph and extended models on same graph
#Increased predictive efficiency after heaviside funciton kicks in
#But, why would we use this, and not just calculate AUC once, on model applied to all data?
AucIntervalPairsHFR90<-cbind(c(1:length(aucListHFR90)),aucListHFR90)
plot(AucIntervalPairsHFR90, xlab=c("Intervals, r, defined by unique event times"), ylab=c("AUC"), col=c("red"))
par(new=TRUE)
plot(AucIntervalPairs, xlab=c("Intervals, r, defined by unique event times"), ylab=c("AUC"))

#Let's try this: HFR categories at r = 40, 80, 120 (event time at tr=168, 358, 652  )

#peek<-ptProcessDat[ptProcessDat$r==120,]
#head(peek)

HFR4080120<-as.integer()
length(HFR4080120)<-length(ptProcessDat$Subject)


for(i in 1:length(ptProcessDat$Subject)){
  HFR4080120[i]<-ifelse(ptProcessDat[i,8]<168, 0, ifelse(ptProcessDat[i,8]<358, 1, 2))
}

dsAddictsEx<-cbind(ptProcessDat, HFR90, HFR4080120)

dsAddictsEx$HFR4080120[16000:16100]
factor(dsAddicts$HFR4080120)

aucListHFR4080120<-length(max(dsAddictsEx$r))-1

for(i in 2:max(dsAddictsEx$r)){
  
  dsSubset<-dsAddictsEx[dsAddictsEx$r<=i,]
  glmSubset<-glm(yir ~ I(as.factor(r)) + Clinic + Clinic*HFR4080120 + Prison + Dose + offset(I(log(dir))), family=poisson(link = "log"), data=dsSubset)
  predSubset<-prediction(glmSubset$linear.predictors, dsSubset$yir)
  perfSubset<-performance(predSubset, measure="auc")
  
  aucListHFR4080120[i-1]<-as.numeric(perfSubset@y.values)
}

AucIntervalPairsHFR4080120<-cbind(c(1:length(aucListHFR4080120)),aucListHFR4080120)


#Results are 100% to be expected; no surprises, so why use this method instead of just comparing model fit? 
#Or just AUC of model over all data?
plot(AucIntervalPairsHFR4080120, xlab=c("Intervals, r, defined by unique event times"), ylab=c("AUC"), col=c("blue"))
par(new=TRUE)
plot(AucIntervalPairsHFR90, xlab=c("Intervals, r, defined by unique event times"), ylab=c("AUC"), col=c("red"))
par(new=TRUE)
plot(AucIntervalPairs, xlab=c("Intervals, r, defined by unique event times"), ylab=c("AUC"))


#...because...let's look at some measures other than AUC!

aucList<-length(max(ptProcessDat$r))-1
accList<-length(max(ptProcessDat$r))-1
fprList<-length(max(ptProcessDat$r))-1
tprList<-length(max(ptProcessDat$r))-1
fnrList<-length(max(ptProcessDat$r))-1
tnrList<-length(max(ptProcessDat$r))-1


for(i in 2:max(ptProcessDat$r)){
  
  dsSubset<-ptProcessDat[ptProcessDat$r<=i,]
  glmSubset<-glm(yir ~ I(as.factor(r)) + Clinic + Prison + Dose + offset(I(log(dir))), family=poisson(link = "log"), data=dsSubset)
  predSubset<-prediction(glmSubset$linear.predictors, dsSubset$yir)
  perfSubsetAUC<-performance(predSubset, measure="auc")
  perfSubsetACC<-performance(predSubset, "acc")
  perfSubsetFPR<-performance(predSubset, measure="fpr")
  perfSubsetTPR<-performance(predSubset, measure="tpr")
  perfSubsetFNR<-performance(predSubset, measure="fnr")
  perfSubsetTNR<-performance(predSubset, measure="tnr")
  
  
  aucList[i-1]<-as.numeric(perfSubsetAUC@y.values)
  accList[i-1]<-max(as.numeric(unlist(perfSubsetACC@y.values)))
  fprList[i-1]<-max(as.numeric(unlist(perfSubsetFPR@y.values)))
  tprList[i-1]<-max(as.numeric(unlist(perfSubsetTPR@y.values)))
  fnrList[i-1]<-max(as.numeric(unlist(perfSubsetFNR@y.values)))
  tnrList[i-1]<-max(as.numeric(unlist(perfSubsetTNR@y.values)))
  
}

AucIntervalPair<-cbind(c(1:length(aucList)),aucList)
plot(AucIntervalPair)

AccIntervalPair<-cbind(c(1:length(accList)),accList)
plot(AccIntervalPair)

FprIntervalPair<-cbind(c(1:length(fprList)),fprList)
plot(FprIntervalPair)

TprIntervalPair<-cbind(c(1:length(tprList)),tprList)
plot(TprIntervalPair)

FnrIntervalPair<-cbind(c(1:length(fnrList)),fnrList)
plot(FnrIntervalPair)

TnrIntervalPair<-cbind(c(1:length(tnrList)),tnrList)
plot(TnrIntervalPair)


#Now for first extended cox model
#HFR90

HF90aucList<-length(max(ptProcessDat$r))-1
HF90accList<-length(max(ptProcessDat$r))-1
HF90fprList<-length(max(ptProcessDat$r))-1
HF90tprList<-length(max(ptProcessDat$r))-1
HF90fnrList<-length(max(ptProcessDat$r))-1
HF90tnrList<-length(max(ptProcessDat$r))-1


for(i in 2:max(dsAddictsEx$r)){
  
  dsSubset<-dsAddictsEx[dsAddictsEx$r<=i,]
  glmSubset<-glm(yir ~ I(as.factor(r)) + Clinic + Clinic*HFR90 + Prison + Dose + offset(I(log(dir))), family=poisson(link = "log"), data=dsSubset)
  predSubset<-prediction(glmSubset$linear.predictors, dsSubset$yir)
  perfSubsetHF90AUC<-performance(predSubset, measure="auc")
  perfSubsetHF90ACC<-performance(predSubset, "acc")
  perfSubsetHF90FPR<-performance(predSubset, measure="fpr")
  perfSubsetHF90TPR<-performance(predSubset, measure="tpr")
  perfSubsetHF90FNR<-performance(predSubset, measure="fnr")
  perfSubsetHF90TNR<-performance(predSubset, measure="tnr")
  
  
  HF90aucList[i-1]<-as.numeric(perfSubsetHF90AUC@y.values)
  HF90accList[i-1]<-max(as.numeric(unlist(perfSubsetHF90ACC@y.values)))
  HF90fprList[i-1]<-max(as.numeric(unlist(perfSubsetHF90FPR@y.values)))
  HF90tprList[i-1]<-max(as.numeric(unlist(perfSubsetHF90TPR@y.values)))
  HF90fnrList[i-1]<-max(as.numeric(unlist(perfSubsetHF90FNR@y.values)))
  HF90tnrList[i-1]<-max(as.numeric(unlist(perfSubsetHF90TNR@y.values)))
  
}



HF90AucIntervalPair<-cbind(c(1:length(aucList)),aucList)
plot(HF90AucIntervalPair, col=c("red"))
par(new=TRUE)
plot(AucIntervalPairs, xlab=c("Intervals, r, defined by unique event times"), ylab=c("AUC"))


HF90AccIntervalPair<-cbind(c(1:length(accList)),accList)
plot(HF90AccIntervalPair)
par(new=TRUE)
plot(AucIntervalPairs, xlab=c("Intervals, r, defined by unique event times"), ylab=c("AUC"))


HF90FprIntervalPair<-cbind(c(1:length(fprList)),fprList)
plot(HF90FprIntervalPair)

HF90TprIntervalPair<-cbind(c(1:length(tprList)),tprList)
plot(HF90TprIntervalPair)

HF90FnrIntervalPair<-cbind(c(1:length(fnrList)),fnrList)
plot(HF90FnrIntervalPair)

HF90TnrIntervalPair<-cbind(c(1:length(tnrList)),tnrList)
plot(HF90TnrIntervalPair)
