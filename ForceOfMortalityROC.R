library(survival)
library(ROCR)
library(foreign)

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

#This uses point process model (poorly?)
Ch4Test1Poisson <- glm(Status ~ Clinic + Prison + Dose + offset(log(SurvivalTime)), (family=poisson), data=dsAddicts)
summary(Ch4Test1Poisson)

#dsAddictsOrd <- dsAddicts[order(dsAddicts$SurvivalTime,-dsAddicts$Status),]
#dsAddictsNew <- dsAddictsOrd[!duplicated(dsAddictsOrd$SurvivalTime),]

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

length(Ch4Test1PoissonNew$coefficients)
CompTable2<-cbind(Ch4Test1PoissonNew$coefficients[141:143], Ch4Test1$coefficients)
colnames(CompTable2)<-c("Poisson coefficients", "Cox coefficients")

#So, this table provides the two different estimates of the regression coefficients. So close, comparing them might be pointless.
CompTable2

#Compare anyway!
#So, ROC of poisson model; prediction and performance object
predPoisson2<-prediction(Ch4Test1PoissonNew$linear.predictors, ptProcessDat$Status)
perfPoisson2<-performance(predPoisson2, measure="tpr", x.measure="fpr")

length(Ch4Test1PoissonNew$linear.predictors[141:143])

#Ditto, coxPH model
predCoxph2<-prediction(Ch4Test1$linear.predictors, dsAddicts$Status)
perfCoxph2<-performance(predCoxph2, measure="tpr", x.measure="fpr")

#Graphs!
plot(perfPoisson2, col=c("red"))
par(new=TRUE)
plot(perfCoxph2)
abline(a=0,b=1)

