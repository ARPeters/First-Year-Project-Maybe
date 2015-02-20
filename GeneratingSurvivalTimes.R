rm(list = ls(all.names = TRUE))
library(survsim)
library(survival)
library(foreign)
library(PermAlgo)
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

str(weibullTest)
weibullTest$loglik[2]


#Functional, but apparently cannot accomodate interactions, or time-dependent covariates. So, different tactic:
#Bender et al, 2005
#Control Data
#X5 and X6 are interactions between 1/3 and 2/4
X1<-rbinom(5000, 1, 0.5)
X2<-rbinom(5000, 1, 0.5)
X3<-rbinom(5000, 1, 0.5)
X4<-rbinom(5000, 1, 0.5)
X5<-X1*X3
X6<-X2*X4
U<-runif(5000)
scale<-1
negUlog<--1*log(U)
betaVector<-c(0.7, 0.7, 0.1, 0.1, 0.7, 0.7)
hf<-scale*(exp(X1*betaVector[1]+X2*betaVector[2]+X3*betaVector[3]+X4*betaVector[4]+X5*betaVector[5]+X6*betaVector[6]))
T<-negUlog/hf
T<-T*100
status<-1
dsNew<-cbind(X1, X2, X3, X4, X5, X6, scale, U, negUlog, hf, T, status)
dsNew<-as.data.frame(dsNew)


survobject<-Surv(time=dsNew$T, event=dsNew$status)
coxphTest<-coxph(survobject~X1 + X2 + X3 + X4 + X5 + X6, data=dsNew)

weibullTest<-survreg(survobject~X1 + X2 + X3 + X4 + X5 + X6, data=dsNew, dist="weibull")


## Example from PermoAlgo package
# Prepare the matrice of covariate (Xmat)
# Here we simulate daily exposures to 2 prescription drugs over a
# year. Drug prescriptions can start any day of follow-up, and their
# duration is a multiple of 7 days. There can be multiple prescriptions
# for each individuals over the year and interuptions of drug use in
# between.
# Additionaly, there is a time-independant binary covarite (sex).
n=500 # subjects
m=365 # days
# Generate the matrix of three covariate, in a 'long' format.
Xmat=matrix(ncol=3, nrow=n*m)
# time-independant binary covariate
Xmat[,1] <- rep(rbinom(n, 1, 0.3), each=m)
# Function to generate an individual time-dependent exposure history
# e.g. generate prescriptions of different durations and doses.
TDhist <- function(m){
  start <- round(runif(1,1,m),0) # individual start date
  duration <- 7 + 7*rpois(1,3) # in weeks
  dose <- round(runif(1,0,10),1)
  vec <- c(rep(0, start-1), rep(dose, duration))
  while (length(vec)<=m){
    intermission <- 21 + 7*rpois(1,3) # in weeks
    duration <- 7 + 7*rpois(1,3) # in weeks
    dose <- round(runif(1,0,10),1)
    vec <- append(vec, c(rep(0, intermission), rep(dose, duration)))}
  return(vec[1:m])}
# create TD var
Xmat[,2] <- do.call("c", lapply(1:n, function(i) TDhist(m)))
Xmat[,3] <- do.call("c", lapply(1:n, function(i) TDhist(m)))

head(Xmat)

# genereate vectors of event and censoring times prior to calling the
# function for the algorithm
eventRandom <- round(rexp(n, 0.012)+1,0)
censorRandom <- round(runif(n, 1,870),0)
# Generate the survival data conditional on the three covariates

data <- permalgorithm(n, m, Xmat, XmatNames=c("sex", "Drug1", "Drug2"),
                      eventRandom = eventRandom, censorRandom=censorRandom, betas=c(log(2),
                                                                                    log(1.04), log(0.99)), groupByD=FALSE )
head(data)
# could use survival library and check whether the data was generated
# properly using coxph(Surv(Start, Stop, Event) ~ sex + Drug1 + Drug2,
# data)

test<-coxph(Surv(Start, Stop, Event) ~ sex + Drug1 + Drug2, data=data)


##So...my turn.
#Generate random survival times, random events; counting process expasion; add time-dependent element and use that for xmat in PermoAlgo.
rm(list = ls(all.names = TRUE))
n=50000
m=365

xmat<-matrix(nrow=n*m, ncol=11)

xmat[,1]<-rep(round(rbinom(n,1, 0.5)), each=m)
xmat[,2]<-rep(round(rbinom(n,1, 0.5)), each=m)
xmat[,3]<-rep(round(rbinom(n,1, 0.5)), each=m)
xmat[,4]<-rep(round(rbinom(n,1, 0.5)), each=m)

#5th column is going to be column number/time interval as if we had applied counting process to a data set
xmat[,5]<-rep.int(1:m, times=n)
xmat[,6]<-log(xmat[,5])*xmat[,1]
xmat[,7]<-log(xmat[,5])*xmat[,2]


xmat[,8]<-xmat[,5]*xmat[,5]*xmat[,1]
xmat[,9]<-xmat[,5]*xmat[,5]*xmat[,2]

xmat[,10]<-xmat[,1]*xmat[,3]
xmat[,11]<-xmat[,2]*xmat[,4]

dsMaster<-as.data.frame(xmat)
head(dsMaster)

colnames(dsMaster)<-c("X1", "X2", "X3", "X4", "intForFt", "logtX1", "logtX2", "tsquaredX1", "tsquaredX2", "X13", "X24")


#Now permute it!

dsTest<-as.matrix(dsMaster[,c(1:4, 6, 7)])
head(dsTest)

dataTest<-permalgorithm(n, m, Xmat=dsTest, XmatNames=c("X1", "X2", "X3", "X4", "logtX1", "logtX2"), betas=c(0.7, 0.7, 0.1, 0.1, 0.7, 0.7))

head(dataTest)
dataTest[121:300,]

test<-coxph(Surv(Start, Stop, Event) ~ X1 + X2 + X3 + X4 + logtX1 + logtX2, data=dataTest)


