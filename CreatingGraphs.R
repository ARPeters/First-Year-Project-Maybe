rm(list = ls(all.names = TRUE))


dsH1000<-read.csv("SimH1000N324TP.csv")
dsH500<-read.csv("SimH500N325TP.csv")

dsInt1000<-read.csv("SimInt1000N325TP.csv")
dsInt500<-read.csv("SimInt500N325TP.csv")

dsLog1000<-read.csv("SimLog1000N325TP.csv")
dsLog500<-read.csv("SimLog500N329TP.csv")

dsT1000<-read.csv("SimT1000N329TP.csv")
dsT500<-read.csv("SimT500N329TP.csv")

#Preparing graphs
#Note: will need to rerun Heaviside 1000; missing time model

ftw<-as.vector(as.numeric(dsH1000[1,2:101]))
AICPropCH1000<-as.vector(as.numeric(dsH1000[3,2:101]))
BICPropCH1000<-as.vector(as.numeric(dsH1000[7,2:101]))
cvPropCH1000<-as.vector(as.numeric(dsH1000[11,2:101]))
AUCPropCH1000<-as.vector(as.numeric(dsH1000[15,2:101]))


ftw<-as.vector(as.numeric(dsH500[1,2:101]))
AICPropCH500<-as.vector(as.numeric(dsH500[3,2:101]))
BICPropCH500<-as.vector(as.numeric(dsH500[7,2:101]))
cvPropCH500<-as.vector(as.numeric(dsH500[11,2:101]))
AUCPropCH500<-as.vector(as.numeric(dsH500[15,2:101]))


plot(ftw, AICPropCH1000, col="red", ylab="Proportion of correct selections", 
     xlab="Regression coefficient of f(t)", ylim=range(c(AICPropCH1000,BICPropCH1000, AUCPropCH1000)),
     main="AIC vs BIC: Heaviside Function, N=1000")
par(new=TRUE)
plot(ftw, BICPropCH1000, xlab="", ylab="", ylim=range(c(AICPropCH1000,BICPropCH1000, AUCPropCH1000)))
legend(x=0, y=1, legend=c("AIC", "BIC"), col=par("red"), fill=FALSE, border=c("Red", "Black"))

plot(ftw, AICPropCH500, col="red", ylab="Proportion of correct selections", 
     xlab="Regression coefficient of f(t)", ylim=range(c(AICPropCH500,BICPropCH500)),
     main="AIC vs BIC: Heaviside Function, N=500")
par(new=TRUE)
plot(ftw, BICPropCH500, xlab="", ylab="", ylim=range(c(AICPropCH500,BICPropCH500)))


plot(ftw, cvPropCH1000)
plot(ftw, AUCPropCH1000)

plot(ftw, AUCPropCH1000, col="blue", ylim=range(c(AICPropCH1000,BICPropCH1000, AUCPropCH1000)),
     xlab="Regression coefficient of f(t)",
     ylab="Proportion of correct selections",
     main="AUC: Heaviside function, N=1000")


plot(ftw, AUCPropCH500, col="blue", ylim=range(c(AICPropCH500,BICPropCH500, AUCPropCH500)),
     xlab="Regression coefficient of f(t)",
     ylab="Proportion of correct selections",
     main="AUC: Heaviside function, N=500")

plot(ftw, cvPropCH1000, col="blue", ylim=range(c(AICPropCH1000,BICPropCH1000, cvPropCH1000)),
     xlab="Regression coefficient of f(t)",
     ylab="Proportion of correct selections",
     main="cv: Heaviside function, N=1000")


plot(ftw, cvPropCH500, col="blue", ylim=range(c(AICPropCH500,BICPropCH500, cvPropCH500)),
     xlab="Regression coefficient of f(t)",
     ylab="Proportion of correct selections",
     main="cv: Heaviside function, N=500")


plot(ftw, AUCPropCH1000, col="black", ylim=range(c(AICPropCH1000, AUCPropCH1000, AUCPropCH500)),
     xlab="Regression coefficient of f(t)",
     ylab="Proportion of correct selections",
     main="AUC: Heaviside Function")
par(new=TRUE)
plot(ftw, AUCPropCH500, col="red", ylim=range(c(AICPropCH1000, AUCPropCH1000, AUCPropCH500)),
     xlab="",
     ylab="",
     main="")
legend(x=.7, y=1, legend=c("N=1000", "N=500"),border=c("black", "red"), fill=FALSE )


###################################################
#Interaction
#Also do not include Time model for some reason

rm(list = ls(all.names = TRUE))

dsInt1000<-read.csv("SimInt1000N325TP.csv")
dsInt500<-read.csv("SimInt500N325TP.csv")

dsInt500[,1:2]

ftw<-as.vector(as.numeric(dsInt1000[1,2:101]))
AICPropCInt1000<-as.vector(as.numeric(dsInt1000[4,2:101]))
BICPropCInt1000<-as.vector(as.numeric(dsInt1000[8,2:101]))
cvPropCInt1000<-as.vector(as.numeric(dsInt1000[12,2:101]))
AUCPropCInt1000<-as.vector(as.numeric(dsInt1000[16,2:101]))


ftw<-as.vector(as.numeric(dsInt500[1,2:101]))
AICPropCInt500<-as.vector(as.numeric(dsInt500[4,2:101]))
BICPropCInt500<-as.vector(as.numeric(dsInt500[8,2:101]))
cvPropCInt500<-as.vector(as.numeric(dsInt500[12,2:101]))
AUCPropCInt500<-as.vector(as.numeric(dsInt500[16,2:101]))


plot(ftw, AICPropCInt1000, col="red", ylab="Proportion of correct selections", 
     xlab="Regression coefficient of Interaction Term", ylim=range(c(AICPropCInt1000,BICPropCInt1000)),
     main="AIC vs BIC: Interaction model, N=1000")
par(new=TRUE)
plot(ftw, BICPropCInt1000, xlab="", ylab="", ylim=range(c(AICPropCInt1000,BICPropCInt1000)))
legend(x=0, y=1, legend=c("AIC", "BIC"), col=par("red"), fill=FALSE, border=c("Red", "Black"))


plot(ftw, AICPropCInt500, col="red", ylab="Proportion of correct selections", 
     xlab="Regression coefficient of Interaction Term", ylim=range(c(AICPropCInt500,BICPropCInt500)),
     main="AIC vs BIC: Interaction Model, N=500")
par(new=TRUE)
plot(ftw, BICPropCInt500, xlab="", ylab="", ylim=range(c(AICPropCInt500,BICPropCInt500)))
legend(x=0, y=1, legend=c("AIC", "BIC"), col=par("red"), fill=FALSE, border=c("Red", "Black"))


###AUC comparison for Int
plot(ftw, AUCPropCInt1000, col="black", ylim=range(c(AICPropCInt1000,AUCPropCInt1000, AUCPropCInt500)),
     xlab="Regression coefficient of f(t)",
     ylab="Proportion of correct selections",
     main="AUC: Interaction")
par(new=TRUE)
plot(ftw, AUCPropCInt500, col="red", ylim=range(c(AICPropCInt1000,AUCPropCInt1000, AUCPropCInt500)),
     xlab="",
     ylab="",
     main="")
legend(x=.7, y=.2, legend=c("N=1000", "N=500"),border=c("black", "red"), fill=FALSE )


###################################################
#f(t)=t
#...The hell?
rm(list = ls(all.names = TRUE))


dsT1000<-read.csv("SimT1000N329TP.csv")
dsT500<-read.csv("SimT500N329TP.csv")

ftwT10<-as.vector(as.numeric(dsT1000[1,2:101]))
AICPropCT1000<-as.vector(as.numeric(dsT1000[6,2:101]))
BICPropCT1000<-as.vector(as.numeric(dsT1000[11,2:101]))
cvPropCT1000<-as.vector(as.numeric(dsT1000[16,2:101]))
AUCPropCT1000<-as.vector(as.numeric(dsT1000[21,2:101]))

dsT500[,1:2]
ftw5<-as.vector(as.numeric(dsT500[1,2:101]))
AICPropCT500<-as.vector(as.numeric(dsT500[6,2:101]))
BICPropCT500<-as.vector(as.numeric(dsT500[11,2:101]))
cvPropCT500<-as.vector(as.numeric(dsT500[16,2:101]))
AUCPropCT500<-as.vector(as.numeric(dsT500[21,2:101]))

dsT1000[,1:2]
dsT500[,1:2]

ftw<-c(1:100)/100

plot(ftw, AICPropCT1000, col="red", ylab="Proportion of correct selections", 
     xlab="Regression coefficient of Time-Covariate Interaction", ylim=range(c(AICPropCT1000,BICPropCT1000)),
     main="AIC vs BIC: f(t)=t, N=1000")
par(new=TRUE)
plot(ftw, BICPropCT1000, xlab="", ylab="", ylim=range(c(AICPropCT1000,BICPropCT1000)))
legend(x=0, y=1, legend=c("AIC", "BIC"), col=par("red"), fill=FALSE, border=c("Red", "Black"))


plot(ftw, AICPropCT500, col="red", ylab="Proportion of correct selections", 
     xlab="Regression coefficient of Time-Covariate Interaction", ylim=range(c(AICPropCT500,BICPropCT500)),
     main="AIC vs BIC: f(t)=t, N=500")
par(new=TRUE)
plot(ftw, BICPropCT500, xlab="", ylab="", ylim=range(c(AICPropCT500,BICPropCT500)))
legend(x=0, y=1, legend=c("AIC", "BIC"), col=par("red"), fill=FALSE, border=c("Red", "Black"))

###################################################
#f(t)=log(t)

rm(list = ls(all.names = TRUE))


dsLog1000<-read.csv("SimLog1000N325TP.csv")
dsLog500<-read.csv("SimLog500N329TP.csv")

dsLog1000[,1:2]
dsLog500[,1:2]

ftw<-as.vector(as.numeric(dsLog1000[1,2:101]))
AICPropCLog1000<-as.vector(as.numeric(dsLog1000[5,2:101]))
BICPropCLog1000<-as.vector(as.numeric(dsLog1000[9,2:101]))
cvPropCLog1000<-as.vector(as.numeric(dsLog1000[13,2:101]))
AUCPropCLog1000<-as.vector(as.numeric(dsLog1000[17,2:101]))


ftw<-as.vector(as.numeric(dsLog500[1,2:101]))
AICPropCLog500<-as.vector(as.numeric(dsLog500[5,2:101]))
BICPropCLog500<-as.vector(as.numeric(dsLog500[10,2:101]))
cvPropCLog500<-as.vector(as.numeric(dsLog500[15,2:101]))
AUCPropCLog500<-as.vector(as.numeric(dsLog500[20,2:101]))


plot(ftw, AICPropCLog1000, col="red", ylab="Proportion of correct selections", 
     xlab="Regression coefficient of f(t)", ylim=range(c(AICPropCLog1000,BICPropCLog1000)),
     main="AIC vs BIC: f(t)=log, N=1000")
par(new=TRUE)
plot(ftw, BICPropCLog1000, xlab="", ylab="", ylim=range(c(AICPropCLog1000,BICPropCLog1000)))
legend(x=0, y=1, legend=c("AIC", "BIC"), col=par("red"), fill=FALSE, border=c("Red", "Black"))


plot(ftw, AICPropCLog500, col="red", ylab="Proportion of correct selections", 
     xlab="Regression coefficient of f(t)", ylim=range(c(AICPropCLog500,BICPropCLog500)),
     main="AIC vs BIC: f(t)=log, N=500")
par(new=TRUE)
plot(ftw, BICPropCLog500, xlab="", ylab="", ylim=range(c(AICPropCLog500,BICPropCLog500)))

###AUC comparison for Log
plot(ftw, AUCPropCLog1000, col="black", ylim=range(c(AICPropCLog1000,AUCPropCLog1000, AUCPropCLog500)),
     xlab="Regression coefficient of f(t)",
     ylab="Proportion of correct selections",
     main="AUC: log(t)")
par(new=TRUE)
plot(ftw, AUCPropCLog500, col="red", ylim=range(c(AICPropCLog1000,AUCPropCLog1000, AUCPropCLog500)),
     xlab="",
     ylab="",
     main="")
legend(x=.7, y=1, legend=c("N=1000", "N=500"),border=c("black", "red"), fill=FALSE )

########cv
plot(ftw, cvPropCLog1000, col="blue", ylim=range(c(AICPropCLog1000,BICPropCLog1000, cvPropCLog1000)),
     xlab="Regression coefficient of f(t)",
     ylab="Proportion of correct selections",
     main="cv: log(t) function, N=1000")


plot(ftw, cvPropCLog500, col="blue", ylim=range(c(AICPropCLog500,BICPropCLog500, cvPropCLog500)),
     xlab="Regression coefficient of f(t)",
     ylab="Proportion of correct selections",
     main="cv: log(t) function, N=500")

