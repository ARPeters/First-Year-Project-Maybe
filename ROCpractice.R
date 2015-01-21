install.packages("pROC")
library(pROC)
install.packages("ROCR")
library(ROCR)

truth<-c(0,0,0,1,1,1,1,1)
predicted<-c(0.1,.5,.3,.8,.9,.4,.9,.5)

pred<- prediction(predicted, truth)
perf<- performance(pred, "tpr", "fpr")


#Tutorial from http://www.r-bloggers.com/roc-curves-and-classification/
db = read.table("http://freakonometrics.free.fr/db.txt",header=TRUE,sep=";")
head(db)
attach(db)


X3bis=rep(NA,length(X3))
X3bis[X3%in%c("A","C","D")]="ACD"
X3bis[X3%in%c("B","E")]="BE"
db$X3bis=as.factor(X3bis)

head(db)

reg=glm(Y~X1+X2+X3bis,family=binomial,data=db)

S=predict(reg, type="response")
?predict()
#Function that gives contingency table (I expect)
roc.curve=function(s,print=FALSE){
  Ps=(S>s)*1
  FP=sum((Ps==1)*(Y==0))/sum(Y==0)
  TP=sum((Ps==1)*(Y==1))/sum(Y==1)
  if(print==TRUE){
    print(table(Observed=Y,Predicted=Ps))
    }
  vect=c(FP,TP)
  names(vect)=c("FPR","TPR")
  return(vect)
  }
threshold = 0.5

#roc.curve for this does not match up with blog example. Continuing anyway, I guess.
roc.curve(threshold, print=TRUE)

#"For convenience"(?)
ROC.curve=Vectorize(roc.curve)

#plotting predicted vs observed
I=(((S>threshold)&(Y==0))|((S<=threshold)&(Y==1)))
plot(S,Y,col=c("red", "blue")[I+1], pch=19, cex=.7, xlab="", ylab="")

#Add vertical line, representing random guessing
abline(a=0, b=1,col="gray")

M.ROC=ROC.curve(seq(0,1,by=.01))
plot(M.ROC[1,], M.ROC[2,], col="grey", lwd=2, type="l", xlab="False Positive Rate", ylab="True Positive Rate")

#Same Blog, differet example
install.packages("tree")
library(tree)

ctr<-tree(Y~X1+X2+X3bis, data=db)
plot(ctr)
text(ctr)

S=predict(ctr)

#Now doing everything again, because we have a new prediction model S.
roc.curve(threshold, print=TRUE)

#"For convenience"(?)
ROC.curve=Vectorize(roc.curve)

#plotting predicted vs observed
I=(((S>threshold)&(Y==0))|((S<=threshold)&(Y==1)))
plot(S,Y,col=c("red", "blue")[I+1], pch=19, cex=.7, xlab="", ylab="")

#abline(v=seuil,col="gray")

M.ROC.tree=ROC.curve(seq(0,1,by=.01))
plot(M.ROC.tree[1,], M.ROC.tree[2,], col="grey", lwd=2, type="l", xlab="False Positive Rate", ylab="True Positive Rate")

#Plotting both on same graph
plot(M.ROC[1,], M.ROC[2,], col="grey", lwd=2, type="l", xlab="False Positive Rate", ylab="True Positive Rate")
lines(M.ROC.tree[1,], M.ROC.tree[2,], type="l", col="grey", lwd=2)

#############################################################################################
#Better Example (I think)
#http://www.r-bloggers.com/a-small-introduction-to-the-rocr-package/

#get ROCR package and simpel dataset
data(ROCR.simple)
ds<-cbind(ROCR.simple$predictions, ROCR.simple$labels)
colnames(ds)<-c("predictions", "labels")
ds<-as.data.frame(ds)
attach(ds)

#Create prediction object
pred<-prediction(predictions, labels)

class(pred)
slotNames(pred)

sn<-slotNames(pred)
sapply(sn, function(x) length(slot(pred, x)))

sapply(sn, function(x) class(slot(pred, x)))

#If interested, can do this same stuff with multiple predictions of the same label; some other time. 

#Creating performance object
perfobject<-performance(pred, measure="tpr", x.measure="fpr")
plot(perfobject)

#Add straight, diagonal line to graph; represents random guessing. 
abline(a=0, b=1)

#Graphing Accuracy (total proportion of correct predictions, (TP+TN)/(P+N)); Cutoffs on x axis, accuracy estimates on y.
accPerf<-performance(pred, measure="acc")
plot(accPerf)

#AUC, Area under the curve. 
aucPerf<-performance(pred, measure="auc")
aucPerf@y.values

#Partial AUC; Area under subset of curve defined by max FPR.
MaxFPR=0.1
paucPerf<-performance(pred, measure="auc", fpr.stop=MaxFPR)
paucPerf@y.values

#Standardize to give auc relative to its own area
spaucPerf<-as.numeric(paucPerf@y.values)/MaxFPR
