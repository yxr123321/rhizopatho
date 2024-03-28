library(meta)
library(metafor)

richness=read.csv("data/all_shannon.csv",header=T)

d1<-escalc(measure="ROM",data=richness,m1i=average.y,sd1i=sd.y,n1i=n.y,m2i=average.x,sd2i=sd.x,n2i=n.x)

total1<-rma(yi,vi, data=d1, method="DL")
funnel(total1,main="Bacterial shannon")

### Begg's test#####
ranktest(total1)[2]


###Trim and fill method
###Trim and fill: A simple funnel-plot-based method 
###of testing and adjusting for publication bias 
###in meta-analysis. Biometrics, 56(2), 455–463. 

summary(total1)
t1<-trimfill(total1,estimator="R0")
funnel(t1)
summary(t1)