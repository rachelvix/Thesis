####Binomial/Poisson N-mixture model to estimate Clapper Rail abundance across Coastal MS -----------------------------------------------------------------------------------------------------
#General approach adapted from Ke?ry & Royle's Applied Hierarchical Modeling in Ecology
#Ch. 6 Binomial N-mixture Models 

#This code incldudes the random effect of site 

##Assumptions
#1. Closed population 
#2. No false-positives 
#3. Independence of detection 
#4. Homogeneity of detection 
#5. Parametric modeling assumptions 

setwd("C:/Users/rva21/Desktop/Analysis/Data")


#TMB installation in this order:
#install.packages("Matrix")
#install.packages("Rcpp")
#install.packages("RcppEigen")
#install.packages("TMB", type="source")
#install.packages("glmmTMB", type="source")


#Load packages
library(TMB)
library(glmmTMB)
library(unmarked)
library(AICcmodavg)
library(ggplot2)
library(nmixgof)


##Read data 
#Clapper rail (CLRA) counts and covariates in wide format 
#262 point count sites in total
#3 visits in year 2021 & 2022 
CLRA.w <- read.table("pc.data.csv", header=T, sep=",")   
head(CLRA.w)
str(CLRA.w)

#select the first year 
covs <- CLRA.w[1:262,]

#Counts with truncation and passive listening only 
CLRAc <- as.matrix(covs[,c("CLRA1", "CLRA2", "CLRA3")])

#Julian day matrix 
daymat <- as.matrix(covs[,c("day.1", "day.2","day.3")])
#scale
daymat <- scale(daymat)

#Wind matrix
covs$wind.1 <- as.factor(covs[,"wind.1"])
covs$wind.2 <- as.factor(covs[,"wind.2"])
covs$wind.3 <- as.factor(covs[,"wind.3"])
windmat <- as.matrix(covs[,c("wind.1", "wind.2","wind.3")]) 

#Visit matrix
covs$v1 <- rep(as.character(1))
covs$v2 <- rep(as.character(2))
covs$v3 <- rep(as.character(3))
visit <- as.matrix(cbind(covs[,c("v1", "v2", "v3")]))

#Treat Point ID as factor
covs$Point_ID <- as.factor(covs$Point_ID)

#Treat observer as factor 
covs$obs.1 <- as.factor(covs$obs.1)
covs$obs.2 <- as.factor(covs$obs.2)
covs$obs.3 <- as.factor(covs$obs.3)
obs <- as.matrix(cbind(covs[,c("obs.1", "obs.2", "obs.3")]))

##===========================================================================================
##Model with just local variables 

#Scale local veg covs
covs[,"SP_Percentage"] <- scale(covs[,"SP_Percentage"])
covs[,"DS_Percentage"] <- scale(covs[,"DS_Percentage"])
covs[,"Height"] <- scale(covs[,"Height"])
covs[,"Fresh"] <- scale(covs[,"Fresh"])

##==========================================================================================
##Create unmarked df------------------------------------------------------------------------ 

umfc.l <- unmarkedFramePCount(y=CLRAc, obsCovs = list(date=(daymat), wind=(windmat), visit=(visit), obs=(obs)), 
                              siteCovs = covs[,c("Point_ID", "DS_Percentage", "SP_Percentage", "Height", "Fresh")])
summary(umfc.l)
summary(apply(CLRAc,1,max,na.rm=TRUE))
#mean counts per point= 3.282
#max counts per point= 20

#Global model 
fmc.l<- pcount(~date+wind+obs+(1|visit) ~DS_Percentage+SP_Percentage+Height+(1|Point_ID), umfc.l)
summary(fmc.l)

#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.743    0.862

#Fixed effects:
#  Estimate     SE     z   P(>|z|)
#(Intercept)      2.831 0.1268 22.33 1.81e-110
#DS_Percentage    0.140 0.0714  1.96  5.06e-02
#SP_Percentage    0.270 0.0656  4.11  3.88e-05
#Height          -0.239 0.0689 -3.47  5.29e-04

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.052    0.228

#Fixed effects:
#  Estimate     SE      z  P(>|z|)
#(Intercept)          -3.372 0.6494 -5.192 2.08e-07
#date                 -0.207 0.0585 -3.531 4.14e-04
#wind1                -0.192 0.0894 -2.151 3.15e-02
#wind2                -0.159 0.0992 -1.604 1.09e-01
#wind3                -0.358 0.1379 -2.598 9.37e-03
#wind4                -0.483 0.1471 -3.285 1.02e-03
#wind5                -0.904 0.3045 -2.971 2.97e-03
#obsMattDunning        1.027 0.6317  1.626 1.04e-01
#obsRachel Anderson    1.332 0.6419  2.075 3.80e-02
#obsSofia Campuzano    0.607 0.6461  0.939 3.48e-01

#AIC: 2845.315 



##=================================================================================
##Model Building ------------------------------------------------------------------
#Detection covs first  
mlist = list()

mlist$d1 = pcount(~date ~1, umfc, engine = "TMB")
mlist$d2 = pcount(~wind ~1, umfc, engine = "TMB")
mlist$d3 = pcount(~(1|visit) ~1, umfc, engine = "TMB")
mlist$d4 = pcount(~date+wind ~1, umfc, engine = "TMB")
mlist$d5 = pcount(~date+(1|visit) ~1, umfc, engine = "TMB")
mlist$d6 = pcount(~wind+(1|visit) ~1, umfc, engine = "TMB")
mlist$d7 = pcount(~date+wind+(1|visit) ~1, umfc, engine = "TMB")
mlist$d8 = pcount(~obs ~1, umfc, engine = "TMB")
mlist$d9 = pcount(~date+obs ~1, umfc, engine = "TMB")
mlist$d10 = pcount(~(1|visit)+obs ~1, umfc, engine = "TMB")
mlist$d11 = pcount(~wind+obs ~1, umfc, engine = "TMB")
mlist$d12 = pcount(~date+wind+obs ~1, umfc, engine = "TMB")
mlist$d13 = pcount(~date+(1|visit)+obs ~1, umfc, engine = "TMB")
mlist$d14 = pcount(~wind+(1|visit)+obs ~1, umfc, engine = "TMB")
mlist$d15 = pcount(~date+wind+(1|visit)+obs ~1, umfc, engine = "TMB")

fl1 <- modSel(fitList(fits = mlist)) # Compare models. Move forward with best model(s) (delta AIC < 2)
#nPars     AIC delta   AICwt cumltvWt
#d15    11 3574.80   0.00 6.0e-01     0.60
#d12    11 3575.65   0.85 3.9e-01     0.99

summary(mlist$d15)
##add abundance covs to best model 
cm1<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage, data=umfc.l, mixture= "P", engine = "TMB")
cm2<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage, data=umfc.l, mixture= "P", engine = "TMB")
cm3<- pcount(~date+wind+(1|visit)+obs ~Height, data=umfc.l, mixture= "P", engine = "TMB")
cm4<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage, data=umfc.l, mixture= "P", engine = "TMB")
cm5<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Height, data=umfc.l, mixture= "P", engine = "TMB")
cm6<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + Height, data=umfc.l, mixture= "P", engine = "TMB")
cm7<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height, data=umfc.l, mixture= "P", engine = "TMB")
cm8<- pcount(~date+wind+(1|visit)+obs ~(1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm9<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + (1|Point_ID),  data=umfc.l, mixture= "P", engine = "TMB")
cm10<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + (1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm11<- pcount(~date+wind+(1|visit)+obs ~Height + (1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm12<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + (1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm13<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Height + (1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm14<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + Height + (1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm15<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height + (1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm16<- pcount(~date+wind+(1|visit)+obs ~Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm17<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm18<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm19<- pcount(~date+wind+(1|visit)+obs ~Height + Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm20<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height +Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm21<- pcount(~date+wind+(1|visit)+obs ~(1|Point_ID)+Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm22<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + (1|Point_ID) +Fresh,  data=umfc.l, mixture= "P", engine = "TMB")
cm23<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + (1|Point_ID) +Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm24<- pcount(~date+wind+(1|visit)+obs ~Height + (1|Point_ID) +Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm25<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + (1|Point_ID)+Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm26<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Height + (1|Point_ID) +Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm27<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + Height + (1|Point_ID) +Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm28<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height + (1|Point_ID) +Fresh, data=umfc.l, mixture= "P", engine = "TMB")

#Create fit list of all models
cmL <- modSel(fitList(cm1, cm2, cm3, cm4, cm5, cm6, cm7, cm8, cm9, cm10, 
                      cm11, cm12, cm13, cm14, cm15, cm16, cm17, cm18, cm19, cm20, cm21, cm22, cm23, cm24, cm25, cm26, cm27, cm28))

cmL

#nPars     AIC  delta    AICwt cumltvWt
#nPars     AIC  delta   AICwt cumltvWt
#cm26    14 2824.96   0.00 6.0e-01     0.60
#cm28    15 2826.41   1.45 2.9e-01     0.89

#Best Models---------------------------------------------------------------------------------
cm26<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Height + (1|Point_ID) +Fresh, data=umfc.l,
              mixture= "P", engine = "TMB")
summary(cm26)

#Random effects:
#Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.678    0.823

#Fixed effects:
#  Estimate     SE     z   P(>|z|)
#(Intercept)      2.759 0.1267 21.77 4.43e-105
#SP_Percentage    0.161 0.0643  2.51  1.21e-02
#Height          -0.222 0.0661 -3.35  8.02e-04
#Fresh           -0.416 0.0875 -4.75  2.04e-06

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.051    0.226

#Fixed effects:
#  Estimate     SE     z  P(>|z|)
#(Intercept)          -3.367 0.6488 -5.19 2.10e-07
#date                 -0.186 0.0563 -3.30 9.83e-04
#wind1                -0.196 0.0894 -2.19 2.84e-02
#wind2                -0.168 0.0990 -1.70 8.97e-02
#wind3                -0.368 0.1358 -2.71 6.79e-03
#wind4                -0.487 0.1458 -3.34 8.33e-04
#wind5                -0.913 0.3034 -3.01 2.62e-03
#obsMattDunning        1.040 0.6313  1.65 9.96e-02
#obsRachel Anderson    1.296 0.6410  2.02 4.32e-02
#obsSofia Campuzano    0.794 0.6467  1.23 2.20e-01

#AIC: 2824.959 


#ZIP
cm26z<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Height + (1|Point_ID) +Fresh, data=umfc.l,
               mixture= "ZIP", engine = "TMB")
summary(cm26z)

#Abundance (log-scale):
#Random effects:
#Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.466    0.682

#Fixed effects:
#  Estimate     SE     z   P(>|z|)
#(Intercept)     3.0423 0.1228 24.78 1.62e-135
#SP_Percentage   0.0778 0.0549  1.42  1.57e-01    ns 
#Height         -0.1687 0.0414 -4.08  4.55e-05
#Fresh          -0.4185 0.0682 -6.14  8.31e-10

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)     0.05    0.224

#Fixed effects:
#  Estimate     SE     z  P(>|z|)
#(Intercept)          -3.378 0.6499 -5.20 2.01e-07
#date                 -0.151 0.0384 -3.94 8.10e-05
#wind1                -0.197 0.0829 -2.38 1.75e-02
#wind2                -0.168 0.0888 -1.89 5.85e-02
#wind3                -0.367 0.1191 -3.08 2.06e-03
#wind4                -0.486 0.1423 -3.42 6.33e-04
#wind5                -0.941 0.2981 -3.16 1.60e-03
#obsMattDunning        1.034 0.6319  1.64 1.02e-01
#obsRachel Anderson    1.238 0.6394  1.94 5.28e-02
#obsSofia Campuzano    0.693 0.6395  1.08 2.78e-01

#Zero-inflation (logit-scale):
#  Estimate    SE     z  P(>|z|)
#-2.3 0.185 -12.5 1.28e-35

#AIC: 2805.827 


##ZIP is best! 

#take out non sig 
cm26z<- pcount(~date+wind+(1|visit)+obs ~Height + (1|Point_ID) +Fresh, data=umfc.l,
               mixture= "ZIP", engine = "TMB")
summary(cm26z)   
#NAs 

##====================================================================================================
##Abundance prediction-----------------------------------------------------------------------------

#Poisson 
renb <- ranef(cm26)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#6110.204 5979.975 6250.025 
 


##Z
renb <- ranef(cm26z)
plot(renb, layout=c(4,3), subset=site%in%1:12, xlim=c(-1, 11), lwd=5)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
##Estimate     2.5%    97.5% 
#6229.211 6085.875 6380.000 




##=========================================================================
##Gof----------------------------------------------------------------------

#Plots 
#Residuals
plot(cm26) 
plot(cm26z) 


#expected vs observed counts 
n.obs <- apply(CLRAc, 1, sum, na.rm=TRUE)
n.predP <- apply(fitted(cm26), 1, sum, na.rm=TRUE )
n.predZ <- apply(fitted(cm26z), 1, sum, na.rm=TRUE )
plot(n.obs, n.predP, frame=F, ylab="Expected Counts", xlab="Observed Counts", main="CLRA local pcount")
abline(0,1)
points(smooth.spline(n.predP ~ n.obs, df=4), type="l", lwd=2, col="blue")
points(smooth.spline(n.predZ ~ n.obs, df=4), type="l", lwd=2, col="red")

#CLRApc.local.ove

#fitted vs residuals 
n.residsZ <-apply(residuals(cm26z),1,sum,na.rm=T)
plot(n.predZ, n.residsZ, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="CLRA stacked pcount NB")
points(smooth.spline(n.residsZ ~ n.predZ, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predZ, df=4), type="l", lwd=2, col="black")
#CLRApc.local.fresid

#fitted vs residuals 
n.residsP <-apply(residuals(cm26),1,sum,na.rm=T)
plot(n.predP, n.residsP, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="CLRA stacked pcount NB")
points(smooth.spline(n.residsP ~ n.predP, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predP, df=4), type="l", lwd=2, col="black")
#CLRApc.local.fresid

##Nmixgof----------------------------------------------------------------------------------
#Chat 
chat(cm26, type="marginal") #0.9515996
chat(cm26z, type="marginal") #0.9398458


#residuals 
residqq(cm26, type="site-sum")
residqq(cm26z, type="site-sum")

residqq(cm26, type="observation")
residqq(cm26z, type="observation")



##Parametric bootstrap------------------------------------------------------------------- 
#Fistats function (AHM book) 
fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2)
  chisq <- sum((observed - expected)^2 / expected)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

#Best model 
#NB 
gofp.l <- parboot(cm26, fitstats, nsim= 1000)
gofp.l
#  t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
#SSE          1841          791.5            63.88            0
#Chisq         736          209.4            20.29            0
#freemanTukey  361           99.2             8.61            0


##==============================================================================
##Cov Plots 

#Take out random effects 
cm26z<- pcount(~date+wind+obs ~SP_Percentage + Height +Fresh, data=umfc.l,
               mixture= "ZIP", engine = "TMB")
summary(cm26z)

#Height 
range(siteCovs(umfc.l)$Height) #~39-188

hplot <- data.frame(Height=seq(min(covs[,c("Height")]), max(covs[,c("Height")]), 
                               length=262), DS_Percentage=0, SP_Percentage=0, Fresh=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm26z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Height,y=Predicted,ymin=lower,ymax=upper))

cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Height,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Height,y=upper), alpha=0.1)

cwwi.plot <- cwwi.plot + geom_line(aes(x=Height,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Vegetation Height (cm)")+
  scale_y_continuous(limits = c(0,130), expand = c(0,0))+
  scale_x_continuous(expand=c(0,0))

h21.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(h21.plot) 


#Fresh
range(siteCovs(umfc.l)$Fresh) #0-100

hplot <- data.frame(Fresh=seq(min(covs[,c("Fresh")]), max(covs[,c("Fresh")]), 
                               length=262), DS_Percentage=0, Height=0, SP_Percentage=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm26z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Fresh,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Fresh,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Fresh,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Fresh,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Freshwater Vegetation (%)")+
scale_y_continuous(limits = c(0,130), expand = c(0,0))+
  scale_x_continuous(expand=c(0,0))

f21.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(t=5, r=8, b=5))

plot(f21.plot) 







###============================================================================
#year 2022
##Read data 
#Clapper rail (CLRA) counts and covariates in wide format 
#262 point count sites in total
#3 visits in year 2021 & 2022 
CLRA.w <- read.table("pc.data.csv", header=T, sep=",")   
head(CLRA.w)
str(CLRA.w)

#select the first year 
covs <- CLRA.w[263:524,]

#Counts with truncation and passive listening only 
CLRAc <- as.matrix(covs[,c("CLRA1", "CLRA2", "CLRA3")])

#Julian day matrix 
daymat <- as.matrix(covs[,c("day.1", "day.2","day.3")])
#scale
daymat <- scale(daymat)

#Wind matrix
covs$wind.1 <- as.factor(covs[,"wind.1"])
covs$wind.2 <- as.factor(covs[,"wind.2"])
covs$wind.3 <- as.factor(covs[,"wind.3"])
windmat <- as.matrix(covs[,c("wind.1", "wind.2","wind.3")]) 

#Visit matrix
covs$v1 <- rep(as.character(1))
covs$v2 <- rep(as.character(2))
covs$v3 <- rep(as.character(3))
visit <- as.matrix(cbind(covs[,c("v1", "v2", "v3")]))

#Treat Point ID as factor
covs$Point_ID <- as.factor(covs$Point_ID)

#Treat observer as factor 
covs$obs.1 <- as.factor(covs$obs.1)
covs$obs.2 <- as.factor(covs$obs.2)
covs$obs.3 <- as.factor(covs$obs.3)
obs <- as.matrix(cbind(covs[,c("obs.1", "obs.2", "obs.3")]))

##===========================================================================================
##Model with just local variables 

#Scale local veg covs
covs[,"SP_Percentage"] <- scale(covs[,"SP_Percentage"])
covs[,"DS_Percentage"] <- scale(covs[,"DS_Percentage"])
covs[,"Height"] <- scale(covs[,"Height"])
covs[,"Fresh"] <- scale(covs[,"Fresh"])

##==========================================================================================
##Create unmarked df------------------------------------------------------------------------ 

umfc.l <- unmarkedFramePCount(y=CLRAc, obsCovs = list(date=(daymat), wind=(windmat), visit=(visit), obs=(obs)), 
                              siteCovs = covs[,c("Point_ID", "DS_Percentage", "SP_Percentage", "Height", "Fresh")])
summary(umfc.l)
summary(apply(CLRAc,1,max,na.rm=TRUE))
#mean counts per point= 3.626
#max counts per point= 14

#Global model 
fmc.l<- pcount(~date+wind+obs+(1|visit) ~DS_Percentage+SP_Percentage+Height+(1|Point_ID), umfc.l)
summary(fmc.l)

#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.743    0.862

#Fixed effects:
#  Estimate     SE     z   P(>|z|)
#(Intercept)      2.831 0.1268 22.33 1.81e-110
#DS_Percentage    0.140 0.0714  1.96  5.06e-02
#SP_Percentage    0.270 0.0656  4.11  3.88e-05
#Height          -0.239 0.0689 -3.47  5.29e-04

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.052    0.228

#Fixed effects:
#  Estimate     SE      z  P(>|z|)
#(Intercept)          -3.372 0.6494 -5.192 2.08e-07
#date                 -0.207 0.0585 -3.531 4.14e-04
#wind1                -0.192 0.0894 -2.151 3.15e-02
#wind2                -0.159 0.0992 -1.604 1.09e-01
#wind3                -0.358 0.1379 -2.598 9.37e-03
#wind4                -0.483 0.1471 -3.285 1.02e-03
#wind5                -0.904 0.3045 -2.971 2.97e-03
#obsMattDunning        1.027 0.6317  1.626 1.04e-01
#obsRachel Anderson    1.332 0.6419  2.075 3.80e-02
#obsSofia Campuzano    0.607 0.6461  0.939 3.48e-01

#AIC: 2845.315 



##=================================================================================
##Model Building ------------------------------------------------------------------
#Detection covs first  
mlist = list()

mlist$d1 = pcount(~date ~1, umfc, engine = "TMB")
mlist$d2 = pcount(~wind ~1, umfc, engine = "TMB")
mlist$d3 = pcount(~(1|visit) ~1, umfc, engine = "TMB")
mlist$d4 = pcount(~date+wind ~1, umfc, engine = "TMB")
mlist$d5 = pcount(~date+(1|visit) ~1, umfc, engine = "TMB")
mlist$d6 = pcount(~wind+(1|visit) ~1, umfc, engine = "TMB")
mlist$d7 = pcount(~date+wind+(1|visit) ~1, umfc, engine = "TMB")
mlist$d8 = pcount(~obs ~1, umfc, engine = "TMB")
mlist$d9 = pcount(~date+obs ~1, umfc, engine = "TMB")
mlist$d10 = pcount(~(1|visit)+obs ~1, umfc, engine = "TMB")
mlist$d11 = pcount(~wind+obs ~1, umfc, engine = "TMB")
mlist$d12 = pcount(~date+wind+obs ~1, umfc, engine = "TMB")
mlist$d13 = pcount(~date+(1|visit)+obs ~1, umfc, engine = "TMB")
mlist$d14 = pcount(~wind+(1|visit)+obs ~1, umfc, engine = "TMB")
mlist$d15 = pcount(~date+wind+(1|visit)+obs ~1, umfc, engine = "TMB")

fl1 <- modSel(fitList(fits = mlist)) # Compare models. Move forward with best model(s) (delta AIC < 2)
#nPars     AIC delta   AICwt cumltvWt
#d15    11 3574.80   0.00 6.0e-01     0.60
#d12    11 3575.65   0.85 3.9e-01     0.99

##add abundance covs to best model 
cm1<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage, data=umfc.l, mixture= "P", engine = "TMB")
cm2<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage, data=umfc.l, mixture= "P", engine = "TMB")
cm3<- pcount(~date+wind+(1|visit)+obs ~Height, data=umfc.l, mixture= "P", engine = "TMB")
cm4<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage, data=umfc.l, mixture= "P", engine = "TMB")
cm5<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Height, data=umfc.l, mixture= "P", engine = "TMB")
cm6<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + Height, data=umfc.l, mixture= "P", engine = "TMB")
cm7<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height, data=umfc.l, mixture= "P", engine = "TMB")
cm8<- pcount(~date+wind+(1|visit)+obs ~(1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm9<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + (1|Point_ID),  data=umfc.l, mixture= "P", engine = "TMB")
cm10<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + (1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm11<- pcount(~date+wind+(1|visit)+obs ~Height + (1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm12<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + (1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm13<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Height + (1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm14<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + Height + (1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm15<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height + (1|Point_ID), data=umfc.l, mixture= "P", engine = "TMB")
cm16<- pcount(~date+wind+(1|visit)+obs ~Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm17<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm18<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm19<- pcount(~date+wind+(1|visit)+obs ~Height + Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm20<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height +Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm21<- pcount(~date+wind+(1|visit)+obs ~(1|Point_ID)+Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm22<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + (1|Point_ID) +Fresh,  data=umfc.l, mixture= "P", engine = "TMB")
cm23<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + (1|Point_ID) +Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm24<- pcount(~date+wind+(1|visit)+obs ~Height + (1|Point_ID) +Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm25<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + (1|Point_ID)+Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm26<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Height + (1|Point_ID) +Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm27<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + Height + (1|Point_ID) +Fresh, data=umfc.l, mixture= "P", engine = "TMB")
cm28<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height + (1|Point_ID) +Fresh, data=umfc.l, mixture= "P", engine = "TMB")

#Create fit list of all models
cmL <- modSel(fitList(cm1, cm2, cm3, cm4, cm5, cm6, cm7, cm8, cm9, cm10, 
                      cm11, cm12, cm13, cm14, cm15, cm16, cm17, cm18, cm19, cm20, cm21, cm22, cm23, cm24, cm25, cm26, cm27, cm28))

cmL

#nPars     AIC  delta   AICwt cumltvWt
#cm28    15 2999.06   0.00  7.2e-01     0.72
#cm27    14 3001.57   2.50  2.1e-01     0.93

#Best Models---------------------------------------------------------------------------------
cm28<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height + (1|Point_ID) +Fresh, 
              data=umfc.l, mixture= "P", engine = "TMB")

summary(cm28)
#Abundance (log-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)     0.83    0.911

#Fixed effects:
#  Estimate     SE     z   P(>|z|)
#(Intercept)      2.792 0.1236 22.58 6.44e-113
#SP_Percentage    0.166 0.0778  2.14  3.25e-02
#DS_Percentage    0.236 0.0854  2.77  5.63e-03
#Height          -0.386 0.0855 -4.52  6.28e-06
#Fresh           -0.423 0.1011 -4.18  2.94e-05

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.005    0.072

#Fixed effects:
#  Estimate      SE      z P(>|z|)
#(Intercept)         -9.6627 39.8126 -0.243  0.8082
#date                -0.1154  0.0546 -2.113  0.0346
#wind1                0.0686  0.0841  0.816  0.4147
#wind2                0.0798  0.0947  0.843  0.3993
#wind3               -0.0602  0.1237 -0.487  0.6265
#wind4               -0.3049  0.1481 -2.058  0.0396
#wind5               -0.8345  0.3516 -2.374  0.0176
#obsJamie Russell     7.0974 39.8126  0.178  0.8585
#obsRachel Anderson   7.1271 39.8123  0.179  0.8579
#obsTyler Connell     6.9454 39.8124  0.174  0.8615

#AIC: 2999.063 


#ZIP
cm28z<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height + (1|Point_ID) +Fresh, 
               data=umfc.l, mixture= "ZIP", engine = "TMB")
summary(cm28z)

#Fixed effects:
#Estimate     SE     z   P(>|z|)
#(Intercept)      2.995 0.1376 21.76 5.04e-105
#SP_Percentage    0.179 0.0729  2.46  1.40e-02
#DS_Percentage    0.257 0.0749  3.44  5.88e-04
#Height          -0.306 0.0703 -4.36  1.29e-05
#Fresh           -0.376 0.1093 -3.44  5.91e-04

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.005    0.072

#Fixed effects:
#  Estimate       SE       z P(>|z|)
#(Intercept)        -12.9997 202.3162 -0.0643  0.9488
#date                -0.0904   0.0615 -1.4700  0.1416
#wind1                0.0524   0.0877  0.5975  0.5502
#wind2                0.0761   0.1006  0.7566  0.4493
#wind3               -0.0908   0.1287 -0.7055  0.4805
#wind4               -0.3231   0.1520 -2.1262  0.0335
#wind5               -0.8743   0.3411 -2.5632  0.0104
#obsJamie Russell    10.4466 202.3192  0.0516  0.9588
#obsRachel Anderson  10.5396 202.3200  0.0521  0.9585
#obsTyler Connell    10.3403 202.3144  0.0511  0.9592

#Zero-inflation (logit-scale):
#  Estimate    SE     z  P(>|z|)
#-2.19 0.237 -9.25 2.24e-20

#AIC: 2983.918 


##ZIP is best! 


##====================================================================================================
##Abundance prediction-----------------------------------------------------------------------------

#Poisson 
renb <- ranef(cm28)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#7442.100 7297.975 7599.025 




##Z
renb <- ranef(cm28z)
plot(renb, layout=c(4,3), subset=site%in%1:12, xlim=c(-1, 11), lwd=5)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
##Estimate     2.5%    97.5% 
#6758.558 6605.000 6921.000 





##=========================================================================
##Gof----------------------------------------------------------------------

#Plots 
#Residuals
plot(cm28) 
plot(cm28z) 


#expected vs observed counts 
n.obs <- apply(CLRAc, 1, sum, na.rm=TRUE)
n.predP <- apply(fitted(cm28), 1, sum, na.rm=TRUE )
n.predZ <- apply(fitted(cm28z), 1, sum, na.rm=TRUE )
plot(n.obs, n.predP, frame=F, ylab="Expected Counts", xlab="Observed Counts", main="CLRA local pcount")
abline(0,1)
points(smooth.spline(n.predP ~ n.obs, df=4), type="l", lwd=2, col="blue")
points(smooth.spline(n.predZ ~ n.obs, df=4), type="l", lwd=2, col="red")

#CLRApc.local.ove

#fitted vs residuals 
n.residsZ <-apply(residuals(cm28z),1,sum,na.rm=T)
plot(n.predZ, n.residsZ, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="CLRA stacked pcount NB")
points(smooth.spline(n.residsZ ~ n.predZ, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predZ, df=4), type="l", lwd=2, col="black")
#CLRApc.local.fresid

#fitted vs residuals 
n.residsP <-apply(residuals(cm28),1,sum,na.rm=T)
plot(n.predP, n.residsP, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="CLRA stacked pcount NB")
points(smooth.spline(n.residsP ~ n.predP, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predP, df=4), type="l", lwd=2, col="black")
#CLRApc.local.fresid

##Nmixgof----------------------------------------------------------------------------------
#Chat 
chat(cm28, type="marginal") #0.9515996
chat(cm28z, type="marginal") #0.9398458


#residuals 
residqq(cm28, type="site-sum")
residqq(cm28z, type="site-sum")

residqq(cm28, type="observation")
residqq(cm28z, type="observation")



##Parametric bootstrap------------------------------------------------------------------- 
#Fistats function (AHM book) 
fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2)
  chisq <- sum((observed - expected)^2 / expected)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

#Best model 
#NB 
gofp.l <- parboot(cm28, fitstats, nsim= 1000)
gofp.l
#  t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
#SSE          1841          791.5            63.88            0
#Chisq         736          209.4            20.29            0
#freemanTukey  361           99.2             8.61            0




##==============================================================================
##Covs Plots 

#Take out random effects
cm28z<- pcount(~date+wind+obs ~SP_Percentage + DS_Percentage + Height +Fresh, 
               data=umfc.l, mixture= "ZIP", engine = "TMB")

#Height 
range(siteCovs(umfc.l)$Height) #~60-186

hplot <- data.frame(Height=seq(min(covs[,c("Height")]), max(covs[,c("Height")]), 
                               length=262), DS_Percentage=0, SP_Percentage=0, Fresh=0, 
                    Point_ID=0, visit=0, date=0, wind=0, observer=0)

h.pred <- predict(cm28z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Height,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Height,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Height,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Height,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Vegetation Height (cm)")+
  scale_y_continuous(limits = c(0,130), expand = c(0,0))+
  scale_x_continuous(expand=c(0,0))

h22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(h22.plot) 


#Fresh
range(siteCovs(umfc.l)$Fresh) #0-95

hplot <- data.frame(Fresh=seq(min(covs[,c("Fresh")]), max(covs[,c("Fresh")]), 
                              length=262), DS_Percentage=0, Height=0, SP_Percentage=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm28z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Fresh,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Fresh,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Fresh,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Fresh,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Freshwater Vegetation (%)")+
scale_y_continuous(limits = c(0,130), expand = c(0,0))+
  scale_x_continuous(expand=c(0,0))

f22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(t=5, r=8, b=5))

plot(f22.plot) 




#SP_Percentage
range(siteCovs(umfc.l)$SP_Percentage) #0-95

hplot <- data.frame(SP_Percentage=seq(min(covs[,c("SP_Percentage")]), max(covs[,c("SP_Percentage")]), 
                              length=262), DS_Percentage=0, Height=0, Fresh=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm28z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=SP_Percentage,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=SP_Percentage,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=SP_Percentage,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=SP_Percentage,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Smooth Cordgrass (%)")+
scale_y_continuous(limits = c(0,130), expand = c(0,0))+
  scale_x_continuous(expand=c(0,0))

sp22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.margin = margin(t=5, r=8, b=5))

plot(sp22.plot) 

#DS_Percentage
range(siteCovs(umfc.l)$DS_Percentage) #0-100

hplot <- data.frame(DS_Percentage=seq(min(covs[,c("DS_Percentage")]), max(covs[,c("DS_Percentage")]), 
                                      length=262), Fresh=0, Height=0, SP_Percentage=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm28z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=DS_Percentage,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=DS_Percentage,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=DS_Percentage,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=DS_Percentage,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Black Needlerush (%)")+
scale_y_continuous(limits = c(0,130), expand = c(0,0))+
  scale_x_continuous(expand=c(0,0))

ds22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.margin = margin(t=5, r=8, b=5))

plot(ds22.plot) 



