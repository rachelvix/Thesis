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

setwd("C:/Users/lab/Desktop/Rachel/Analysis")


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
#SESP counts and covariates in wide format 
#262 point count sites in total
#3 visits in year 2021 & 2022 
SESP.w <- read.table("pc.data.csv", header=T, sep=",")   
head(SESP.w)
str(SESP.w)

#select the first year 
covs <- SESP.w[1:262,]

#Counts with truncation and passive listening only 
SESPc <- as.matrix(covs[,c("SESP1", "SESP2", "SESP3")])

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
covs$nobs.1 <- as.factor(covs$nobs.1)
covs$nobs.2 <- as.factor(covs$nobs.2)
covs$nobs.3 <- as.factor(covs$nobs.3)
obs <- as.matrix(cbind(covs[,c("nobs.1", "nobs.2", "nobs.3")]))

##===========================================================================================
##Model with just local variables 

#Scale local veg covs
covs[,"SP_Percentage"] <- scale(covs[,"SP_Percentage"])
covs[,"DS_Percentage"] <- scale(covs[,"DS_Percentage"])
covs[,"Height"] <- scale(covs[,"Height"])
covs[,"Fresh"] <- scale(covs[,"Fresh"])

##==========================================================================================
##Create unmarked df------------------------------------------------------------------------ 

umc.l <- unmarkedFramePCount(y=SESPc, obsCovs = list(date=(daymat), wind=(windmat), visit=(visit), obs=(obs)), 
                              siteCovs = covs[,c("Point_ID", "DS_Percentage", "SP_Percentage", "Height", "Fresh")])
summary(umc.l)
summary(apply(SESPc,1,max,na.rm=TRUE))
#mean counts per point= 1.095
#max counts per point= 15

#Global model 
fmc.l<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage+SP_Percentage+Height+(1|Point_ID), umc.l)
summary(fmc.l)

##=================================================================================
##Model Building ------------------------------------------------------------------
#Detection covs first  
mlist = list()

mlist$d1 = pcount(~date ~1, umc.l, engine = "TMB")
mlist$d2 = pcount(~wind ~1, umc.l, engine = "TMB")
mlist$d3 = pcount(~(1|visit) ~1, umc.l, engine = "TMB")
mlist$d4 = pcount(~date+wind ~1, umc.l, engine = "TMB")
mlist$d5 = pcount(~date+(1|visit) ~1, umc.l, engine = "TMB")
mlist$d6 = pcount(~wind+(1|visit) ~1, umc.l, engine = "TMB")
mlist$d7 = pcount(~date+wind+(1|visit) ~1, umc.l, engine = "TMB")
mlist$d8 = pcount(~obs ~1, umc.l, engine = "TMB")
mlist$d9 = pcount(~date+obs ~1, umc.l, engine = "TMB")
mlist$d10 = pcount(~(1|visit)+obs ~1, umc.l, engine = "TMB")
mlist$d11 = pcount(~wind+obs ~1, umc.l, engine = "TMB")
mlist$d12 = pcount(~date+wind+obs ~1, umc.l, engine = "TMB")
mlist$d13 = pcount(~date+(1|visit)+obs ~1, umc.l, engine = "TMB")
mlist$d14 = pcount(~wind+(1|visit)+obs ~1, umc.l, engine = "TMB")
mlist$d15 = pcount(~date+wind+(1|visit)+obs ~1, umc.l, engine = "TMB")


fl1 <- modSel(fitList(fits = mlist)) # Compare models. Move forward with best model(s) (delta AIC < 2)
#nPars     AIC delta   AICwt cumltvWt
#d12    11 1725.54  0.00 6.9e-01     0.69
#d15    11 1727.54  2.00 2.6e-01     0.95

##add abundance covs to best model 
cm1<- pcount(~date+wind+obs ~SP_Percentage, data=umc.l, mixture= "P", engine = "TMB")
cm2<- pcount(~date+wind+obs ~DS_Percentage, data=umc.l, mixture= "P", engine = "TMB")
cm3<- pcount(~date+wind+obs ~Height, data=umc.l, mixture= "P", engine = "TMB")
cm4<- pcount(~date+wind+obs ~SP_Percentage + DS_Percentage, data=umc.l, mixture= "P", engine = "TMB")
cm5<- pcount(~date+wind+obs ~SP_Percentage + Height, data=umc.l, mixture= "P", engine = "TMB")
cm6<- pcount(~date+wind+obs ~DS_Percentage + Height, data=umc.l, mixture= "P", engine = "TMB")
cm7<- pcount(~date+wind+obs ~SP_Percentage + DS_Percentage + Height, data=umc.l, mixture= "P", engine = "TMB")
cm8<- pcount(~date+wind+obs ~(1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm9<- pcount(~date+wind+obs ~SP_Percentage + (1|Point_ID),  data=umc.l, mixture= "P", engine = "TMB")
cm10<- pcount(~date+wind+obs ~DS_Percentage + (1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm11<- pcount(~date+wind+obs ~Height + (1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm12<- pcount(~date+wind+obs ~SP_Percentage + DS_Percentage + (1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm13<- pcount(~date+wind+obs ~SP_Percentage + Height + (1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm14<- pcount(~date+wind+obs ~DS_Percentage + Height + (1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm15<- pcount(~date+wind+obs ~SP_Percentage + DS_Percentage + Height + (1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm16<- pcount(~date+wind+obs ~Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm17<- pcount(~date+wind+obs ~SP_Percentage + Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm18<- pcount(~date+wind+obs ~DS_Percentage + Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm19<- pcount(~date+wind+obs ~Height + Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm20<- pcount(~date+wind+obs ~SP_Percentage + DS_Percentage + Height +Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm21<- pcount(~date+wind+obs ~(1|Point_ID)+Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm22<- pcount(~date+wind+obs ~SP_Percentage + (1|Point_ID) +Fresh,  data=umc.l, mixture= "P", engine = "TMB")
cm23<- pcount(~date+wind+obs ~DS_Percentage + (1|Point_ID) +Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm24<- pcount(~date+wind+obs ~Height + (1|Point_ID) +Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm25<- pcount(~date+wind+obs ~SP_Percentage + DS_Percentage + (1|Point_ID)+Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm26<- pcount(~date+wind+obs ~SP_Percentage + Height + (1|Point_ID) +Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm27<- pcount(~date+wind+obs ~DS_Percentage + Height + (1|Point_ID) +Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm28<- pcount(~date+wind+obs ~SP_Percentage + DS_Percentage + Height + (1|Point_ID) +Fresh, data=umc.l, mixture= "P", engine = "TMB")

#Create fit list of all models
cmL <- modSel(fitList(cm1, cm2, cm3, cm4, cm5, cm6, cm7, cm8, cm9, cm10, 
                      cm11, cm12, cm13, cm14, cm15, cm16, cm17, cm18, cm19, cm20, cm21, cm22, cm23, cm24, cm25, cm26, cm27, cm28))

cmL

#nPars     AIC  delta    AICwt cumltvWt
#cm24    13 1178.73   0.00  3.0e-01     0.30
#cm26    14 1178.94   0.20  2.7e-01     0.57
#cm28    15 1179.15   0.41  2.4e-01     0.81
#cm27    14 1179.69   0.96  1.9e-01     1.00
#cm22    13 1195.74  17.01  6.1e-05     1.00

#Best Models---------------------------------------------------------------------------------

summary(cm26)
summary(cm24)


#ZIP
cm24z<- pcount(~date+wind+obs ~Height + (1|Point_ID) +Fresh, data=umc.l,
               mixture= "ZIP", engine = "TMB")
summary(cm24z)

#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.913    0.956

#Fixed effects:
#  Estimate    SE     z  P(>|z|)
#(Intercept)    1.118 0.392  2.85 4.31e-03
#Height        -0.446 0.186 -2.40 1.65e-02
#Fresh         -1.929 0.381 -5.07 4.02e-07

#Detection (logit-scale):
#  Estimate    SE       z  P(>|z|)
#(Intercept)          -2.18627 0.301 -7.2555 4.00e-13
#date                  0.08126 0.103  0.7873 4.31e-01
#wind1                 0.01480 0.188  0.0787 9.37e-01
#wind2                -0.00358 0.163 -0.0220 9.82e-01
#wind3                 0.27254 0.200  1.3608 1.74e-01
#wind4                -0.17356 0.221 -0.7862 4.32e-01
#wind5                -0.54640 0.345 -1.5859 1.13e-01
#obsAllister Rutledge  0.45034 0.388  1.1605 2.46e-01
#obsBrittany Holliker  0.61217 0.314  1.9484 5.14e-02
#obsMattDunning        0.82609 0.733  1.1272 2.60e-01

#Zero-inflation (logit-scale):
#  Estimate    SE    z P(>|z|)
#0.193 0.191 1.01   0.313

#AIC: 1173.239 

##ZIP is best! 



##====================================================================================================
##Abundance prediction-----------------------------------------------------------------------------

#Poisson 
renb <- ranef(cm24)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#1259.423 1198.975 1313.000  

 

##Z
renb <- ranef(cm24z)
plot(renb, layout=c(4,3), subset=site%in%1:12, xlim=c(-1, 11), lwd=5)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
##Estimate     2.5%    97.5% 
# 881.312  832.000  927.000 



##=========================================================================
##Gof----------------------------------------------------------------------

#Plots 
#Residuals
plot(cm24) #saved as SESPpc.local.p.resid
plot(cm24z) #saved as SESPpc.local.zip.resid


#expected vs observed counts 
n.obs <- apply(SESPc, 1, sum, na.rm=TRUE)
n.predP <- apply(fitted(cm24), 1, sum, na.rm=TRUE )
n.predZ <- apply(fitted(cm24z), 1, sum, na.rm=TRUE )
plot(n.obs, n.predP, frame=F, ylab="Expected Counts", xlab="Observed Counts", main="SESP local pcount")
abline(0,1)
points(smooth.spline(n.predP ~ n.obs, df=4), type="l", lwd=2, col="blue")
points(smooth.spline(n.predZ ~ n.obs, df=4), type="l", lwd=2, col="red")

#Z looks worse! 

#fitted vs residuals 
n.residsZ <-apply(residuals(cm24z),1,sum,na.rm=T)
plot(n.predZ, n.residsZ, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="SESP stacked pcount NB")
points(smooth.spline(n.residsZ ~ n.predZ, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predZ, df=4), type="l", lwd=2, col="black")
#SESPpc.local.fresid

#fitted vs residuals 
n.residsP <-apply(residuals(cm24),1,sum,na.rm=T)
plot(n.predP, n.residsP, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="SESP stacked pcount NB")
points(smooth.spline(n.residsP ~ n.predP, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predP, df=4), type="l", lwd=2, col="black")
#SESPpc.local.fresid

##Nmixgof----------------------------------------------------------------------------------
#Chat 
chat(cm24, type="marginal") #0.9515996
chat(cm24z, type="marginal") #0.9398458


#residuals 
residqq(cm24, type="site-sum")
residqq(cm24z, type="site-sum")

residqq(cm24, type="observation")
residqq(cm24z, type="observation")



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
gofp.l <- parboot(cm24, fitstats, nsim= 1000)
gofp.l
#  t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
#SSE          1841          791.5            63.88            0
#Chisq         736          209.4            20.29            0
#freemanTukey  361           99.2             8.61            0


##===============================================================================
##Plot Covs 

#take out random effects
#ZIP
cm24z<- pcount(~date+wind+obs ~Height +Fresh, data=umc.l,
               mixture= "ZIP", engine = "TMB")


#Height 
range(siteCovs(umc.l)$Height) #~39-188

hplot <- data.frame(Height=seq(min(covs[,c("Height")]), max(covs[,c("Height")]), 
                               length=262), DS_Percentage=0, SP_Percentage=0, Fresh=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm24z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Height,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Height,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Height,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Height,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Vegetation Height (cm)")+
  scale_y_continuous(limits = c(0,10), expand = c(0,0))+
  scale_x_continuous(expand=c(0,0))

h21.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(h21.plot) 


#Fresh
range(siteCovs(umc.l)$Fresh) #0-100

hplot <- data.frame(Fresh=seq(min(covs[,c("Fresh")]), max(covs[,c("Fresh")]), 
                              length=262), DS_Percentage=0, Height=0, SP_Percentage=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm24z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Fresh,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Fresh,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Fresh,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Fresh,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Freshwater Vegetation (%)")+
  scale_y_continuous(limits = c(0,10), expand = c(0,0), breaks= c(0, 2.5, 5.0, 7.5, 10),
                     labels = c("0", "3", "5", "8", "10"))+
  scale_x_continuous(expand=c(0,0))

f21.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.margin = margin(t=5, r=10, b=5))

plot(f21.plot) 












##===============================================================================
#select year 2022
covs <- SESP.w[263:524,]

#Counts with truncation and passive listening only 
SESPc <- as.matrix(covs[,c("SESP1", "SESP2", "SESP3")])

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
covs$nobs.1 <- as.factor(covs$nobs.1)
covs$nobs.2 <- as.factor(covs$nobs.2)
covs$nobs.3 <- as.factor(covs$nobs.3)
obs <- as.matrix(cbind(covs[,c("nobs.1", "nobs.2", "nobs.3")]))

##===========================================================================================
##Model with just local variables 

#Scale local veg covs
covs[,"SP_Percentage"] <- scale(covs[,"SP_Percentage"])
covs[,"DS_Percentage"] <- scale(covs[,"DS_Percentage"])
covs[,"Height"] <- scale(covs[,"Height"])
covs[,"Fresh"] <- scale(covs[,"Fresh"])

##==========================================================================================
##Create unmarked df------------------------------------------------------------------------ 

umc.l <- unmarkedFramePCount(y=SESPc, obsCovs = list(date=(daymat), wind=(windmat), visit=(visit), obs=(obs)), 
                             siteCovs = covs[,c("Point_ID", "DS_Percentage", "SP_Percentage", "Height", "Fresh")])
summary(umc.l)
summary(apply(SESPc,1,max,na.rm=TRUE))
#mean counts per point= 1.095
#max counts per point= 15

#Global model 
fmc.l<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage+SP_Percentage+Height+(1|Point_ID), umc.l)
summary(fmc.l)

##=================================================================================
##Model Building ------------------------------------------------------------------
#Detection covs first  
mlist = list()

mlist$d1 = pcount(~date ~1, umc.l, engine = "TMB")
mlist$d2 = pcount(~wind ~1, umc.l, engine = "TMB")
mlist$d3 = pcount(~(1|visit) ~1, umc.l, engine = "TMB")
mlist$d4 = pcount(~date+wind ~1, umc.l, engine = "TMB")
mlist$d5 = pcount(~date+(1|visit) ~1, umc.l, engine = "TMB")
mlist$d6 = pcount(~wind+(1|visit) ~1, umc.l, engine = "TMB")
mlist$d7 = pcount(~date+wind+(1|visit) ~1, umc.l, engine = "TMB")
mlist$d8 = pcount(~obs ~1, umc.l, engine = "TMB")
mlist$d9 = pcount(~date+obs ~1, umc.l, engine = "TMB")
mlist$d10 = pcount(~(1|visit)+obs ~1, umc.l, engine = "TMB")
mlist$d11 = pcount(~wind+obs ~1, umc.l, engine = "TMB")
mlist$d12 = pcount(~date+wind+obs ~1, umc.l, engine = "TMB")
mlist$d13 = pcount(~date+(1|visit)+obs ~1, umc.l, engine = "TMB")
mlist$d14 = pcount(~wind+(1|visit)+obs ~1, umc.l, engine = "TMB")
mlist$d15 = pcount(~date+wind+(1|visit)+obs ~1, umc.l, engine = "TMB")


fl1 <- modSel(fitList(fits = mlist)) # Compare models. Move forward with best model(s) (delta AIC < 2)
#nPars     AIC delta   AICwt cumltvWt
#d15    13 1081.40   0.00 6.8e-01     0.68
#d12    13 1082.90   1.49 3.2e-01     1.00

##add abundance covs to best model 
cm1<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage, data=umc.l, mixture= "P", engine = "TMB")
cm2<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage, data=umc.l, mixture= "P", engine = "TMB")
cm3<- pcount(~date+wind+(1|visit)+obs ~Height, data=umc.l, mixture= "P", engine = "TMB")
cm4<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage, data=umc.l, mixture= "P", engine = "TMB")
cm5<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Height, data=umc.l, mixture= "P", engine = "TMB")
cm6<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + Height, data=umc.l, mixture= "P", engine = "TMB")
cm7<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height, data=umc.l, mixture= "P", engine = "TMB")
cm8<- pcount(~date+wind+(1|visit)+obs ~(1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm9<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + (1|Point_ID),  data=umc.l, mixture= "P", engine = "TMB")
cm10<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + (1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm11<- pcount(~date+wind+(1|visit)+obs ~Height + (1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm12<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + (1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm13<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Height + (1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm14<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + Height + (1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm15<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height + (1|Point_ID), data=umc.l, mixture= "P", engine = "TMB")
cm16<- pcount(~date+wind+(1|visit)+obs ~Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm17<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm18<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm19<- pcount(~date+wind+(1|visit)+obs ~Height + Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm20<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height +Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm21<- pcount(~date+wind+(1|visit)+obs ~(1|Point_ID)+Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm22<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + (1|Point_ID) +Fresh,  data=umc.l, mixture= "P", engine = "TMB")
cm23<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + (1|Point_ID) +Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm24<- pcount(~date+wind+(1|visit)+obs ~Height + (1|Point_ID) +Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm25<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + (1|Point_ID)+Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm26<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + Height + (1|Point_ID) +Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm27<- pcount(~date+wind+(1|visit)+obs ~DS_Percentage + Height + (1|Point_ID) +Fresh, data=umc.l, mixture= "P", engine = "TMB")
cm28<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height + (1|Point_ID) +Fresh, data=umc.l, mixture= "P", engine = "TMB")

#Create fit list of all models
cmL <- modSel(fitList(cm1, cm2, cm3, cm4, cm5, cm6, cm7, cm8, cm9, cm10, 
                      cm11, cm12, cm13, cm14, cm15, cm16, cm17, cm18, cm19, cm20, cm21, cm22, cm23, cm24, cm25, cm26, cm27, cm28))

cmL

#nPars     AIC  delta    AICwt cumltvWt
#cm15    16  707.71   0.00 2.9e-01     0.29
#cm28    17  708.36   0.66 2.1e-01     0.49
#cm26    16  709.69   1.99 1.1e-01     0.60

#Best Models---------------------------------------------------------------------------------

summary(cm28)
summary(cm26)
summary(cm15)
#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    9.027    3.004

#Fixed effects:
#  Estimate    SE     z P(>|z|)
#(Intercept)     -1.801 0.647 -2.78 0.00540
#SP_Percentage    0.730 0.257  2.84 0.00445
#DS_Percentage    0.755 0.311  2.43 0.01518
#Height          -0.918 0.364 -2.52 0.01180

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)     0.02    0.142

#Fixed effects:
#  Estimate     SE       z  P(>|z|)
#(Intercept)           -4.891  0.391 -12.507 6.83e-36
#date                  -0.548  0.169  -3.250 1.16e-03
#wind1                  0.208  0.253   0.821 4.12e-01
#wind2                  0.569  0.282   2.017 4.37e-02
#wind3                  0.249  0.293   0.850 3.95e-01
#wind4                  0.733  0.412   1.779 7.53e-02
#wind5                 -7.053 64.249  -0.110 9.13e-01
#obsDanielle Canning    1.599  0.545   2.932 3.37e-03
#obsFred Hall           2.338  0.413   5.661 1.50e-08
#obsMiranda Mordue      1.720  0.433   3.975 7.05e-05
#obsSofia Campuzano    -3.001 27.012  -0.111 9.12e-01
#obsTyler Connell       0.692  1.151   0.601 5.48e-01

#AIC: 707.7052 

#ZIP
cm15z<- pcount(~date+wind+(1|visit)+obs ~SP_Percentage + DS_Percentage + Height + (1|Point_ID), data=umc.l,
               mixture= "ZIP", engine = "TMB")
summary(cm15z) #AIC: 724.8177 




##P is best! 



##====================================================================================================
##Abundance prediction-----------------------------------------------------------------------------

#Poisson 
renb <- ranef(cm15)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#1194.435 1133.000 1254.000 
  



##Z
renb <- ranef(cm15z)
plot(renb, layout=c(4,3), subset=site%in%1:12, xlim=c(-1, 11), lwd=5)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
##Estimate     2.5%    97.5% 
# 665.974  617.000  718.000  



##=========================================================================
##Gof----------------------------------------------------------------------

#Plots 
#Residuals
plot(cm15) #saved as SESPpc.local.p.resid
plot(cm15z) #saved as SESPpc.local.zip.resid


#expected vs observed counts 
n.obs <- apply(SESPc, 1, sum, na.rm=TRUE)
n.predP <- apply(fitted(cm15), 1, sum, na.rm=TRUE )
n.predZ <- apply(fitted(cm15z), 1, sum, na.rm=TRUE )
plot(n.obs, n.predP, frame=F, ylab="Expected Counts", xlab="Observed Counts", main="SESP local pcount")
abline(0,1)
points(smooth.spline(n.predP ~ n.obs, df=4), type="l", lwd=2, col="blue")
points(smooth.spline(n.predZ ~ n.obs, df=4), type="l", lwd=2, col="red")


#fitted vs residuals 
n.residsZ <-apply(residuals(cm15z),1,sum,na.rm=T)
plot(n.predZ, n.residsZ, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="SESP stacked pcount NB")
points(smooth.spline(n.residsZ ~ n.predZ, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predZ, df=4), type="l", lwd=2, col="black")
#SESPpc.local.fresid

#fitted vs residuals 
n.residsP <-apply(residuals(cm15),1,sum,na.rm=T)
plot(n.predP, n.residsP, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="SESP stacked pcount NB")
points(smooth.spline(n.residsP ~ n.predP, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predP, df=4), type="l", lwd=2, col="black")
#SESPpc.local.fresid

##Nmixgof----------------------------------------------------------------------------------
#Chat 
chat(cm15, type="marginal") #0.9515996
chat(cm15z, type="marginal") #0.9398458


#residuals 
residqq(cm15, type="site-sum")
residqq(cm15z, type="site-sum")

residqq(cm15, type="observation")
residqq(cm15z, type="observation")



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
gofp.l <- parboot(cm15, fitstats, nsim= 1000)
gofp.l
#  t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
#SSE          1841          791.5            63.88            0
#Chisq         736          209.4            20.29            0
#freemanTukey  361           99.2             8.61            0


##===============================================================================
##Plot Covs 

#take out random effects

cm15<- pcount(~date+wind+obs ~SP_Percentage + DS_Percentage + Height,
              data=umc.l, mixture= "P", engine = "TMB")



#Height 
range(siteCovs(umc.l)$Height) #~39-188

hplot <- data.frame(Height=seq(min(covs[,c("Height")]), max(covs[,c("Height")]), 
                               length=262), DS_Percentage=0, SP_Percentage=0, Fresh=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm15, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Height,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Height,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Height,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Height,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Vegetation Height (cm)")+
  scale_y_continuous(limits=c(0,4), breaks=seq(0,4,by=1), expand = c(0,0))+
  scale_x_continuous(expand=c(0,0))

h22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(h22.plot) 


#DS_Percentage
range(siteCovs(umc.l)$DS_Percentage) #0-100

hplot <- data.frame(DS_Percentage=seq(min(covs[,c("DS_Percentage")]), max(covs[,c("DS_Percentage")]), 
                              length=262), DS_Percentage=0, Height=0, SP_Percentage=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm15, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=DS_Percentage,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=DS_Percentage,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=DS_Percentage,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=DS_Percentage,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Black Needlerush (%)")+
  scale_y_continuous(limits = c(0,1250), expand = c(0,0))+
  scale_x_continuous(expand=c(0,0))

ds22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.margin = margin(t=5, r=8, b=5))

plot(ds22.plot) 


#SP_Percentage
range(siteCovs(umc.l)$SP_Percentage) #0-100

hplot <- data.frame(SP_Percentage=seq(min(covs[,c("SP_Percentage")]), max(covs[,c("SP_Percentage")]), 
                                      length=262), DS_Percentage=0, Height=0, SP_Percentage=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm15, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=SP_Percentage,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=SP_Percentage,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=SP_Percentage,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=SP_Percentage,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Smooth Cordgrass (%)")+
scale_y_continuous(limits = c(0,200), expand = c(0,0))+
  scale_x_continuous(expand=c(0,0))

sp22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.margin = margin(t=5, r=8, b=5))

plot(sp22.plot) 



