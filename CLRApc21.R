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

#Scale site covs
covs[,"EH"] <- scale(covs[,"EH"])
covs[,"Updev"] <- scale(covs[,"Updev"])
covs[,"Estnegdist"] <- scale(covs[,"Estnegdist"])
covs[,"Watdist"] <- scale(covs[,"Watdist"])
covs[,"Palfor"] <- scale(covs[,"Palfor"])
##===================================================================================
##Create unmarked df-----------------------------------------------------------------

umfc <- unmarkedFramePCount(y=CLRAc, obsCovs = list(date=(daymat), wind=(windmat), visit=(visit), obs=(obs)), 
                            siteCovs = covs[,c("EH", "Watdist", "Estnegdist", "Updev", "Palfor", "Point_ID")])

summary(umfc)
summary(apply(CLRAc,1,max,na.rm=TRUE))
#mean counts per point= 3.282
#max counts per point= 20

#Null model 
nullm <- pcount(~1~1, data = umfc, engine = "TMB")
backTransform(nullm, type="det") #0.275

#fit global model 
fmc<- pcount(~date+wind+obs+(1|Point_ID)+(1|visit)
             ~EH+Watdist+Estnegdist+Updev+Palfor+(1|Point_ID), umfc, engine = "TMB")

summary(fmc)



#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.626    0.791

#Fixed effects:
#  Estimate     SE      z   P(>|z|)
#(Intercept)   2.7658 0.1126 24.562 3.26e-133
#EH           -0.4482 0.0813 -5.509  3.60e-08
#Watdist      -0.4507 0.0945 -4.767  1.87e-06
#Estnegdist   -0.2559 0.0838 -3.055  2.25e-03
#Updev         0.2414 0.0825  2.927  3.42e-03
#Palfor       -0.0345 0.0819 -0.421  6.74e-01

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.047    0.217

#Fixed effects:
#  Estimate     SE      z  P(>|z|)
#(Intercept)   -2.299 0.1725 -13.33 1.55e-40
#date          -0.191 0.0543  -3.53 4.22e-04
#wind1         -0.263 0.0872  -3.02 2.57e-03
#wind2         -0.217 0.0963  -2.26 2.41e-02
#wind3         -0.411 0.1322  -3.11 1.87e-03
#wind4         -0.496 0.1425  -3.48 4.96e-04
#wind5         -0.866 0.3041  -2.85 4.39e-03

#AIC: 2821.781 

##===============================================================================================
##Model Building---------------------------------------------------------------------------------
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
#d15    11 3241.86   0.00 8.7e-01     0.87
#d13     6 3245.69   3.83 1.3e-01     1.00

##add abundance covs to best detection sub model 
cm1<- pcount(~date+wind+(1|visit)+obs ~EH, data=umfc, mixture= "P", engine = "TMB")
cm2<- pcount(~date+wind+(1|visit)+obs ~Watdist, data=umfc, mixture= "P", engine = "TMB")
cm3<- pcount(~date+wind+(1|visit)+obs ~Estnegdist, data=umfc, mixture= "P", engine = "TMB")
cm4<- pcount(~date+wind+(1|visit)+obs ~Updev, data=umfc, mixture= "P", engine = "TMB")
cm5<- pcount(~date+wind+(1|visit)+obs ~Palfor, data=umfc, mixture= "P", engine = "TMB")
cm6<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist, data=umfc, mixture= "P", engine = "TMB")
cm7<- pcount(~date+wind+(1|visit)+obs ~EH + Estnegdist, data=umfc, mixture= "P", engine = "TMB")
cm8<- pcount(~date+wind+(1|visit)+obs ~EH + Updev, data=umfc, mixture= "P", engine = "TMB")
cm9<- pcount(~date+wind+(1|visit)+obs ~EH + Palfor,  data=umfc, mixture= "P", engine = "TMB")
cm10<- pcount(~date+wind+(1|visit)+obs ~Watdist + Estnegdist, data=umfc, mixture= "P", engine = "TMB")
cm11<- pcount(~date+wind+(1|visit)+obs ~Watdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm12<- pcount(~date+wind+(1|visit)+obs ~Watdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm13<- pcount(~date+wind+(1|visit)+obs ~Estnegdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm14<- pcount(~date+wind+(1|visit)+obs ~Estnegdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm15<- pcount(~date+wind+(1|visit)+obs ~Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm16<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Estnegdist, data=umfc, mixture= "P", engine = "TMB")
cm17<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm18<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm19<- pcount(~date+wind+(1|visit)+obs ~EH + Estnegdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm20<- pcount(~date+wind+(1|visit)+obs ~EH + Estnegdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm21<- pcount(~date+wind+(1|visit)+obs ~EH + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm22<- pcount(~date+wind+(1|visit)+obs ~Watdist + Estnegdist + Updev,  data=umfc, mixture= "P", engine = "TMB")
cm23<- pcount(~date+wind+(1|visit)+obs ~Watdist + Estnegdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm24<- pcount(~date+wind+(1|visit)+obs ~Watdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm25<- pcount(~date+wind+(1|visit)+obs ~Estnegdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm26<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Estnegdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm27<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Estnegdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm28<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm29<- pcount(~date+wind+(1|visit)+obs ~EH + Estnegdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm30<- pcount(~date+wind+(1|visit)+obs ~Watdist + Estnegdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm31<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Estnegdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm32<- pcount(~date+wind+(1|visit)+obs ~(1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm33<- pcount(~date+wind+(1|visit)+obs ~EH + (1|Point_ID),  data=umfc, mixture= "P", engine = "TMB")
cm34<- pcount(~date+wind+(1|visit)+obs ~Watdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm35<- pcount(~date+wind+(1|visit)+obs ~Estnegdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm36<- pcount(~date+wind+(1|visit)+obs ~Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm37<- pcount(~date+wind+(1|visit)+obs ~Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm38<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm39<- pcount(~date+wind+(1|visit)+obs ~EH + Estnegdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm40<- pcount(~date+wind+(1|visit)+obs ~EH + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm41<- pcount(~date+wind+(1|visit)+obs ~EH + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm42<- pcount(~date+wind+(1|visit)+obs ~Watdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm43<- pcount(~date+wind+(1|visit)+obs ~Watdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm44<- pcount(~date+wind+(1|visit)+obs ~Watdist + Estnegdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm45<- pcount(~date+wind+(1|visit)+obs ~Estnegdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm46<- pcount(~date+wind+(1|visit)+obs ~Estnegdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm47<- pcount(~date+wind+(1|visit)+obs ~Palfor + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm48<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Estnegdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm49<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm50<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm51<- pcount(~date+wind+(1|visit)+obs ~EH + Estnegdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm52<- pcount(~date+wind+(1|visit)+obs ~EH + Estnegdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm53<- pcount(~date+wind+(1|visit)+obs ~EH + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm54<- pcount(~date+wind+(1|visit)+obs ~Watdist + Estnegdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm55<- pcount(~date+wind+(1|visit)+obs ~Watdist + Estnegdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm56<- pcount(~date+wind+(1|visit)+obs ~EH + Estnegdist + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm57<- pcount(~date+wind+(1|visit)+obs ~Watdist + Estnegdist + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm58<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Estnegdist + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")

#Create fit list of all models
cmL <- modSel(fitList(cm1, cm2, cm3, cm4, cm5, cm6, cm7, cm8, cm9, cm10, 
                      cm11, cm12, cm13, cm14, cm15, cm16, cm17, cm18, cm19, cm20,
                      cm21, cm22, cm23, cm24, cm25, cm26, cm27, cm28, cm29, cm30, cm31, cm32, cm33, cm34,
                      cm35, cm36,cm37, cm38, cm39, cm40, cm41, cm42, cm43, cm45, cm46, cm47, 
                      cm48, cm49, cm50, cm51, cm52, cm53, cm54, cm55, cm56, cm57, cm58))
cmL

#cm58    16 2803.38   0.00 8.5e-01     0.85
#cm48    14 2808.27   4.89 7.4e-02     0.92

##Best Models---------------------------------------------------------------------

summary(cm58)
#Abundance (log-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.554    0.744

#Fixed effects:
#  Estimate     SE     z   P(>|z|)
#(Intercept)   2.9011 0.1265 22.94 1.83e-116
#EH           -0.4584 0.0803 -5.71  1.16e-08
#Watdist      -0.4170 0.0907 -4.60  4.26e-06
#Estnegdist   -0.2405 0.0810 -2.97  2.98e-03
#Updev         0.1783 0.0799  2.23  2.56e-02
#Palfor        0.0153 0.0805  0.19  8.49e-01        non-sig 

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)     0.05    0.223

#Fixed effects:
#  Estimate     SE      z  P(>|z|)
#(Intercept)          -3.320 0.6244 -5.317 1.05e-07
#date                 -0.222 0.0549 -4.050 5.11e-05
#wind1                -0.193 0.0888 -2.179 2.93e-02
#wind2                -0.177 0.0981 -1.802 7.15e-02
#wind3                -0.381 0.1338 -2.848 4.39e-03
#wind4                -0.495 0.1444 -3.427 6.11e-04
#wind5                -0.895 0.3037 -2.946 3.22e-03
#obsMattDunning        1.057 0.6157  1.717 8.60e-02
#obsRachel Anderson    1.053 0.6220  1.693 9.04e-02
#obsSofia Campuzano    0.457 0.6245  0.731 4.65e-01

#AIC: 2803.379

#take out non-sig 
cm58<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Estnegdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")

summary(cm58)
#
#Abundance (log-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.554    0.744

#Fixed effects:
#  Estimate     SE     z   P(>|z|)
#(Intercept)    2.901 0.1260 23.02 2.81e-117
#EH            -0.456 0.0793 -5.75  8.93e-09
#Watdist       -0.420 0.0895 -4.69  2.75e-06
#Estnegdist    -0.243 0.0801 -3.03  2.46e-03
#Updev          0.188 0.0628  2.99  2.79e-03

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)     0.05    0.223

#Fixed effects:
#  Estimate     SE      z  P(>|z|)
#(Intercept)          -3.320 0.6240 -5.320 1.04e-07
#date                 -0.221 0.0542 -4.075 4.59e-05
#wind1                -0.194 0.0888 -2.180 2.92e-02
#wind2                -0.175 0.0975 -1.793 7.30e-02
#wind3                -0.380 0.1335 -2.842 4.48e-03
#wind4                -0.494 0.1443 -3.422 6.22e-04
#wind5                -0.892 0.3035 -2.941 3.27e-03
#obsMattDunning        1.054 0.6151  1.713 8.67e-02
#obsRachel Anderson    1.054 0.6217  1.696 9.00e-02
#obsSofia Campuzano    0.459 0.6241  0.735 4.62e-01

#AIC: 2801.415 

#ZIP without nonsig 
cm58Z<- pcount(~date+wind+(1|visit)+obs ~EH+Watdist+Estnegdist+Updev+(1|Point_ID),
               umfc, engine = "TMB", mixture="ZIP")

summary(cm58Z) #AIC: 2784.232 
#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.309    0.556

#Fixed effects:
#  Estimate     SE     z   P(>|z|)
#(Intercept)   3.1864 0.1412 22.57 8.03e-113
#EH           -0.3908 0.0819 -4.77  1.84e-06
#Watdist      -0.4120 0.0818 -5.04  4.72e-07
#Estnegdist   -0.2872 0.0745 -3.86  1.15e-04
#Updev         0.0978 0.0537  1.82  6.87e-02    #not sig 

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.049    0.221

#Fixed effects:
#  Estimate     SE      z  P(>|z|)
#(Intercept)          -3.322 0.6190 -5.368 7.97e-08
#date                 -0.197 0.0536 -3.680 2.33e-04
#wind1                -0.195 0.0881 -2.216 2.67e-02
#wind2                -0.181 0.0955 -1.891 5.86e-02
#wind3                -0.404 0.1322 -3.053 2.27e-03
#wind4                -0.491 0.1424 -3.445 5.70e-04
#wind5                -0.895 0.2984 -3.000 2.70e-03
#obsMattDunning        1.032 0.6109  1.690 9.10e-02
#obsRachel Anderson    0.982 0.6167  1.592 1.11e-01
#obsSofia Campuzano    0.366 0.6175  0.592 5.54e-01

#Zero-inflation (logit-scale):
#  Estimate    SE     z  P(>|z|)
#-2.12 0.363 -5.84 5.12e-09

#AIC: 2784.232


#ZIP is best! 




##=====================================================================================================
##Abundance prediction --------------------------------------------------------------------------------

#Poisson 
renb <- ranef(cm58)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#6775.594 6636.000 6920.000 


#ZIP 
renb <- ranef(cm58Z)
plot(renb, layout=c(4,3), subset=site%in%1:12, xlim=c(-1, 11), lwd=5)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#6682.357 6531.000 6843.025 



##GoF===============================================================================
##Plots-----------------------------------------------------------------------------

#Residuals
plot(cm58) 
plot(cm58Z) #CLRA.pc21.zresid


#expected vs observed counts 
n.obs <- apply(CLRAc, 1, sum, na.rm=TRUE)
n.predP <- apply(fitted(cm58), 1, sum, na.rm=TRUE )
n.predZ <- apply(fitted(cm58Z), 1, sum, na.rm=TRUE )

plot(n.obs, n.predP, frame=F, ylab="Expected Counts", xlab="Observed Counts", main="CLRA pc21 + random effect")
abline(0,1)
points(smooth.spline(n.predP ~ n.obs, df=4), type="l", lwd=2, col="blue")
points(smooth.spline(n.predZ ~ n.obs, df=4), type="l", lwd=2, col="green") 
#CLRA.pc21.ove

#P is better than Z 

#fitted vs residuals 
n.residsP <-apply(residuals(cm58),1,sum,na.rm=T)
plot(n.predP, n.residsP, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="CLRA pc21 + random effect")
points(smooth.spline(n.residsP ~ n.predP, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predP, df=4), type="l", lwd=2, col="black")
#CLRApc.nb.fit

#fitted vs residuals 
n.residsZ <-apply(residuals(cm58Z),1,sum,na.rm=T)
plot(n.predZ, n.residsZ, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="CLRA pc21 + random effect")
points(smooth.spline(n.residsZ ~ n.predZ, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predZ, df=4), type="l", lwd=2, col="black")
#CLRA.pc21.fit

##Parametric Bootstrap---------------------------------------------------------------------------- 
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



gof.p <- parboot(cm58, fitstats, nsim= 1000)
gof.p

#Parametric Bootstrap Statistics:
#Parametric Bootstrap Statistics:
#t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
#SSE          1836          733.0            60.99            0
#Chisq         760          199.4            21.85            0
#freemanTukey  363           90.4             9.46            0



#compute c-hat 
gof.p@t0[2] / mean(gof.p@t.star[,2])
#Chisq 
#1.355801 
 

gofZ <- parboot(fmcZ2, fitstats, nsim= 1000)
gofZ

#Parametric Bootstrap Statistics:
#t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
#SSE          1855          770.9             96.9        0.002
#Chisq         769          220.4             26.0        0.000
#freemanTukey  358           91.8             10.9        0.000

#t_B quantiles:
#  0% 2.5%  25%  50%  75% 97.5% 100%
#SSE          864  963 1031 1072 1126  1256 2755
#Chisq        488  507  531  546  563   604  709
#freemanTukey 241  249  259  265  272   291  341


##Nmixgof-------------------------------------------------------------------------------

#P
#chat 
chat(fmc2, type = "marginal")
#0.9963905

#Residuals against covariates
residcov(fmc2)

#Residual against fitted 
residfit(fmc2, type="marginal")
residfit(fmc2, type="site-sum")
residfit(fmc2, type="observation")

#Qq plot of randomized quantile residuals against standard normal quantiles
CLRA21s.qq <- residqq(cm58Z, "site-sum")
CLRA21obs.qq <- residqq(cm58Z, "observation") #CLRA21pc.rand_qq.obs


#three types of randomized quantile residuals for binomial N-mixture models
rqresiduals(fmc2, type = "marginal")
rqresiduals(fmc2, type="site-sum")



#ZIP 
#Qq plot of randomized quantile residuals against standard normal quantiles
residqq(cm58Z, "site-sum") #CLRA.pc21.ss.qq
residqq(cm58Z, "observation")#CLRA.pc21.ob.qq





##=============================================================================
##Cov plots 

#Take out random effects
cm58Z<- pcount(~date+wind+obs ~EH+Watdist+Estnegdist+Updev,
               umfc, engine = "TMB", mixture="ZIP")


#EH
range(siteCovs(umfc)$EH) 

ehplot <- data.frame(EH=seq(min(covs[,c("EH")]), max(covs[,c("EH")]), 
                           length=262),Watdist=0, Estnegdist=0, Updev=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

eh.pred <- predict(cm58Z, type="state", newdata=ehplot, appendData=TRUE)

eh.plot <- data.frame(eh.pred,ehplot)

cwwi.plot <- ggplot(data=eh.plot,aes(x=EH,y=Predicted,ymin=lower,ymax=upper))

cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=EH,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=EH,y=upper), alpha=0.1) 

cwwi.plot <- cwwi.plot + geom_line(aes(x=EH,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Landscape Heterogeneity")+
  scale_y_continuous(limits = c(0,70), expand = c(0,0))+
  scale_x_continuous(breaks= c(-2, -1, 0, 1), labels = c("0.0", "0.2", "0.4", "0.6"), expand=c(0,0))

eh21.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(eh21.plot) 


#Watdist
range(siteCovs(umfc)$Watdist) 

hplot <- data.frame(Watdist=seq(min(covs[,c("Watdist")]), max(covs[,c("Watdist")])), length=262,
                    EH=0, Estnegdist=0, Updev=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm58Z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Watdist,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Watdist,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Watdist,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Watdist,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Distance to Open Water (km)")+
  scale_y_continuous(limits = c(0,70), expand = c(0,0))+
  scale_x_continuous(breaks= c(0, 2.5, 5), labels = c("0", "3", "6"),expand = c(0,0))

w21.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(w21.plot) 


#Estnegdist
range(siteCovs(umfc)$Estnegdist) 

hplot <- data.frame(Estnegdist=seq(min(covs[,c("Estnegdist")]), max(covs[,c("Estnegdist")])), length=262,
                    EH=0, Watdist=0, Updev=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm58Z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Estnegdist,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Estnegdist,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Estnegdist,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Estnegdist,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Distance to Estuarine Marsh (km)")+
  scale_y_continuous(limits = c(0,70), expand = c(0,0))+
  scale_x_continuous(breaks= c(-2, 0, 2), labels = c("-1", "0", "1"), expand = c(0,0))

e21.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(e21.plot) 



