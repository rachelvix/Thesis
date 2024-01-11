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
fmc<- pcount(~date+wind+obs+(1|visit)
             ~EH+Watdist+Estnegdist+Updev+Palfor+(1|Point_ID), umfc, engine = "TMB")

summary(fmc)

#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)     0.72    0.848

#Fixed effects:
#  Estimate     SE       z   P(>|z|)
#(Intercept)  2.87575 0.1009 28.4880 1.65e-178
#EH          -0.44582 0.0876 -5.0892  3.60e-07
#Watdist     -0.53024 0.1025 -5.1756  2.27e-07
#Estnegdist  -0.44970 0.0920 -4.8906  1.01e-06
#Updev        0.38312 0.0885  4.3276  1.51e-05
#Palfor      -0.00836 0.0874 -0.0956  9.24e-01    ns 

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.006    0.075

#Fixed effects:
#  Estimate      SE      z P(>|z|)
#(Intercept)         -9.7915 35.7948 -0.274  0.7844
#date                -0.1253  0.0508 -2.464  0.0137
#wind1                0.0768  0.0832  0.924  0.3557
#wind2                0.0968  0.0948  1.021  0.3070
#wind3               -0.0457  0.1224 -0.373  0.7090
#wind4               -0.3407  0.1516 -2.248  0.0246
#wind5               -0.9113  0.3551 -2.566  0.0103
#obsJamie Russell     7.3003 35.7948  0.204  0.8384
#obsRachel Anderson   7.0266 35.7947  0.196  0.8444
#obsTyler Connell     6.9562 35.7948  0.194  0.8459

#AIC: 2963.42 

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
#d15    11 3574.80   0.00 6.0e-01     0.60
#d12    11 3575.65   0.85 3.9e-01     0.99

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

#cm58    16 2963.42   0.00  1.0e+00     1.00
#cm50    14 2984.60  21.18  2.5e-05     1.00

##Best Models---------------------------------------------------------------------

summary(cm58)
#Random effects:
#Groups        Name Variance Std.Dev.
#Point_ID (Intercept)     0.72    0.848

#Fixed effects:
#  Estimate     SE       z   P(>|z|)
#(Intercept)  2.87575 0.1009 28.4880 1.65e-178
#EH          -0.44582 0.0876 -5.0892  3.60e-07
#Watdist     -0.53024 0.1025 -5.1756  2.27e-07
#Estnegdist  -0.44970 0.0920 -4.8906  1.01e-06
#Updev        0.38312 0.0885  4.3276  1.51e-05
#Palfor      -0.00836 0.0874 -0.0956  9.24e-01   ns

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.006    0.075

#Fixed effects:
#  Estimate      SE      z P(>|z|)
#(Intercept)         -9.7915 35.7948 -0.274  0.7844
#date                -0.1253  0.0508 -2.464  0.0137
#wind1                0.0768  0.0832  0.924  0.3557
#wind2                0.0968  0.0948  1.021  0.3070
#wind3               -0.0457  0.1224 -0.373  0.7090
#wind4               -0.3407  0.1516 -2.248  0.0246
#wind5               -0.9113  0.3551 -2.566  0.0103
#obsJamie Russell     7.3003 35.7948  0.204  0.8384
#obsRachel Anderson   7.0266 35.7947  0.196  0.8444
#obsTyler Connell     6.9562 35.7948  0.194  0.8459

#AIC: 2963.42 

#take out non-sig 
cm58<- pcount(~date+wind+(1|visit)+obs 
              ~EH + Watdist + Estnegdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")

summary(cm58)
#
#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)     0.72    0.849

#Fixed effects:
#  Estimate     SE     z   P(>|z|)
#(Intercept)    2.876 0.1009 28.49 1.47e-178
#EH            -0.447 0.0869 -5.14  2.68e-07
#Watdist       -0.528 0.1008 -5.24  1.59e-07
#Estnegdist    -0.449 0.0910 -4.93  8.12e-07
#Updev          0.378 0.0701  5.39  7.02e-08

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.006    0.076

#Fixed effects:
#  Estimate      SE      z P(>|z|)
#(Intercept)         -8.2614 16.6895 -0.495  0.6206
#date                -0.1252  0.0509 -2.461  0.0139
#wind1                0.0762  0.0831  0.917  0.3589
#wind2                0.0962  0.0947  1.017  0.3093
#wind3               -0.0463  0.1222 -0.379  0.7049
#wind4               -0.3418  0.1514 -2.258  0.0240
#wind5               -0.9127  0.3550 -2.571  0.0102
#obsJamie Russell     5.7713 16.6894  0.346  0.7295
#obsRachel Anderson   5.4953 16.6892  0.329  0.7419
#obsTyler Connell     5.4258 16.6894  0.325  0.7451

#AIC: 2961.435  

#ZIP without nonsig 
cm58Z<- pcount(~date+wind+(1|visit)+obs 
               ~EH+Watdist+Estnegdist+Updev+(1|Point_ID),
               umfc, engine = "TMB", mixture="ZIP")

summary(cm58Z) #AIC: 2968.788 

#P is best 

cmL <- modSel(fitList(cm58, cm58Z)) 


##=====================================================================================================
##Abundance prediction --------------------------------------------------------------------------------

#Poisson 
renb <- ranef(cm58)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#7952.174 7800.950 8101.025 



#ZIP 
renb <- ranef(cm58Z)
plot(renb, layout=c(4,3), subset=site%in%1:12, xlim=c(-1, 11), lwd=5)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#7331.033 7164.000 7490.050 
 





##GoF===============================================================================
##Plots-----------------------------------------------------------------------------

#Residuals
plot(cm58) 
plot(cm58Z)


#expected vs observed counts 
n.obs <- apply(CLRAc, 1, sum, na.rm=TRUE)
n.predP <- apply(fitted(cm58), 1, sum, na.rm=TRUE )
n.predZ <- apply(fitted(cm58Z), 1, sum, na.rm=TRUE )

plot(n.obs, n.predP, frame=F, ylab="Expected Counts", xlab="Observed Counts", main="CLRA pc22")
abline(0,1)
points(smooth.spline(n.predP ~ n.obs, df=4), type="l", lwd=2, col="blue")
points(smooth.spline(n.predZ ~ n.obs, df=4), type="l", lwd=2, col="green") 
#CLRA.pc22.ove

#P is better than Z 

#fitted vs residuals 
n.residsP <-apply(residuals(cm58),1,sum,na.rm=T)
plot(n.predP, n.residsP, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="CLRA pc21 + random effect")
points(smooth.spline(n.residsP ~ n.predP, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predP, df=4), type="l", lwd=2, col="black")
#CLRA.pc22.fit

#fitted vs residuals 
n.residsZ <-apply(residuals(cm58Z),1,sum,na.rm=T)
plot(n.predZ, n.residsZ, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="CLRA pc21 + random effect")
points(smooth.spline(n.residsZ ~ n.predZ, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predZ, df=4), type="l", lwd=2, col="black")


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
#t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
#SSE          2587           1438            69.83            0
#Chisq         917            391            21.45            0
#freemanTukey  397            146             8.96            0

#t_B quantiles:
#  0% 2.5%  25%  50%  75% 97.5% 100%
#SSE          980 1028 1098 1145 1192  1296 1456
#Chisq        470  490  512  526  539   574  603
#freemanTukey 229  234  245  251  257   269  283


#compute c-hat 
gof.p@t0[2] / mean(gof.p@t.star[,2])
#Chisq 
#1.742535 
 


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
residqq(cm58Z, "site-sum") #CLRA.pc22.ss.qq
residqq(cm58Z, "observation") #CLRA.pc22.obs.qq


#three types of randomized quantile residuals for binomial N-mixture models
rqresiduals(cm58, type = "marginal")
rqresiduals(cm58, type="site-sum")



#ZIP 
#Qq plot of randomized quantile residuals against standard normal quantiles
residqq(cm58Z, "site-sum") #CLRA.pc21.ss.qq
residqq(cm58Z, "observation")#CLRA.pc21.ob.qq







##==============================================================================
##Plot Covs

#Take out random effects
cm58<- pcount(~date+wind+obs 
              ~EH + Watdist + Estnegdist + Updev, data=umfc, mixture= "P", engine = "TMB")

#EH
range(siteCovs(umfc)$EH) 

hplot <- data.frame(EH=seq(min(covs[,c("EH")]), max(covs[,c("EH")]), 
                                      length=262),Watdist=0, Estnegdist=0, Updev=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm58, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=EH,y=Predicted,ymin=lower,ymax=upper))

cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=EH,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=EH,y=upper), alpha=0.1)

cwwi.plot <- cwwi.plot + geom_line(aes(x=EH,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Landscape Heterogeneity")+
  scale_y_continuous(limits = c(0,70), expand = c(0,0))+
  scale_x_continuous(breaks= c(-2, -1, 0, 1), labels = c("0.0", "0.2", "0.4", "0.6"), expand=c(0,0))

eh22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(eh22.plot) 



#Watdist
range(siteCovs(umfc)$Watdist) 

hplot <- data.frame(Watdist=seq(min(covs[,c("Watdist")]), max(covs[,c("Watdist")])), length=262,
                                EH=0, Estnegdist=0, Updev=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm58, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Watdist,y=Predicted,ymin=lower,ymax=upper))

cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Watdist,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Watdist,y=upper), alpha=0.1)

cwwi.plot <- cwwi.plot + geom_line(aes(x=Watdist,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Distance to Open Water (km)")+
  scale_y_continuous(limits = c(0,70), expand = c(0,0))+
  scale_x_continuous(breaks= c(0, 2.5, 5), labels = c("0", "3", "6"),expand=c(0,0))

w22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(w22.plot) 





#Estnegdist
range(siteCovs(umfc)$Estnegdist) 

hplot <- data.frame(Estnegdist=seq(min(covs[,c("Estnegdist")]), max(covs[,c("Estnegdist")])), length=262,
                    EH=0, Watdist=0, Updev=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm58, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Estnegdist,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Estnegdist,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Estnegdist,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Estnegdist,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Distance to Estuarine Marsh (km)")+
  scale_y_continuous(limits = c(0,70), expand = c(0,0))+
  scale_x_continuous(breaks= c(-2, 0, 2), labels = c("-1", "0", "1"), expand = c(0,0))

e22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(e22.plot) 


#Updev
range(siteCovs(umfc)$Updev) 

hplot <- data.frame(Updev=seq(min(covs[,c("Updev")]), max(covs[,c("Updev")])), length=262,
                    EH=0, Watdist=0, Estnegdist=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm58, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Updev,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Updev,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Updev,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Updev,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Distance to Upland/Developed (km)")+
  scale_y_continuous(limits = c(0,70), expand = c(0,0))+
  scale_x_continuous(breaks= c(-1, 0, 1, 2, 3, 4), labels = c("0", "1", "2", "3", "4", "5"), expand = c(0,0))

u22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(u22.plot) 





