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
#LEBI counts and covariates in wide format 
#262 point count sites in total
#3 visits in year 2021 & 2022 
LEBI.w <- read.table("pc.data.csv", header=T, sep=",")   
head(LEBI.w)
str(LEBI.w)

#select the first year 
covs <- LEBI.w[1:262,]

#Counts with truncation and passive listening only 
LEBIc <- as.matrix(covs[,c("LEBI1", "LEBI2", "LEBI3")])

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

umfc <- unmarkedFramePCount(y=LEBIc, obsCovs = list(date=(daymat), wind=(windmat), visit=(visit), obs=(obs)), 
                            siteCovs = covs[,c("EH", "Watdist", "Estnegdist", "Updev", "Palfor", "Point_ID")])

summary(umfc)
summary(apply(LEBIc,1,max,na.rm=TRUE))
#mean counts per point= 0.3511
#max counts per point= 3.0

#Null model 
nullm <- pcount(~1~1, data = umfc, engine = "TMB")
backTransform(nullm, type="det") #0.083

#fit global model 
fmc<- pcount(~date+wind+obs+(1|visit)
             ~EH+Watdist+Estnegdist+Updev+Palfor+(1|Point_ID), umfc, engine = "TMB")

summary(fmc)


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
#d15    11 654.14  0.00 5.9e-01     0.59
#d13     6 656.68  2.54 1.6e-01     0.75

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

#convergence issues 

##Best Models---------------------------------------------------------------------
cm58<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Estnegdist + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")

summary(cm58)

#Abundance (log-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.728    0.853

#Fixed effects:
#  Estimate    SE      z P(>|z|)
#(Intercept)    2.065 0.729  2.831 0.00463
#EH            -0.560 0.175 -3.199 0.00138
#Watdist       -0.429 0.213 -2.010 0.04443
#Estnegdist    -0.134 0.179 -0.751 0.45291
#Updev         -0.123 0.206 -0.599 0.54893
#Palfor        -0.128 0.209 -0.613 0.53994

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.062     0.25

#Fixed effects:
#  Estimate      SE       z P(>|z|)
#(Intercept)          -9.925  62.119 -0.1598  0.8731
#date                  0.214   0.130  1.6521  0.0985
#wind1                -0.435   0.282 -1.5407  0.1234
#wind2                -0.582   0.315 -1.8472  0.0647
#wind3                -0.869   0.513 -1.6944  0.0902
#wind4                -1.397   0.762 -1.8332  0.0668
#wind5               -10.052 150.371 -0.0668  0.9467
#obsMattDunning        5.697  62.114  0.0917  0.9269
#obsRachel Anderson    5.964  62.115  0.0960  0.9235
#obsSofia Campuzano    5.316  62.115  0.0856  0.9318

#AIC: 639.6473

#P without nonsig 
cm58<- pcount(~date+wind+(1|visit)+obs ~EH+Watdist+(1|Point_ID), 
              umfc, engine = "TMB", mixture="P")
summary(cm58) #AIC: 2819.957 

#Abundance (log-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.747    0.864

#Fixed effects:
#  Estimate    SE     z  P(>|z|)
#(Intercept)    2.025 0.693  2.92 0.003493
#EH            -0.515 0.154 -3.34 0.000851
#Watdist       -0.313 0.168 -1.86 0.062440    not sig 

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.059    0.242

#Fixed effects:
#  Estimate      SE       z P(>|z|)
#(Intercept)          -7.467  20.020 -0.3730  0.7092
#date                  0.221   0.127  1.7366  0.0825
#wind1                -0.456   0.282 -1.6196  0.1053
#wind2                -0.642   0.310 -2.0685  0.0386
#wind3                -0.962   0.506 -1.9012  0.0573
#wind4                -1.549   0.751 -2.0632  0.0391
#wind5               -10.691 174.587 -0.0612  0.9512
#obsMattDunning        3.376  20.009  0.1687  0.8660
#obsRachel Anderson    3.547  20.010  0.1773  0.8593
#obsSofia Campuzano    2.939  20.010  0.1469  0.8832

#AIC: 636.1343

#Need to take out watdist 
cm58<- pcount(~date+wind+obs+(1|visit) ~EH+(1|Point_ID), 
              umfc, engine = "TMB", mixture="P")
summary(cm58) #AIC: 637.9593 
#Abundance (log-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.807    0.898

#Fixed effects:
#  Estimate    SE     z P(>|z|)
#(Intercept)     2.04 0.699  2.92 0.00345
#EH             -0.36 0.130 -2.78 0.00544

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.051    0.225

#Fixed effects:
#  Estimate       SE        z P(>|z|)
#(Intercept)         -13.328  342.403 -0.03893  0.9689
#date                  0.229    0.127  1.80603  0.0709
#wind1                -0.463    0.284 -1.63325  0.1024
#wind2                -0.624    0.312 -2.00161  0.0453
#wind3                -0.974    0.508 -1.91504  0.0555
#wind4                -1.557    0.756 -2.06002  0.0394
#wind5               -18.281 7716.318 -0.00237  0.9981
#obsMattDunning        9.161  342.401  0.02675  0.9787
#obsRachel Anderson    9.492  342.401  0.02772  0.9779
#obsSofia Campuzano    8.770  342.401  0.02561  0.9796

#AIC: 637.9593 


#ZIP without nonsig 
cm58Z<- pcount(~date+wind+(1|visit)+obs ~EH+(1|Point_ID),
               umfc, engine = "TMB", mixture="ZIP")

summary(cm58Z) #AIC: 639.7237 

#Poisson wins 



##=====================================================================================================
##Abundance prediction --------------------------------------------------------------------------------

#Poisson 
renb <- ranef(cm58)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
# 2487.256 2393.975 2578.025 



#ZIP 
renb <- ranef(cm58Z)
plot(renb, layout=c(4,3), subset=site%in%1:12, xlim=c(-1, 11), lwd=5)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#4001.082 3744.950 4267.125 
 




##GoF===============================================================================
##Plots-----------------------------------------------------------------------------

#Residuals
plot(cm58) #saved as LEBI.pc21.rand.resid
plot(cm58Z)


#expected vs observed counts 
n.obs <- apply(LEBIc, 1, sum, na.rm=TRUE)
n.predP <- apply(fitted(cm58), 1, sum, na.rm=TRUE )
n.predZ <- apply(fitted(cm58Z), 1, sum, na.rm=TRUE )

plot(n.obs, n.predP, frame=F, ylab="Expected Counts", xlab="Observed Counts", main="LEBI pc21")
abline(0,1)
points(smooth.spline(n.predP ~ n.obs, df=4), type="l", lwd=2, col="blue")
points(smooth.spline(n.predZ ~ n.obs, df=4), type="l", lwd=2, col="green")  
#LEBI.pc21.ove

#P is better than Z 

#fitted vs residuals 
n.residsP <-apply(residuals(cm58),1,sum,na.rm=T)
plot(n.predP, n.residsP, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="LEBI pc + random effects")
points(smooth.spline(n.residsP ~ n.predP, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predP, df=4), type="l", lwd=2, col="black")
#LEBI.pc21.rand.fit

#fitted vs residuals 
n.resids <-apply(residuals(cm58Z),1,sum,na.rm=T)
plot(n.predZ, n.resids, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="LEBI pcount+re P")
points(smooth.spline(n.resids ~ n.predZ, df=4), type="l", lwd=2, col="red")
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
#SSE          104.8           9.89             11.1        0.191
#Chisq        439.6        -199.36            116.2        0.961
#freemanTukey  99.4          -5.23             12.8        0.686

#t_B quantiles:
#  0% 2.5% 25% 50% 75% 97.5% 100%
#SSE           60   73  87  94 102   117  134
#Chisq        290  431 551 636 738   840  914
#freemanTukey  62   79  97 105 113   128  142

#compute c-hat 
gof.p@t0[2] / mean(gof.p@t.star[,2])
#Chisq 
#0.6880042  
 


##Nmixgof-------------------------------------------------------------------------------

#P
#chat 
chat(cm58, type = "marginal")
#[1] 0.5898854


#Residuals against covariates
residcov(fmc2)

#Residual against fitted 
residfit(fmc2, type="marginal")
residfit(fmc2, type="site-sum")
residfit(fmc2, type="observation")

#Qq plot of randomized quantile residuals against standard normal quantiles
residqq(cm58, "site-sum") #LEBI.pc21.ss.qq
residqq(cm58, "observation") 


#three types of randomized quantile residuals for binomial N-mixture models
rqresiduals(fmc2, type = "marginal")
rqresiduals(fmc2, type="site-sum")


##=============================================================================
#Plot Covs

#take out random effects
cm58<- pcount(~date+wind+obs ~EH, 
              umfc, engine = "TMB", mixture="P")
#EH

hplot <- data.frame(EH=seq(min(covs[,c("EH")]), max(covs[,c("EH")])),length=262,
                    Updev=0, Watdist=0, Estnegdist=0, Palfor=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm58, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=EH,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=EH,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=EH,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=EH,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Landscape Heterogeneity")+
scale_x_continuous(breaks= c(-2, -1, 0, 1), labels = c("0", "0.025", "0.050", "0.075"), expand = c(0,0))+
scale_y_continuous(limits = c(0,25), expand = c(0,0))

eh21.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


plot(eh21.plot)







##===============================================================================
#2022 data 

##Read data 
#LEBI counts and covariates in wide format 
#262 point count sites in total
#3 visits in year 2021 & 2022 
LEBI.w <- read.table("pc.data.csv", header=T, sep=",")   
head(LEBI.w)
str(LEBI.w)

#select the second year 
covs <- LEBI.w[263:524,]

#Counts with truncation and passive listening only 
LEBIc <- as.matrix(covs[,c("LEBI1", "LEBI2", "LEBI3")])

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

umfc <- unmarkedFramePCount(y=LEBIc, obsCovs = list(date=(daymat), wind=(windmat), visit=(visit), obs=(obs)), 
                            siteCovs = covs[,c("EH", "Watdist", "Estnegdist", "Updev", "Palfor", "Point_ID")])

summary(umfc)
summary(apply(LEBIc,1,max,na.rm=TRUE))
#mean counts per point= 0.4885
#max counts per point= 4

#Null model 
nullm <- pcount(~1~1, data = umfc, engine = "TMB")
backTransform(nullm, type="det") #0.275

#fit global model 
fmc<- pcount(~date+wind+obs+(1|visit)
             ~EH+Watdist+Estnegdist+Updev+Palfor+(1|Point_ID), umfc, engine = "TMB")

summary(fmc)

#AAbundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.337    0.581

#Fixed effects:
#  Estimate    SE      z P(>|z|)
#(Intercept)   2.8084 0.881  3.187 0.00144
#EH           -0.3017 0.130 -2.316 0.02057
#Watdist      -0.5229 0.175 -2.990 0.00279
#Estnegdist   -0.3159 0.134 -2.350 0.01877
#Updev         0.0952 0.146  0.652 0.51423    ns
#Palfor       -0.1550 0.139 -1.118 0.26360    ns 

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.262    0.512

#Fixed effects:
#  Estimate       SE       z  P(>|z|)
#(Intercept)         -17.560 1047.248 -0.0168 0.986622
#date                 -0.372    0.112 -3.3116 0.000927
#wind1                 0.102    0.208  0.4903 0.623954
#wind2                -0.346    0.280 -1.2384 0.215579
#wind3                -0.336    0.357 -0.9411 0.346648
#wind4                -1.557    0.754 -2.0641 0.039008
#wind5               -11.042  215.501 -0.0512 0.959133
#obsJamie Russell     12.350 1047.246  0.0118 0.990591
#obsRachel Anderson   13.266 1047.246  0.0127 0.989893
#obsTyler Connell     12.552 1047.246  0.0120 0.990437

#AIC: 795.8747   

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
#d15    11 806.36  0.00 8.9e-01     0.89
#d13     6 810.69  4.33 1.0e-01     1.00

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

# nPars    AIC delta   AICwt cumltvWt
#cm48    14 793.21  0.00 0.37596     0.38
#cm58    16 795.87  2.66 0.09922     0.48 

##Best Models---------------------------------------------------------------------
#changed cm48 to fmc 
fmc<- pcount(~date+wind+(1|visit)+obs ~EH + Watdist + Estnegdist + (1|Point_ID), 
                    data=umfc, mixture= "P", engine = "TMB")
summary(fmc)
#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.331    0.575

#Fixed effects:
#  Estimate    SE     z P(>|z|)
#(Intercept)    2.751 1.853  1.48 0.13759
#EH            -0.320 0.119 -2.69 0.00718
#Watdist       -0.483 0.168 -2.88 0.00397
#Estnegdist    -0.282 0.130 -2.17 0.02983     

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.265    0.515

#Fixed effects:
#  Estimate     SE      z  P(>|z|)
#(Intercept)         -9.2933 17.867 -0.520 0.602973
#date                -0.3712  0.106 -3.491 0.000482
#wind1                0.0722  0.207  0.349 0.726755
#wind2               -0.3826  0.278 -1.377 0.168570
#wind3               -0.3619  0.353 -1.024 0.305656
#wind4               -1.5937  0.751 -2.121 0.033943
#wind5               -8.5376 59.524 -0.143 0.885949
#obsJamie Russell     4.2243 17.765  0.238 0.812049
#obsRachel Anderson   5.0557 17.764  0.285 0.775955
#obsTyler Connell     4.3638 17.765  0.246 0.805963

#AIC: 793.2104 



#ZIP without nonsig 
fmcZ<- pcount(~date+wind+(1|visit)+obs ~EH+Watdist+Estnegdist+(1|Point_ID),
               umfc, engine = "TMB", mixture="ZIP")

summary(fmcZ) #AIC: 791.7366 

#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.005    0.073

#Fixed effects:
#  Estimate    SE     z  P(>|z|)
#(Intercept)    3.554 0.791  4.50 6.95e-06
#EH            -0.300 0.114 -2.65 8.15e-03
#Watdist       -0.490 0.161 -3.04 2.38e-03
#Estnegdist    -0.307 0.124 -2.47 1.34e-02

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.286    0.535

#Fixed effects:
#  Estimate      SE       z  P(>|z|)
#(Intercept)        -15.3624 316.649 -0.0485 0.961305
#date                -0.3977   0.109 -3.6538 0.000258
#wind1                0.0367   0.201  0.1827 0.855003
#wind2               -0.4627   0.277 -1.6725 0.094428
#wind3               -0.4877   0.350 -1.3926 0.163731
#wind4               -1.6981   0.748 -2.2689 0.023277
#wind5              -14.0565 898.923 -0.0156 0.987524
#obsJamie Russell     9.9685 316.648  0.0315 0.974886
#obsRachel Anderson  10.7733 316.648  0.0340 0.972859
#obsTyler Connell    10.1605 316.648  0.0321 0.974402

#Zero-inflation (logit-scale):
#  Estimate    SE     z P(>|z|)
#-1.09 0.459 -2.37  0.0176

#AIC: 791.7366

#ZIP wins 


##=====================================================================================================
##Abundance prediction --------------------------------------------------------------------------------

#Poisson 
renb <- ranef(fmc)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
# 4598.594 4472.000 4732.000 




#ZIP 
renb <- ranef(fmcZ)
plot(renb, layout=c(4,3), subset=site%in%1:12, xlim=c(-1, 11), lwd=5)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#5775.443 5426.000 6108.000 
 





##GoF===============================================================================
##Plots-----------------------------------------------------------------------------

#Residuals
plot(fmc) #saved as LEBI.pc22.rand.resid
plot(fmcZ)


#expected vs observed counts 
n.obs <- apply(LEBIc, 1, sum, na.rm=TRUE)
n.predP <- apply(fitted(fmc), 1, sum, na.rm=TRUE )
n.predZ <- apply(fitted(fmcZ), 1, sum, na.rm=TRUE )

plot(n.obs, n.predP, frame=F, ylab="Expected Counts", xlab="Observed Counts", main="LEBI pc 22")
abline(0,1)
points(smooth.spline(n.predP ~ n.obs, df=4), type="l", lwd=2, col="blue")
points(smooth.spline(n.predZ ~ n.obs, df=4), type="l", lwd=2, col="green")  
#LEBI.pc22.ove

#P is better than Z 

#fitted vs residuals 
n.residsP <-apply(residuals(fmc),1,sum,na.rm=T)
plot(n.predP, n.residsP, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="LEBI pc + random effects")
points(smooth.spline(n.residsP ~ n.predP, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predP, df=4), type="l", lwd=2, col="black")
#LEBI.pc22.rand.fit

#fitted vs residuals 
n.resids <-apply(residuals(fmcZ),1,sum,na.rm=T)
plot(n.predZ, n.resids, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="LEBI pc 22")
points(smooth.spline(n.resids ~ n.predZ, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predZ, df=4), type="l", lwd=2, col="black")
#LEBI.pc22.fit

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



gof.z <- parboot(fmcZ, fitstats, nsim= 1000)
gof.z

#Parametric Bootstrap Statistics:
#t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
#SSE          188          16.22             21.1        0.235
#Chisq        832          65.81             75.8        0.159
#freemanTukey 161           6.76             11.5        0.282

#t_B quantiles:
#  0% 2.5% 25% 50% 75% 97.5% 100%
#SSE          113  130 158 171 187   214  250
#Chisq        463  583 727 769 809   907 1044
#freemanTukey 111  132 147 155 163   176  196


#compute c-hat 
gof.z@t0[2] / mean(gof.z@t.star[,2])
#Chisq 
#0.7191204 



##Nmixgof-------------------------------------------------------------------------------

#P
#chat 
chat(cm58, type = "marginal")
#[1] 0.5898854


#Residuals against covariates
residcov(fmc2)

#Residual against fitted 
residfit(fmc2, type="marginal")
residfit(fmc2, type="site-sum")
residfit(fmc2, type="observation")

#Qq plot of randomized quantile residuals against standard normal quantiles
residqq(fmcZ, "site-sum") #LEBI.pc22.ss.qq
residqq(fmcZ, "observation") #LEBI.pc22.obs.qq


#three types of randomized quantile residuals for binomial N-mixture models
rqresiduals(fmc2, type = "marginal")
rqresiduals(fmc2, type="site-sum")


##=============================================================================
#Plot Covs

#take out random effects
fmcZ<- pcount(~date+wind+obs ~EH+Watdist+Estnegdist,
              umfc, engine = "TMB", mixture="ZIP")
#EH

hplot <- data.frame(EH=seq(min(covs[,c("EH")]), max(covs[,c("EH")])),length=262,
                    Updev=0, Watdist=0, Estnegdist=0, Palfor=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(fmcZ, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=EH,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=EH,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=EH,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=EH,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Landscape Heterogeneity")+
 scale_x_continuous(breaks= c(-2, -1, 0, 1), labels = c("0", "0.025", "0.050", "0.075"),expand = c(0,0))+
  scale_y_continuous(limits = c(0,200), expand = c(0,0))

eh22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(eh22.plot)

#Watdist
range(siteCovs(umfc)$Watdist) 

hplot <- data.frame(Watdist=seq(min(covs[,c("Watdist")]), max(covs[,c("Watdist")])), length=262,
                    EH=0, Estnegdist=0, Updev=0,Palfor=0,
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(fmcZ, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Watdist,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Watdist,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Watdist,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Watdist,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Distance to Open Water (km)")+
scale_x_continuous(breaks= c(0, 2.5, 5), labels = c("0", "3", "6"),expand = c(0,0))+
scale_y_continuous(limits = c(0,200), expand = c(0,0))

w22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(w22.plot) 


#Estnegdist
range(siteCovs(umfc)$Estnegdist) 

hplot <- data.frame(Estnegdist=seq(min(covs[,c("Estnegdist")]), max(covs[,c("Estnegdist")])), length=262,
                    EH=0, Watdist=0, Updev=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(fmcZ, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Estnegdist,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Estnegdist,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Estnegdist,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Estnegdist,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Distance to Estuarine Marsh (km)")+
  scale_x_continuous(breaks= c(-2, 0, 2), labels = c("-1", "0", "1"),expand = c(0,0))+
  scale_y_continuous(limits = c(0,200), expand = c(0,0))

e22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(e22.plot) 

