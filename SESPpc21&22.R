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

#Scale site covs
covs[,"EH"] <- scale(covs[,"EH"])
covs[,"Updev"] <- scale(covs[,"Updev"])
covs[,"Estnegdist"] <- scale(covs[,"Estnegdist"])
covs[,"Watdist"] <- scale(covs[,"Watdist"])
covs[,"Palfor"] <- scale(covs[,"Palfor"])




##===================================================================================
##Create unmarked df-----------------------------------------------------------------

umfc <- unmarkedFramePCount(y=SESPc, obsCovs = list(date=(daymat), wind=(windmat), visit=(visit), obs=(obs)), 
                            siteCovs = covs[,c("EH", "Watdist", "Estnegdist", "Updev", "Palfor", "Point_ID")])

summary(umfc)
summary(apply(SESPc,1,max,na.rm=TRUE))
#mean counts per point= 1.095
#max counts per point= 15

#Null model 
nullm <- pcount(~1~1, data = umfc, engine = "TMB")
backTransform(nullm, type="det") #0.566

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
#d12    11 1725.54  0.00 6.9e-01     0.69
#d15    11 1727.54  2.00 2.6e-01     0.95

##add abundance covs to best detection sub model 
cm1<- pcount(~date+wind+obs ~EH, data=umfc, mixture= "P", engine = "TMB")
cm2<- pcount(~date+wind+obs ~Watdist, data=umfc, mixture= "P", engine = "TMB")
cm3<- pcount(~date+wind+obs ~Estnegdist, data=umfc, mixture= "P", engine = "TMB")
cm4<- pcount(~date+wind+obs ~Updev, data=umfc, mixture= "P", engine = "TMB")
cm5<- pcount(~date+wind+obs ~Palfor, data=umfc, mixture= "P", engine = "TMB")
cm6<- pcount(~date+wind+obs ~EH + Watdist, data=umfc, mixture= "P", engine = "TMB")
cm7<- pcount(~date+wind+obs ~EH + Estnegdist, data=umfc, mixture= "P", engine = "TMB")
cm8<- pcount(~date+wind+obs ~EH + Updev, data=umfc, mixture= "P", engine = "TMB")
cm9<- pcount(~date+wind+obs ~EH + Palfor,  data=umfc, mixture= "P", engine = "TMB")
cm10<- pcount(~date+wind+obs ~Watdist + Estnegdist, data=umfc, mixture= "P", engine = "TMB")
cm11<- pcount(~date+wind+obs ~Watdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm12<- pcount(~date+wind+obs ~Watdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm13<- pcount(~date+wind+obs ~Estnegdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm14<- pcount(~date+wind+obs ~Estnegdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm15<- pcount(~date+wind+obs ~Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm16<- pcount(~date+wind+obs ~EH + Watdist + Estnegdist, data=umfc, mixture= "P", engine = "TMB")
cm17<- pcount(~date+wind+obs ~EH + Watdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm18<- pcount(~date+wind+obs ~EH + Watdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm19<- pcount(~date+wind+obs ~EH + Estnegdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm20<- pcount(~date+wind+obs ~EH + Estnegdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm21<- pcount(~date+wind+obs ~EH + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm22<- pcount(~date+wind+obs ~Watdist + Estnegdist + Updev,  data=umfc, mixture= "P", engine = "TMB")
cm23<- pcount(~date+wind+obs ~Watdist + Estnegdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm24<- pcount(~date+wind+obs ~Watdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm25<- pcount(~date+wind+obs ~Estnegdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm26<- pcount(~date+wind+obs ~EH + Watdist + Estnegdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm27<- pcount(~date+wind+obs ~EH + Watdist + Estnegdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm28<- pcount(~date+wind+obs ~EH + Watdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm29<- pcount(~date+wind+obs ~EH + Estnegdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm30<- pcount(~date+wind+obs ~Watdist + Estnegdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm31<- pcount(~date+wind+obs ~EH + Watdist + Estnegdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm32<- pcount(~date+wind+obs ~(1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm33<- pcount(~date+wind+obs ~EH + (1|Point_ID),  data=umfc, mixture= "P", engine = "TMB")
cm34<- pcount(~date+wind+obs ~Watdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm35<- pcount(~date+wind+obs ~Estnegdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm36<- pcount(~date+wind+obs ~Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm37<- pcount(~date+wind+obs ~Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm38<- pcount(~date+wind+obs ~EH + Watdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm39<- pcount(~date+wind+obs ~EH + Estnegdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm40<- pcount(~date+wind+obs ~EH + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm41<- pcount(~date+wind+obs ~EH + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm42<- pcount(~date+wind+obs ~Watdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm43<- pcount(~date+wind+obs ~Watdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm44<- pcount(~date+wind+obs ~Watdist + Estnegdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm45<- pcount(~date+wind+obs ~Estnegdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm46<- pcount(~date+wind+obs ~Estnegdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm47<- pcount(~date+wind+obs ~Palfor + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm48<- pcount(~date+wind+obs ~EH + Watdist + Estnegdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm49<- pcount(~date+wind+obs ~EH + Watdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm50<- pcount(~date+wind+obs ~EH + Watdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm51<- pcount(~date+wind+obs ~EH + Estnegdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm52<- pcount(~date+wind+obs ~EH + Estnegdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm53<- pcount(~date+wind+obs ~EH + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm54<- pcount(~date+wind+obs ~Watdist + Estnegdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm55<- pcount(~date+wind+obs ~Watdist + Estnegdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm56<- pcount(~date+wind+obs ~EH + Estnegdist + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm57<- pcount(~date+wind+obs ~Watdist + Estnegdist + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm58<- pcount(~date+wind+obs ~EH + Watdist + Estnegdist + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")

#Create fit list of all models
cmL <- modSel(fitList(cm1, cm2, cm3, cm4, cm5, cm6, cm7, cm8, cm9, cm10, 
                      cm11, cm12, cm13, cm14, cm15, cm16, cm17, cm18, cm19, cm20,
                      cm21, cm22, cm23, cm24, cm25, cm26, cm27, cm28, cm29, cm30, cm31, cm32, cm33, cm34,
                      cm35, cm36,cm37, cm38, cm39, cm40, cm41, cm42, cm43, cm45, cm46, cm47, 
                      cm48, cm49, cm50, cm51, cm52, cm53, cm54, cm55, cm56, cm57, cm58))
cmL

#cm56    15 1139.64   0.00  5.7e-01     0.57
#cm58    16 1141.61   1.97  2.1e-01     0.79

##Best Models---------------------------------------------------------------------
cm56<- pcount(~date+wind+obs 
              ~EH + Estnegdist + Updev + Palfor + (1|Point_ID), 
              data=umfc, mixture= "P", engine = "TMB")
summary(cm56)

#
#Abundance (log-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    3.296    1.815

#Fixed effects:
#  Estimate    SE     z  P(>|z|)
#(Intercept)   -0.617 0.374 -1.65 0.099187
#EH            -0.767 0.227 -3.38 0.000735
#Estnegdist    -0.693 0.187 -3.71 0.000206
#Updev          0.670 0.210  3.20 0.001396
#Palfor         0.424 0.203  2.09 0.036544

#Detection (logit-scale):
#  Estimate    SE       z  P(>|z|)
#(Intercept)          -2.03605 0.292 -6.9737 3.09e-12
#date                  0.09656 0.112  0.8609 3.89e-01
#wind1                 0.01015 0.190  0.0535 9.57e-01
#wind2                -0.00812 0.163 -0.0498 9.60e-01
#wind3                 0.22741 0.200  1.1395 2.55e-01
#wind4                -0.20636 0.223 -0.9235 3.56e-01
#wind5                -0.57396 0.346 -1.6572 9.75e-02
#obsAllister Rutledge -0.49228 0.411 -1.1977 2.31e-01
#obsBrittany Holliker  0.28073 0.324  0.8674 3.86e-01
#obsMattDunning        0.94823 0.804  1.1794 2.38e-01

#AIC: 1139.636


#ZIP  
cm56Z<- pcount(~date+wind+obs ~EH + Estnegdist + Updev + Palfor + (1|Point_ID), 
               data=umfc, mixture= "ZIP", engine = "TMB")


summary(cm56Z) 
#AIC: 1141.644 
#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    3.297    1.816

#Fixed effects:
#  Estimate    SE     z  P(>|z|)
#(Intercept)   -0.617 0.374 -1.65 0.099274
#EH            -0.767 0.227 -3.38 0.000736
#Estnegdist    -0.693 0.187 -3.71 0.000207
#Updev          0.670 0.210  3.20 0.001394
#Palfor         0.424 0.203  2.09 0.036561

#Detection (logit-scale):
#  Estimate    SE       z  P(>|z|)
#(Intercept)          -2.03622 0.292 -6.9694 3.18e-12
#date                  0.09675 0.112  0.8626 3.88e-01
#wind1                 0.00997 0.190  0.0526 9.58e-01
#wind2                -0.00841 0.163 -0.0515 9.59e-01
#wind3                 0.22707 0.200  1.1378 2.55e-01
#wind4                -0.20637 0.223 -0.9236 3.56e-01
#wind5                -0.57385 0.346 -1.6570 9.75e-02
#obsAllister Rutledge -0.49250 0.411 -1.1982 2.31e-01
#obsBrittany Holliker  0.28082 0.324  0.8676 3.86e-01
#obsMattDunning        0.94531 0.803  1.1767 2.39e-01

#Zero-inflation (logit-scale):
#  Estimate   SE     z P(>|z|)
#-9.15 16.1 -0.57   0.569

#AIC: 1141.644 





#Poisson wins! 
#But ZIP is better for GOF qq plots...

cmL <- modSel(fitList(cm56, cm56Z))
#delta AIC isnt good for ZIP though 

##=====================================================================================================
##Abundance prediction --------------------------------------------------------------------------------

#Poisson 
renb <- ranef(cm56)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#1528.184 1466.975 1590.025 



#ZIP 
renb <- ranef(cm56Z)
plot(renb, layout=c(4,3), subset=site%in%1:12, xlim=c(-1, 11), lwd=5)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#1528.949 1468.000 1598.000  




##GoF===============================================================================
##Plots-----------------------------------------------------------------------------

#Residuals
plot(cm56) #SESP.pc21.resid
plot(cm56Z)


#expected vs observed counts 
n.obs <- apply(SESPc, 1, sum, na.rm=TRUE)
n.predP <- apply(fitted(cm56), 1, sum, na.rm=TRUE )
n.predZ <- apply(fitted(cm56Z), 1, sum, na.rm=TRUE )

plot(n.obs, n.predP, frame=F, ylab="Expected Counts", xlab="Observed Counts", main="SESP pcount + random effect")
abline(0,1)
points(smooth.spline(n.predP ~ n.obs, df=4), type="l", lwd=2, col="blue")
points(smooth.spline(n.predZ ~ n.obs, df=4), type="l", lwd=2, col="green")  


#fitted vs residuals 
n.residsP <-apply(residuals(cm56),1,sum,na.rm=T)
plot(n.predP, n.residsP, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="SESP pcount+re NB")
points(smooth.spline(n.residsP ~ n.predP, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predP, df=4), type="l", lwd=2, col="black")
#SESP.pc21.fit

#fitted vs residuals 
n.residsZ <-apply(residuals(cm56Z),1,sum,na.rm=T)
plot(n.predZ, n.residsZ, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="SESP pcount+re P")
points(smooth.spline(n.residsZ ~ n.predZ, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predZ, df=4), type="l", lwd=2, col="black")
#SESPpc.p.fit

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



gof.p <- parboot(cm56, fitstats, nsim= 1000)
gof.p

#Parametric Bootstrap Statistics:
#t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
#SSE          386.9          70.53            42.93       0.0609
#Chisq        244.7          -5.82            33.26       0.5285
#freemanTukey  90.5           3.10             7.39       0.3307

#compute c-hat 
gof.p@t0[2] / mean(gof.p@t.star[,2])
#Chisq 
#0.9767711 





##Nmixgof-------------------------------------------------------------------------------

#P
#chat 
chat(cm56, type = "marginal")


#Residuals against covariates
residcov(cm56)

#Residual against fitted 
residfit(cm56, type="marginal")
residfit(cm56, type="site-sum")
residfit(cm56, type="observation")

#Qq plot of randomized quantile residuals against standard normal quantiles
residqq(cm56, "site-sum") #SESP.pc21.ss.qq
residqq(cm56, "observation") #SESP21.pc.obs.qq


#three types of randomized quantile residuals for binomial N-mixture models
rqresiduals(cm56, type = "marginal")
rqresiduals(cm56, type="site-sum")



#ZIP 
#Qq plot of randomized quantile residuals against standard normal quantiles
residqq(cm56Z, "site-sum")
residqq(cm56Z, "observation")


##=============================================================================
#Plot Covs

#take out random effects
cm56Z<- pcount(~date+wind+obs 
              ~EH + Estnegdist + Updev + Palfor, 
              data=umfc, mixture= "P", engine = "TMB")
#EH

hplot <- data.frame(EH=seq(min(covs[,c("EH")]), max(covs[,c("EH")])),length=262,
                    Updev=0, Watdist=0, Estnegdist=0, Palfor=0, 
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm56Z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=EH,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=EH,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=EH,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=EH,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Landscape Heterogeneity")+
scale_x_continuous(breaks= c(-2, -1, 0, 1), labels = c("0", "0.025", "0.050", "0.075"), expand = c(0,0))+
  scale_y_continuous(limits = c(0,10), expand = c(0,0), breaks= c(0, 2.5, 5.0, 7.5, 10),
                     labels = c("0", "3", "5", "8", "10"))


eh21.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(eh21.plot)





#Estnegdist
range(siteCovs(umfc)$Estnegdist) 

hplot <- data.frame(Estnegdist=seq(min(covs[,c("Estnegdist")]), max(covs[,c("Estnegdist")])), length=262,
                    EH=0, Watdist=0, Updev=0, Palfor=0,
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm56, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Estnegdist,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Estnegdist,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Estnegdist,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Estnegdist,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Distance to Estuarine Marsh (km)")+
  scale_x_continuous(breaks= c(-2, 0, 2), labels = c("-1", "0", "1"), expand = c(0,0))+
  scale_y_continuous(limits = c(0,10), expand = c(0,0))

e21.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(e21.plot) 


#Updev
range(siteCovs(umfc)$Updev) 

hplot <- data.frame(Updev=seq(min(covs[,c("Updev")]), max(covs[,c("Updev")])), length=262,
                    EH=0, Watdist=0, Estnegdist=0, Palfor=0,
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm56, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Updev,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Updev,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Updev,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Updev,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Distance to Upland/Developed (km)")+
  scale_x_continuous(breaks= c(-1, 0, 1, 2, 3, 4), labels = c("0", "1", "2", "3", "4", "5"), expand = c(0,0))+
  scale_y_continuous(limits = c(0,10), expand = c(0,0), breaks= c(0, 2.5, 5.0, 7.5, 10),
                     labels = c("0", "3", "5", "8", "10"))

u21.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(u21.plot) 




#Palfor
range(siteCovs(umfc)$Palfor) 

hplot <- data.frame(Palfor=seq(min(covs[,c("Palfor")]), max(covs[,c("Palfor")])), length=262,
                    EH=0, Watdist=0, Estnegdist=0, Updev=0,
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm56, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Palfor,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Palfor,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Palfor,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Palfor,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Distance to Palustrine Wetlands (km)")+
scale_x_continuous(breaks= c(-1, 0, 1, 2, 3, 4), labels = c("0", "1", "2", "3", "4","5"), expand = c(0,0))+
scale_y_continuous(limits = c(0,10), expand = c(0,0))

p21.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(p21.plot) 













#===================================================================================================
#2022 


##Read data 
#SESP counts and covariates in wide format 
#262 point count sites in total
#3 visits in year 2021 & 2022 
SESP.w <- read.table("pc.data.csv", header=T, sep=",")   
head(SESP.w)
str(SESP.w)

#select the first year 
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

#Scale site covs
covs[,"EH"] <- scale(covs[,"EH"])
covs[,"Updev"] <- scale(covs[,"Updev"])
covs[,"Estnegdist"] <- scale(covs[,"Estnegdist"])
covs[,"Watdist"] <- scale(covs[,"Watdist"])
covs[,"Palfor"] <- scale(covs[,"Palfor"])



##===================================================================================
##Create unmarked df-----------------------------------------------------------------

umfc <- unmarkedFramePCount(y=SESPc, obsCovs = list(date=(daymat), wind=(windmat), visit=(visit), obs=(obs)), 
                            siteCovs = covs[,c("EH", "Watdist", "Estnegdist", "Updev", "Palfor", "Point_ID")])

summary(umfc)
summary(apply(SESPc,1,max,na.rm=TRUE))
#mean counts per point= 0.5992 
#max counts per point= 10

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
#Point_ID (Intercept)    6.521    2.554

#Fixed effects:
#  Estimate    SE      z P(>|z|)
#(Intercept)   -1.161 0.597 -1.945 0.05172   ns
#EH            -0.476 0.354 -1.344 0.17893   ns
#Watdist        0.296 0.311  0.950 0.34187   ns
#Estnegdist    -0.458 0.337 -1.362 0.17324   ns
#Updev         -0.104 0.301 -0.347 0.72849   ns
#Palfor         0.859 0.263  3.269 0.00108

#Detection (logit-scale):
#  Random effects:
#  Groups        Name Variance Std.Dev.
#visit (Intercept)    0.008    0.087

#Fixed effects:
#  Estimate     SE      z  P(>|z|)
#(Intercept)           -5.007  0.673 -7.443 9.82e-14
#date                  -0.559  0.187 -2.997 2.73e-03
#wind1                  0.210  0.251  0.837 4.02e-01
#wind2                  0.635  0.284  2.234 2.55e-02
#wind3                  0.209  0.291  0.719 4.72e-01
#wind4                  0.827  0.405  2.040 4.14e-02
#wind5                 -7.987 96.218 -0.083 9.34e-01
#obsDanielle Canning    1.214  0.757  1.603 1.09e-01
#obsFred Hall           2.048  0.678  3.020 2.53e-03
#obsMiranda Mordue      1.933  0.658  2.940 3.29e-03
#obsSofia Campuzano    -5.477 50.737 -0.108 9.14e-01
#obsTyler Connell       0.312  1.254  0.249 8.03e-01

#AIC: 707.9624 

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
#d15    13 1081.40   0.00 6.8e-01     0.68
#d12    13 1082.90   1.49 3.2e-01     1.00

##add abundance covs to best detection sub model 
cm1<- pcount(~date+wind+obs ~EH, data=umfc, mixture= "P", engine = "TMB")
cm2<- pcount(~date+wind+obs ~Watdist, data=umfc, mixture= "P", engine = "TMB")
cm3<- pcount(~date+wind+obs ~Estnegdist, data=umfc, mixture= "P", engine = "TMB")
cm4<- pcount(~date+wind+obs ~Updev, data=umfc, mixture= "P", engine = "TMB")
cm5<- pcount(~date+wind+obs ~Palfor, data=umfc, mixture= "P", engine = "TMB")
cm6<- pcount(~date+wind+obs ~EH + Watdist, data=umfc, mixture= "P", engine = "TMB")
cm7<- pcount(~date+wind+obs ~EH + Estnegdist, data=umfc, mixture= "P", engine = "TMB")
cm8<- pcount(~date+wind+obs ~EH + Updev, data=umfc, mixture= "P", engine = "TMB")
cm9<- pcount(~date+wind+obs ~EH + Palfor,  data=umfc, mixture= "P", engine = "TMB")
cm10<- pcount(~date+wind+obs ~Watdist + Estnegdist, data=umfc, mixture= "P", engine = "TMB")
cm11<- pcount(~date+wind+obs ~Watdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm12<- pcount(~date+wind+obs ~Watdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm13<- pcount(~date+wind+obs ~Estnegdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm14<- pcount(~date+wind+obs ~Estnegdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm15<- pcount(~date+wind+obs ~Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm16<- pcount(~date+wind+obs ~EH + Watdist + Estnegdist, data=umfc, mixture= "P", engine = "TMB")
cm17<- pcount(~date+wind+obs ~EH + Watdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm18<- pcount(~date+wind+obs ~EH + Watdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm19<- pcount(~date+wind+obs ~EH + Estnegdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm20<- pcount(~date+wind+obs ~EH + Estnegdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm21<- pcount(~date+wind+obs ~EH + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm22<- pcount(~date+wind+obs ~Watdist + Estnegdist + Updev,  data=umfc, mixture= "P", engine = "TMB")
cm23<- pcount(~date+wind+obs ~Watdist + Estnegdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm24<- pcount(~date+wind+obs ~Watdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm25<- pcount(~date+wind+obs ~Estnegdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm26<- pcount(~date+wind+obs ~EH + Watdist + Estnegdist + Updev, data=umfc, mixture= "P", engine = "TMB")
cm27<- pcount(~date+wind+obs ~EH + Watdist + Estnegdist + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm28<- pcount(~date+wind+obs ~EH + Watdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm29<- pcount(~date+wind+obs ~EH + Estnegdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm30<- pcount(~date+wind+obs ~Watdist + Estnegdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm31<- pcount(~date+wind+obs ~EH + Watdist + Estnegdist + Updev + Palfor, data=umfc, mixture= "P", engine = "TMB")
cm32<- pcount(~date+wind+obs ~(1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm33<- pcount(~date+wind+obs ~EH + (1|Point_ID),  data=umfc, mixture= "P", engine = "TMB")
cm34<- pcount(~date+wind+obs ~Watdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm35<- pcount(~date+wind+obs ~Estnegdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm36<- pcount(~date+wind+obs ~Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm37<- pcount(~date+wind+obs ~Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm38<- pcount(~date+wind+obs ~EH + Watdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm39<- pcount(~date+wind+obs ~EH + Estnegdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm40<- pcount(~date+wind+obs ~EH + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm41<- pcount(~date+wind+obs ~EH + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm42<- pcount(~date+wind+obs ~Watdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm43<- pcount(~date+wind+obs ~Watdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm44<- pcount(~date+wind+obs ~Watdist + Estnegdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm45<- pcount(~date+wind+obs ~Estnegdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm46<- pcount(~date+wind+obs ~Estnegdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm47<- pcount(~date+wind+obs ~Palfor + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm48<- pcount(~date+wind+obs ~EH + Watdist + Estnegdist + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm49<- pcount(~date+wind+obs ~EH + Watdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm50<- pcount(~date+wind+obs ~EH + Watdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm51<- pcount(~date+wind+obs ~EH + Estnegdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm52<- pcount(~date+wind+obs ~EH + Estnegdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm53<- pcount(~date+wind+obs ~EH + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm54<- pcount(~date+wind+obs ~Watdist + Estnegdist + Updev + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm55<- pcount(~date+wind+obs ~Watdist + Estnegdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm56<- pcount(~date+wind+obs ~EH + Estnegdist + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm57<- pcount(~date+wind+obs ~Watdist + Estnegdist + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm58<- pcount(~date+wind+obs ~EH + Watdist + Estnegdist + Updev + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")

#Create fit list of all models
cmL <- modSel(fitList(cm1, cm2, cm3, cm4, cm5, cm6, cm7, cm8, cm9, cm10, 
                      cm11, cm12, cm13, cm14, cm15, cm16, cm17, cm18, cm19, cm20,
                      cm21, cm22, cm23, cm24, cm25, cm26, cm27, cm28, cm29, cm30, cm31, cm32, cm33, cm34,
                      cm35, cm36,cm37, cm38, cm39, cm40, cm41, cm42, cm43, cm45, cm46, cm47, 
                      cm48, cm49, cm50, cm51, cm52, cm53, cm54, cm55, cm56, cm57, cm58))
cmL

#cm52    16  703.12   0.00 2.2e-01     0.22
#cm55    16  703.98   0.85 1.4e-01     0.36
#cm45    15  704.04   0.91 1.4e-01     0.50
#cm49    16  704.18   1.05 1.3e-01     0.63
#cm43    15  704.91   1.79 9.0e-02     0.72
#cm56    17  704.97   1.85 8.7e-02     0.81

##Best Models---------------------------------------------------------------------
cm52<- pcount(~date+wind+obs ~EH + Estnegdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm55<- pcount(~date+wind+obs ~Watdist + Estnegdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")
cm45<- pcount(~date+wind+obs ~Estnegdist + Palfor + (1|Point_ID), data=umfc, mixture= "P", engine = "TMB")

summary(cm52)
summary(cm55)
summary(cm45)

#
#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    7.717    2.778

#Fixed effects:
#  Estimate    SE     z P(>|z|)
#(Intercept)   -1.337 0.634 -2.11 0.03505
#Estnegdist    -0.909 0.258 -3.53 0.00042
#Palfor         0.764 0.239  3.20 0.00137

#Detection (logit-scale):
#  Estimate      SE       z  P(>|z|)
#(Intercept)           -4.963   0.642 -7.7260 1.11e-14
#date                  -0.571   0.182 -3.1281 1.76e-03
#wind1                  0.229   0.246  0.9312 3.52e-01
#wind2                  0.668   0.273  2.4473 1.44e-02
#wind3                  0.213   0.288  0.7377 4.61e-01
#wind4                  0.855   0.395  2.1640 3.05e-02
#wind5                 -8.127 106.927 -0.0760 9.39e-01
#obsDanielle Canning    1.083   0.694  1.5595 1.19e-01
#obsFred Hall           2.032   0.658  3.0891 2.01e-03
#obsMiranda Mordue      1.865   0.630  2.9607 3.07e-03
#obsSofia Campuzano    -5.932  67.874 -0.0874 9.30e-01
#obsTyler Connell       0.191   1.225  0.1558 8.76e-01

#AIC: 704.0387 


#ZIP  
cm45Z<- pcount(~date+wind+obs ~Estnegdist + Palfor + (1|Point_ID), 
               data=umfc, mixture= "ZIP", engine = "TMB")



summary(cm45Z)  
#AIC: 704.7686 
#Abundance (log-scale):
#Random effects:
#  Groups        Name Variance Std.Dev.
#Point_ID (Intercept)    0.498    0.705

#Fixed effects:
#  Estimate    SE     z  P(>|z|)
#(Intercept)    1.641 0.345  4.75 1.99e-06
#Estnegdist    -0.788 0.142 -5.56 2.64e-08
#Palfor         0.573 0.127  4.50 6.73e-06

#Detection (logit-scale):
#  Estimate      SE       z  P(>|z|)
#(Intercept)          -4.2232   0.482 -8.7582 1.98e-18
#date                 -0.6435   0.138 -4.6636 3.11e-06
#wind1                 0.1991   0.231  0.8617 3.89e-01
#wind2                 0.7946   0.265  2.9983 2.72e-03
#wind3                 0.3318   0.287  1.1551 2.48e-01
#wind4                 0.8581   0.373  2.3019 2.13e-02
#wind5                -8.7175 197.051 -0.0442 9.65e-01
#obsDanielle Canning   0.0805   0.500  0.1610 8.72e-01
#obsFred Hall          1.1390   0.439  2.5916 9.55e-03
#obsMiranda Mordue     1.2449   0.450  2.7647 5.70e-03
#obsSofia Campuzano  -11.0796 479.280 -0.0231 9.82e-01
#obsTyler Connell     -0.4588   1.116 -0.4113 6.81e-01

#Zero-inflation (logit-scale):
#  Estimate   SE    z P(>|z|)
#0.518 0.29 1.79  0.0739

#AIC: 704.7686 


#very close, so compare below  


##=====================================================================================================
##Abundance prediction --------------------------------------------------------------------------------

#Poisson 
renb <- ranef(cm45)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
#1314.453 1253.975 1377.000 




#ZIP 
renb <- ranef(cm45Z)
plot(renb, layout=c(4,3), subset=site%in%1:12, xlim=c(-1, 11), lwd=5)

N.total.post <- predict(renb, func=sum, nsim=1000)

hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")

c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025), quantile(N.total.post, prob=0.975))
#Estimate     2.5%    97.5% 
# 615.338  563.000  673.025  




##GoF===============================================================================
##Plots-----------------------------------------------------------------------------

#Residuals
plot(cm45) 
plot(cm45Z)


#expected vs observed counts 
n.obs <- apply(SESPc, 1, sum, na.rm=TRUE)
n.predP <- apply(fitted(cm45), 1, sum, na.rm=TRUE )
n.predZ <- apply(fitted(cm45Z), 1, sum, na.rm=TRUE )

plot(n.obs, n.predP, frame=F, ylab="Expected Counts", xlab="Observed Counts", main="SESP pcount + random effect")
abline(0,1)
points(smooth.spline(n.predP ~ n.obs, df=4), type="l", lwd=2, col="blue")
points(smooth.spline(n.predZ ~ n.obs, df=4), type="l", lwd=2, col="green")  
#Error in smooth.spline(n.predP ~ n.obs, df = 4) : 
#'tol' must be strictly positive and finite

#fitted vs residuals 
n.residsP <-apply(residuals(cm45),1,sum,na.rm=T)
plot(n.predP, n.residsP, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="SESP pcount+re NB")
points(smooth.spline(n.residsP ~ n.predP, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predP, df=4), type="l", lwd=2, col="black")
#SESP.pc21.fit

#fitted vs residuals 
n.residsZ <-apply(residuals(cm45Z),1,sum,na.rm=T)
plot(n.predZ, n.residsZ, frame=F, ylab="Residuals", xlab="Fitted Values (P)", main="SESP pcount+re P")
points(smooth.spline(n.residsZ ~ n.predZ, df=4), type="l", lwd=2, col="red")
points(smooth.spline(rep(0,262) ~ n.predZ, df=4), type="l", lwd=2, col="black")
##SESP.pc22.zfit

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



gof.p <- parboot(cm45, fitstats, nsim= 1000)
gof.p



gofZ <- parboot(cm45Z, fitstats, nsim= 1000)
gofZ

#Parametric Bootstrap Statistics:
#t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
#SSE           823          307.1            435.2       0.0154
#Chisq        1019           90.1            253.7       0.3564
#freemanTukey  183           18.5             41.5       0.3003

#t_B quantiles:
#  0% 2.5% 25% 50%  75% 97.5%  100%
#SSE          39   94 403 509  609   787 11717
#Chisq        72  149 870 965 1057  1281  2041
#freemanTukey 22   40 156 171  187   218   506

#compute c-hat 
gofZ@t0[2] / mean(gofZ@t.star[,2])
#Chisq 
#1.097001  

##Nmixgof-------------------------------------------------------------------------------

#P
#chat 
chat(cm45, type = "marginal")


#Residuals against covariates
residcov(cm45)

#Residual against fitted 
residfit(cm45, type="marginal")
residfit(cm45, type="site-sum")
residfit(cm45, type="observation")

#Qq plot of randomized quantile residuals against standard normal quantiles
residqq(cm45, "site-sum") #SESP.pc21.ss.qq
residqq(cm45, "observation") #SESP21.pc.obs.qq


#three types of randomized quantile residuals for binomial N-mixture models
rqresiduals(fmc2, type = "marginal")
rqresiduals(fmc2, type="site-sum")



#ZIP 
#Qq plot of randomized quantile residuals against standard normal quantiles
residqq(cm45Z, "site-sum")
residqq(cm45Z, "observation") #SESP.pc22.z.obs.qq



##=============================================================================
#Plot Covs

#take out random effects
cm45Z<- pcount(~date+wind+obs ~Estnegdist + Palfor, 
               data=umfc, mixture= "ZIP", engine = "TMB")





#Estnegdist
range(siteCovs(umfc)$Estnegdist) 

hplot <- data.frame(Estnegdist=seq(min(covs[,c("Estnegdist")]), max(covs[,c("Estnegdist")])), length=262,
                    EH=0, Watdist=0, Updev=0, Palfor=0,
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm45Z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Estnegdist,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Estnegdist,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Estnegdist,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Estnegdist,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Distance to Estuarine Marsh (km)")+
  scale_x_continuous(breaks= c(-2, 0, 2), labels = c("-1", "0", "1"),expand = c(0,0))+
  scale_y_continuous(limits = c(0,40), expand = c(0,0))

e22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(e22.plot) 



#Palfor
range(siteCovs(umfc)$Palfor) 

hplot <- data.frame(Palfor=seq(min(covs[,c("Palfor")]), max(covs[,c("Palfor")])), length=262,
                    EH=0, Watdist=0, Estnegdist=0, Updev=0,
                    Point_ID=0, visit=0, date=0, wind=0, obs=0)

h.pred <- predict(cm45Z, type="state", newdata=hplot, appendData=TRUE)

h.plot <- data.frame(h.pred,hplot)

cwwi.plot <- ggplot(data=h.plot,aes(x=Palfor,y=Predicted,ymin=lower,ymax=upper))
cwwi.plot <- cwwi.plot + geom_ribbon(aes(x=Palfor,y=lower), alpha=0.1) +
  geom_ribbon(aes(x=Palfor,y=upper), alpha=0.1)
cwwi.plot <- cwwi.plot + geom_line(aes(x=Palfor,y=Predicted))

cwwi.plot <- cwwi.plot + theme_bw() + ylab("") +xlab("Distance to Palustrine Wetlands (km)")+
  scale_x_continuous(breaks= c(-1, 0, 1, 2, 3, 4), labels = c("0", "1", "2", "3", "4","5"), expand = c(0,0))+
  scale_y_continuous(limits = c(0,40), expand = c(0,0))

p22.plot <- cwwi.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(p22.plot) 










