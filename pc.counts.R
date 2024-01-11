##Formatting counts for pcount analysis 

setwd("C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/PCOuntAnalysis/Data")

library(tidyverse)
library(reshape2)

##Bring in your total dataset and take a look 
pointc <- read.csv("C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/Data/allQuery_CLEAN.csv")
str(pointc) 
names(pointc)

#create column of 1s for detections 
pointc$Count <- rep(c(1))


##CLRA--------------------------------------------------------------------------
#Need to set count to 1 for species of interest and all the rest to 0
#This keeps all of the sites instead of subsetting out 
pointc$Count[which(pointc$Species != "CLRA")] <- 0
pointc$Count[is.na(pointc$Species)] <- 0
pointc$Count <- as.numeric(pointc$Count)

# Truncate distances based on histogram of distance data 
pointc$Count[which(pointc$Distance >= 230 )] <- 0

#aggregate counts by species per site, visit, and year 
#First subset by year 
#2021 
c21 <- subset(pointc, pointc$Survey_Year==2021)
#2022
c22 <- subset(pointc, pointc$Survey_Year==2022)

#Now select for point_id, count, visit, and year 
clra21 <- c21[,c("Point_ID", "Count", "Visit_Number", "Survey_Year")]
clra22 <- c22[,c("Point_ID", "Count", "Visit_Number", "Survey_Year")]

#aggregate counts by point_ID
#Aggregate sum of % by point ID 
sumc21<- aggregate(Count ~ Point_ID+Visit_Number, clra21, FUN=sum)
sumc22<- aggregate(Count ~ Point_ID+Visit_Number, clra22, FUN=sum)

#dcast is going to "cast" the count data in a wide format by visit number  
c21.w <- dcast(sumc21, Point_ID~Visit_Number)
c22.w <- dcast(sumc22, Point_ID~Visit_Number)

#This is essentialy the format you want for every obervational covarite as well (Each site is a row and the three columns are each  visit populated with the variable of interest, here it was just the variable of counts)


#bind 
clracounts <- rbind(c21.w, c22.w)
write.csv(clracounts,"C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/PCOuntAnalysis/Data\\CLRAc.csv", row.names = TRUE)
#with truncation 
write.csv(clracounts,"C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/PCOuntAnalysis/Data\\CLRActrunc.csv", row.names = TRUE)


##SESP--------------------------------------------------------------------------
#Need to set count to 1 for species of interest and all the rest to 0
#This keeps all of the sites instead of subsetting out 
pointc$Count[which(pointc$Species != "SESP")] <- 0
pointc$Count[is.na(pointc$Species)] <- 0
pointc$Count <- as.numeric(pointc$Count)

#aggregate counts by species per site, visit, and year 
#First subset by year 
#2021 
s21 <- subset(pointc, pointc$Survey_Year==2021)
#2022
s22 <- subset(pointc, pointc$Survey_Year==2022)

#Now select for point_id, count, visit, and year 
sesp21 <- s21[,c("Point_ID", "Count", "Visit_Number", "Survey_Year")]
sesp22 <- s22[,c("Point_ID", "Count", "Visit_Number", "Survey_Year")]

#aggregate counts by point_ID
#Aggregate sum of % by point ID 
sums21<- aggregate(Count ~ Point_ID+Visit_Number, sesp21, FUN=sum)
sums22<- aggregate(Count ~ Point_ID+Visit_Number, sesp22, FUN=sum)

#dcast is going to "cast" the count data in a wide format by visit number  
s21.w <- dcast(sums21, Point_ID~Visit_Number)
s22.w <- dcast(sums22, Point_ID~Visit_Number)

#This is essentialy the format you want for every obervational covarite as well (Each site is a row and the three columns are each  visit populated with the variable of interest, here it was just the variable of counts)


#bind 
SESPcounts <- rbind(s21.w, s22.w)
write.csv(SESPcounts,"C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/PCOuntAnalysis/Data\\SESPc.csv", row.names = TRUE)

##LEBI--------------------------------------------------------------------------
#Need to set count to 1 for species of interest and all the rest to 0
#This keeps all of the sites instead of subsetting out 
pointc$Count[which(pointc$Species != "LEBI")] <- 0
pointc$Count[is.na(pointc$Species)] <- 0
pointc$Count <- as.numeric(pointc$Count)

#aggregate counts by species per site, visit, and year 
#First subset by year 
#2021 
l21 <- subset(pointc, pointc$Survey_Year==2021)
#2022
l22 <- subset(pointc, pointc$Survey_Year==2022)

#Now select for point_id, count, visit, and year 
lebi21 <- l21[,c("Point_ID", "Count", "Visit_Number", "Survey_Year")]
lebi22 <- l22[,c("Point_ID", "Count", "Visit_Number", "Survey_Year")]

#aggregate counts by point_ID
#Aggregate sum of % by point ID 
suml21<- aggregate(Count ~ Point_ID+Visit_Number, lebi21, FUN=sum)
suml22<- aggregate(Count ~ Point_ID+Visit_Number, lebi22, FUN=sum)

#dcast is going to "cast" the count data in a wide format by visit number  
l21.w <- dcast(suml21, Point_ID~Visit_Number)
l22.w <- dcast(suml22, Point_ID~Visit_Number)
#This is essentialy the format you want for every obervational covarite as well (Each site is a row and the three columns are each  visit populated with the variable of interest, here it was just the variable of counts)


#bind 
LEBIcounts <- rbind(l21.w, l22.w)

#Export
write.csv(LEBIcounts,"C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/PCOuntAnalysis/Data\\LEBIpc.csv", row.names = TRUE)







#------------------------------------------------------------------------------
#Passive listening period only 

#Subset for passive listening period 
pass <- pointc[,c("Point_ID","Visit_Number", "Species", "Count","Pass1", "Pass2", "Pass3", "Pass4", "Pass5", "Survey_Year", "Distance")]

#Need to correct for blank cells by creating column that identifies no detections in passive period 
pass$Detect <- ifelse(pass$Pass1 == "" & pass$Pass2 == "" & pass$Pass3 == "" & pass$Pass4 == "" & pass$Pass5 == "" , "No", "Yes")


##CLRA--------------------------------------------------------------------------
#Need to set count to 1 for species of interest and all the rest to 0
#This keeps all of the sites instead of subsetting out 
pass$Count[which(pass$Species != "CLRA")] <- 0

#Need Count to be 0 for species = na & when detect = no 
pass$Count[is.na(pass$Species)] <- 0
pass$Count[which(pass$Detect == "No")] <- 0

#Make sure count is treated as numeric 
pass$Count <- as.numeric(pass$Count)

# Truncate distances based on histogram of distance data 
pass$Count[which(pass$Distance >= 230 )] <- 0

#aggregate counts by species per site, visit, and year 
#First subset by year 
#2021 
c21 <- subset(pass, pass$Survey_Year==2021)
#2022
c22 <- subset(pass, pass$Survey_Year==2022)

#Now select for point_id, count, visit, and year 
clra21 <- c21[,c("Point_ID", "Count", "Visit_Number", "Survey_Year")]
clra22 <- c22[,c("Point_ID", "Count", "Visit_Number", "Survey_Year")]

#aggregate counts by point_ID
#Aggregate sum of % by point ID 
sumc21<- aggregate(Count ~ Point_ID+Visit_Number, clra21, FUN=sum)
sumc22<- aggregate(Count ~ Point_ID+Visit_Number, clra22, FUN=sum)

#dcast is going to "cast" the count data in a wide format by visit number  
c21.w <- dcast(sumc21, Point_ID~Visit_Number)
c22.w <- dcast(sumc22, Point_ID~Visit_Number)

#This is essentialy the format you want for every obervational covarite as well (Each site is a row and the three columns are each  visit populated with the variable of interest, here it was just the variable of counts)

#bind 
clracounts <- rbind(c21.w, c22.w)
#write.csv(clracounts,"C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/PCOuntAnalysis/Data\\CLRAc.tp.csv", row.names = TRUE)

#With coorections 
write.csv(clracounts,"C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/PCOuntAnalysis/Data\\CLRAc.b.csv", row.names = TRUE)


##SESP--------------------------------------------------------------------------
#Need to set count to 1 for species of interest and all the rest to 0
#This keeps all of the sites instead of subsetting out 
pass$Count[which(pass$Species != "SESP")] <- 0

#Need Count to be 0 for species = na & when detect = no 
pass$Count[is.na(pass$Species)] <- 0
pass$Count[which(pass$Detect == "No")] <- 0

# Truncate distances based on histogram of distance data 
pass$Count[which(pass$Distance >= 100 )] <- 0

#treat as numeric 
pass$Count <- as.numeric(pass$Count)


#aggregate counts by species per site, visit, and year 
#First subset by year 
#2021 
s21 <- subset(pass, pass$Survey_Year==2021)
#2022
s22 <- subset(pass, pass$Survey_Year==2022)

#Now select for point_id, count, visit, and year 
SESP21 <- s21[,c("Point_ID", "Count", "Visit_Number", "Survey_Year")]
SESP22 <- s22[,c("Point_ID", "Count", "Visit_Number", "Survey_Year")]

#aggregate counts by point_ID
#Aggregate sum of % by point ID 
sums21<- aggregate(Count ~ Point_ID+Visit_Number, SESP21, FUN=sum)
sums22<- aggregate(Count ~ Point_ID+Visit_Number, SESP22, FUN=sum)

#dcast is going to "cast" the count data in a wide format by visit number  
s21.w <- dcast(sums21, Point_ID~Visit_Number)
s22.w <- dcast(sums22, Point_ID~Visit_Number)

#This is essentialy the format you want for every obervational covarite as well (Each site is a row and the three columns are each  visit populated with the variable of interest, here it was just the variable of counts)


#bind 
SESPcounts <- rbind(s21.w, s22.w)
#write.csv(SESPcounts,"C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/PCOuntAnalysis/Data\\SESPc.tp.csv", row.names = TRUE)

#With corrections 
write.csv(SESPcounts,"C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/PCOuntAnalysis/Data\\SESPc.b.csv", row.names = TRUE)

##LEBI--------------------------------------------------------------------------
#Need to set count to 1 for species of interest and all the rest to 0
#This keeps all of the sites instead of subsetting out 
pass$Count[which(pass$Species != "LEBI")] <- 0

#Need Count to be 0 for species = na & when detect = no 
pass$Count[is.na(pass$Species)] <- 0
pass$Count[which(pass$Detect == "No")] <- 0

# Truncate distances based on histogram of distance data 
pass$Count[which(pass$Distance >= 175 )] <- 0

#TReat as numeric 
pass$Count <- as.numeric(pass$Count)

#aggregate counts by species per site, visit, and year 
#First subset by year 
#2021 
l21 <- subset(pass, pass$Survey_Year==2021)
#2022
l22 <- subset(pass, pass$Survey_Year==2022)

#Now select for point_id, count, visit, and year 
LEBI21 <- l21[,c("Point_ID", "Count", "Visit_Number", "Survey_Year")]
LEBI22 <- l22[,c("Point_ID", "Count", "Visit_Number", "Survey_Year")]

#aggregate counts by point_ID
#Aggregate sum of % by point ID 
suml21<- aggregate(Count ~ Point_ID+Visit_Number, LEBI21, FUN=sum)
suml22<- aggregate(Count ~ Point_ID+Visit_Number, LEBI22, FUN=sum)

#dcast is going to "cast" the count data in a wide format by visit number  
l21.w <- dcast(suml21, Point_ID~Visit_Number)
l22.w <- dcast(suml22, Point_ID~Visit_Number)

#This is essentialy the format you want for every obervational covarite as well (Each site is a row and the three columns are each  visit populated with the variable of interest, here it was just the variable of counts)


#bind 
LEBIcounts <- rbind(l21.w, l22.w)
#write.csv(LEBIcounts,"C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/PCOuntAnalysis/Data\\LEBIc.tp.csv", row.names = TRUE)

#With corrections 
write.csv(LEBIcounts,"C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/PCOuntAnalysis/Data\\LEBIc.b.csv", row.names = TRUE)



