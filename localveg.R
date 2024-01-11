###Local veg variables

setwd("C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/")

###Overall Objectives: ------------------------------------------------------------
##1. Percent cover of SPAL, JURO, and freshwater plants at each point for both years
##2. Percent cover (density) of all spp at each point for both years
##3. Average height of all spp at each point for both years
##4. Percent cover of open water at each point for both years 

##1----------------------------------------------------------------------------------
##Percent cover of JURO 
#Import domspp dataset
domspp <- read.csv("C:/Users/Rachel Anderson/Desktop/Analysis/Data/domspp.clean.csv")
str(domspp) 
names(domspp)
View(domspp)

#delete first column and extra rows (found later; data was entered twice for GBR01_2 & HCR06_2)
domspp<- domspp[,-1]

##Percent cover of JURO at each point for both years 
#First, need all sites with zero for percentage if not JURO 
domspp$DS_Percentage[which(domspp$Species != "JURO")] <- "0"
domspp$DS_Percentage <- as.numeric(domspp$DS_Percentage)

#Now convert all to JURO 
domspp$Species[which(domspp$Species != "JURO")] <- "JURO"

#2021 
domj21 <- subset(domspp, domspp$Year==21)

#2022
domj22 <- subset(domspp, domspp$Year==22)

#Aggregate sum of % by point ID 
sumj21 <- aggregate(DS_Percentage ~ Point_ID, domj21, FUN=sum)

#Make sure no errors
over21 <- subset(sumj21, sumj21$DS_Percentage > 100)

#Do the same for 2022
sumj22 <- aggregate(DS_Percentage ~ Point_ID, domj22, FUN=sum)

#Make sure no errors
over22 <- subset(sumj22, sumj22$DS_Percentage > 100)
#n = 266 because NRI was surveyed that year 

#combine and export 
j2122 <- rbind(sumj21, sumj22[-c(151,152,153,154,155),]) #don't need NRI for now

write.csv(j2122,"C:/Users/Rachel Anderson/Desktop/Analysis/Data\\JUROper.csv", row.names = TRUE)

#---------------------------------------------------------------------------------
##Do the same procedure for % SPAL
#Import domspp dataset
domspp <- read.csv("C:/Users/Rachel Anderson/Desktop/Analysis/Data/domspp.clean.csv")

#delete first column and extra rows (found later; data was entered twice for GBR01_2 & HCR06_2)
domspp<- domspp[,-1]

##Percent cover of SPAL at each point for both years 
#First, need all sites with zero for percentage if not SPAL 
domspp$DS_Percentage[which(domspp$Species != "SPAL")] <- "0"
domspp$DS_Percentage <- as.numeric(domspp$DS_Percentage)

#Now convert all to SPAL
domspp$Species[which(domspp$Species != "SPAL")] <- "SPAL"

#2021 
doms21 <- subset(domspp, domspp$Year==21)

#2022
doms22 <- subset(domspp, domspp$Year==22)

#Aggregate sum of % by point ID 
sums21 <- aggregate(DS_Percentage ~ Point_ID, doms21, FUN=sum)

#Make sure no errors
overs21 <- subset(sums21, sums21$DS_Percentage > 100)

#Do the same for 2022
sums22 <- aggregate(DS_Percentage ~ Point_ID, doms22, FUN=sum)
#n = 266 because NRI was surveyed that year 

#Make sure no errors
overs22 <- subset(sums22, sums22$DS_Percentage > 100)

#combine and export 
s2122 <- rbind(sums21, sums22[-c(151,152,153,154,155),]) #don't need NRI for now

write.csv(s2122,"C:/Users/Rachel Anderson/Desktop/Analysis/Data\\SPALper.csv", row.names = TRUE)

#-------------------------------------------------------------------------------
##Percent fresh 
#Need to decide which plants are indicative of freshwater 
#SALA, SALT, POCO, CLMA, THPA, ALPH

#Import domspp dataset
domspp <- read.csv("C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/Data/domspp.clean.csv")

#delete first column and extra rows (found later; data was entered twice for GBR01_2 & HCR06_2)
domspp<- domspp[,-1]

##Percent cover of freshwater spp at each point for both years 
#Need to set this up a little differently since there are multiple species   
domspp$Species[which(domspp$Species == "SALA")] <- "FRESH"
domspp$Species[which(domspp$Species == "SALT")] <- "FRESH"
domspp$Species[which(domspp$Species == "POCO")] <- "FRESH"
domspp$Species[which(domspp$Species == "CLMA")] <- "FRESH"
domspp$Species[which(domspp$Species == "THPA")] <- "FRESH"
domspp$Species[which(domspp$Species == "ALPH")] <- "FRESH"

#Fill in zeros for non FRESH
domspp$DS_Percentage[which(domspp$Species != "FRESH")] <- "0"
domspp$DS_Percentage <- as.numeric(domspp$DS_Percentage)

#Now convert all to FRESH
domspp$Species[which(domspp$Species != "FRESH")] <- "FRESH"

#2021 
domf21 <- subset(domspp, domspp$Year==21)

#2022
domf22 <- subset(domspp, domspp$Year==22)

#Aggregate sum of % by point ID 
sumf21 <- aggregate(DS_Percentage ~ Point_ID, domf21, FUN=sum)

#Make sure no errors
overf21 <- subset(sumf21, sumf21$DS_Percentage > 100)

#Do the same for 2022
sumf22 <- aggregate(DS_Percentage ~ Point_ID, domf22, FUN=sum)
#n = 266 because NRI was surveyed that year 

#Make sure no errors
overf22 <- subset(sumf22, sumf22$DS_Percentage > 100)

#combine and export 
f2122 <- rbind(sumf21, sumf22[-c(151,152,153,154,155),]) #don't need NRI for now

write.csv(f2122,"C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/Data\\FRESHper.csv", row.names = TRUE)

##2-------------------------------------------------------------------------------
#Percent cover of all spp for both years 

#Import domspp dataset
domspp <- read.csv("C:/Users/Rachel Anderson/Desktop/Analysis7.10.23/Data/domspp.clean.csv")

#delete first column and extra rows (found later; data was entered twice for GBR01_2 & HCR06_2)
domspp<- domspp[,-1]
domspp<- domspp[,-2]#don't need date either


#2021 
domspp21 <- subset(domspp, domspp$Year==21)
length(unique(domspp21$Point_ID))

#2022
domspp22 <- subset(domspp, domspp$Year==22)
length(unique(domspp22$Point_ID)) #this includes NRI

#Aggregate sum of % by point ID 
sumspp21 <- aggregate(DS_Percentage ~ Point_ID, domspp21, FUN=sum)
sumspp22 <- aggregate(DS_Percentage ~ Point_ID, domspp22, FUN=sum)

#Check min and max 
min(sumspp21$DS_Percentage)
max(sumspp21$DS_Percentage)

min(sumspp22$DS_Percentage)
max(sumspp22$DS_Percentage)

which(sumspp21$DS_Percentage > 100)
gb <- subset(domspp, domspp$Point_ID=="GBR01_1")



##3-------------------------------------------------------------------------------
#Average height of all spp at each point for both years

#Import height dataset 
height <- read.csv("C:/Users/Rachel Anderson/Desktop/Analysis/Data/Heights.clean.csv")

#delete first column and extra rows (found later; data was entered twice for GBR01_2 (2021) & HCR06_2 (2022))
height<- height[,-1]
height<- height[-c(1469, 1470,1471,1472,1473,1474,1475,1476,1478,1479,1480,1481,1482,1483, 2795, 2796,2797,2798,2799,2800,2801,2802,2803,2804,2805),]

#2021 
hspp21 <- subset(height, height$Year==21)
length(unique(hspp21$Point_ID))

#2022
hspp22 <- subset(height, height$Year==22)
length(unique(hspp22$Point_ID))

#Aggregate avg height of all spp by point ID 
avgspp21 <- aggregate(Height ~ Point_ID, hspp21, FUN=mean)
avgspp22 <- aggregate(Height ~ Point_ID, hspp22, FUN=mean)

#Check min and max 
min(avgspp21$Height)
max(avgspp21$Height)

min(avgspp22$Height)
max(avgspp22$Height)


#combine and export 
avgspp <- rbind(avgspp21, avgspp22)
write.csv(avgspp,"C:/Users/Rachel Anderson/Desktop/Analysis/Data\\avgspp.csv", row.names = TRUE)


##4-------------------------------------------------------------------------------
#Percent open water 
comm <- read.csv("C:/Users/Rachel Anderson/Desktop/Analysis/Data/commper_clean.csv")

#delete first column and extra rows (found later; data was entered twice for GBR01_2 (2021) & HCR06_2 (2022))
comm<- comm[,-1]
comm<- comm[-c(497,498,1280)] 

##Percent cover of open water at each point for both years 
#First, need all sites with zero for percentage if not SPAL 
comm$Hab_Percentage[which(comm$CommHab != "OpenWater")] <- "0"
comm$Hab_Percentage <- as.numeric(comm$Hab_Percentage)

#Now convert all to SPAL
comm$CommHab[which(comm$CommHab != "OpenWater")] <- "OpenWater"

#2021 
comm21 <- subset(comm, comm$Year==21)
length(unique(comm21$Point_ID))

#2022
comm22 <- subset(comm,comm$Year==22)
length(unique(comm22$Point_ID))

min(comm21$Hab_Percentage)
max(comm21$Hab_Percentage)


min(comm22$Hab_Percentage)
max(comm22$Hab_Percentage)

#combine and export 
openwater <- rbind(comm21, comm22[-c(1349,1350,1351,1352, 1363, 1364, 1365),]) #don't need NRI for now
write.csv(openwater,"C:/Users/Rachel Anderson/Desktop/Analysis/Data\\openwater.csv", row.names = TRUE )
