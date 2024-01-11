#install.packages("readxl")
library(readxl)

#Set your working directory 
setwd("C:/Users/Rachel Anderson/Desktop/BMB Analysis/Data")

#bring in your total dataset and take a look 
total <- read_excel(file.choose())
str(total) 
head(total)

#fix dates to treat as a character instead of a numeric value 
total$Survey_Date <- as.character(total$Survey_Date)

#calculate Julian Day for each date 
total$Day <- as.POSIXlt(total$Survey_Date)$yday + 1 


#select for species of interest. In this case, I pulled out all of my clapper rail detections. You can change "CLRA" to whatever species you want to do the analysis for
CLRA <- total[total$Species== "22",] 
#22 was the ID in my database, so I changed it to "CLRA" instead
#change 22 
CLRA$Species[CLRA$Species=="22"] <- "CLRA"
#write.csv(CLRA,"C:/Users/Rachel Anderson/Desktop/BMB Analysis/Data\\CLRA.csv", row.names = TRUE)

#Create data frame with the columns of interest 
#SurveyPCountObs_Point_ID is the site ID 
#F_Observer is "Focal observer', which was also a numeric ID in my database 
CLRA <- CLRA[,c("Species", "Distance", "Day", "SurveyPCountObs_Point_ID", "Visit_Number","Survey_Year","F_Observer")]

#treat Observer as a factor 
newyear$F_Observer<- as.factor("F_Observer")

#subset to only include data from 2021
#You can choose whatever years you want 
newyear <- subset(CLRA, Survey_Year=="2021")
#so you now have a datframe with the columns of interest for a particular species and year(s)

###Now you want to reshape the format from long to wide for unmarked-----------------------
library(reshape2)

#dcast is going to "cast" the count data in a wide format by visit number  
CLRA.w <- dcast(newyear, SurveyPCountObs_Point_ID~Species+Visit_Number)
#This is essentialy the format you want for every obervational covarite as well (Each site is a row and the three columns are each  visit populated with the variable of interest, here it was just the variable of counts)


###reshape julian day from long to wide 
#dates can get a little weird so we had to add some extra code to get it formatted correctly 
duplicate <-newyear[!duplicated(newyear[,c("SurveyPCountObs_Point_ID","Visit_Number")]),]
#duplicated makes sure you're taking dates that were duplicates (different point counts that were visited on the same day due to different teams surveying)

ram <- melt(duplicate[,c(3,4,5)],id=2:3) #melt the data frame by the value of day  

CLRA.t<- reshape(data = ram[,-3], idvar="SurveyPCountObs_Point_ID", v.names = "value",timevar = "Visit_Number", direction = "wide")
View(CLRA.t) #Reshape takes the value of day and populates by visit number into a wide format 
#export dataframe to import into unmarked
#write.csv(CLRA.t,"C:/Users/Rachel Anderson/Desktop/BMB Analysis/Data\\CLRA.time.csv", row.names = TRUE)

######################################Reshape new dates 
year <- subset(total, Survey_Year=="2021")
#write.csv(year,"C:/Users/Rachel Anderson/Desktop/BMB Analysis/Data\\year.csv", row.names = TRUE)
yearcounts <- year[,c("Species","Day", "SurveyPCountObs_Point_ID", "Visit_Number")]

ram3 <- melt(yearcounts[,c(2,3,4)],id=2:3)

year.w<- reshape(data = ram3[,-3], idvar="SurveyPCountObs_Point_ID", v.names = "value",timevar = "Visit_Number", direction = "wide")
View(year.w)
#write.csv(year.w,"C:/Users/Rachel Anderson/Desktop/BMB Analysis/Data\\year.w.csv", row.names = TRUE)

#########################################Reshape Wind & sky covariates 
Wns <- year[,c("Species","WindB","SurveyPCountObs_Point_ID", "Visit_Number")]
ram4 <- melt(Wns[,c(2,3,4)],id=2:3)
WnS.w<- reshape(data =ram4[,-3], idvar="SurveyPCountObs_Point_ID", v.names = "value",timevar = "Visit_Number", direction = "wide") 
View(WnS.w)
write.csv(WnS.w,"C:/Users/Rachel Anderson/Desktop/BMB Analysis/Data\\WnS.w.csv", row.names = TRUE)


Sns <- year[,c("Species","Sky","SurveyPCountObs_Point_ID", "Visit_Number")]
ram5 <- melt(Sns[,c(2,3,4)],id=2:3)
SnS.w<- reshape(data =ram5[,-3], idvar="SurveyPCountObs_Point_ID", v.names = "value",timevar = "Visit_Number", direction = "wide") 

#write.csv(SnS.w,"C:/Users/Rachel Anderson/Desktop/BMB Analysis/Data\\SnS.w.csv", row.names = TRUE)

#########################################Reshape Observer
Obs <- year[,c("Species","F_Observer","SurveyPCountObs_Point_ID", "Visit_Number")]
ram6 <- melt(Obs[,c(2,3,4)],id=2:3)

CLRA.O<- reshape(data = ram6[,-3], idvar="SurveyPCountObs_Point_ID", v.names = "value",timevar = "Visit_Number", direction = "wide")

write.csv(CLRA.t,"C:/Users/Rachel Anderson/Desktop/BMB Analysis/Data\\CLRA.time.csv", row.names = TRUE)




