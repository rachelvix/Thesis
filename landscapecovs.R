####Landscape predictor analysis 
###Moving window for landscape predictors 
# 1. Proportion of estuarine emergent marsh in 1 km buffer
# 2. Proportion of palustrine emergent marsh in 1 km buffer

##Note: I found it easier to import raster into ArcGIS to visualize map and further extract point values; 
#This script is for projecting the 10 m ccap data and conducting the moving window analysis only 

#----------------------------------------------------------------------------------------
#Load packages 
library(terra)
#library(dismo)
library(sf)
#library(maptools)
#library(rgeos)
#library(dplyr)
library(tidyr)
library(ggplot2)

#####Bring in landcover data 
ccap <- ("C:/Users/Rachel Anderson//Desktop/MarshGIS/ccap10m/ccap10m.tif")

#Import landcover data file as a raster file 
ccap.r <- rast(ccap)
ccap.r #NAD83

# Plot:
plot(ccap.r)

#Check how many classes do you have
unique(ccap.r) 
#15 classes

#####Read in SHP file for extent and crop 
shape <- st_read("C:/Users/Rachel Anderson//Desktop/MarshGIS/ccap10m/marsh_project_sa.shp")
#Shape file is in mstm, so we need to reproject 

#treat as vector 
shape <- vect(shape)

#project shape file crs of ccap
repro <- project(shape, crs(ccap.r))

#Check to see if this worked
repro

#set ccap extent to match shape file 
extent <- as.vector(ext(repro))
ccap.ext <- crop(ccap.r, extent)

#Check this too 
ext(repro)
ext(ccap.ext)
#This doesn't match exactly but close...Idk if that's a problem


#now plot again 
plot(ccap.ext)



#create a vector of landcover values for area of interest:
unique(ccap.ext)
ccap_class <- c(0,2,8,11,12,13,14,15,17,18,19,20,21,22,23)


#------------------------------------------Prop Est Marsh--------------------------------------------
# create a vector specifying how you want to reclassify each cover type. 
# this vector should be the same length as your above vector  
ccap_est <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)

#Combine vectors into a matrix 
reclass.mat <- cbind(ccap_class, ccap_est)

#View it to make sure that you assigned the reclassification categories correctly
# double check the classification legend. 
View(reclass.mat)
#looks good

#Now, reclassify your raster file based on this reclassification matrix
ccap_reclass <- classify(ccap.ext, reclass.mat)

#plot 
plot(ccap_reclass)

#Make sure the values are correct
unique(ccap_reclass)

# Now that we have the reclassified raster for our landcover type of 
# interest, we want to do the moving window analysis for each 10 m raster cell. 
# An important thing to consider is how wide of a buffer you are interested 
# in / how wide of a buffer will impact your species of interest. 
# Some people choose a few different buffer radii and look for the ones 
# that are most closely correlated with their species (you would use a 
# pairplot like the onle Carlos sent you).  
# I might recommend trying 200 m, 500 m, and 1000 m at least. 
# Or base it on something known about the species already. 

# set buffer radius of interest:
buffer.radius <- (1000)

# We are going to create a matrix of focal weights to use as the 
# moving window.  I am using a circle here, so it will be a buffer
# that is a circle with a 100m radius.  
fm.1000m <- focalMat(ccap_reclass, buffer.radius, type = "circle")
fm.1000m 

#Convert nonzero values to 1 - so the matrix becomes binary
fm.1000m <- ifelse(fm.1000m > 0, 1, 0)
fm.1000m # Now you can see they are all 0s and 1s


#### Conduct the moving window analysis  
# w is your matrix of weights, fun = "sum" means that the focal 
# function will add up all of the values within that moving window. 
# since we made the raster binary (developed = 1, nondeveloped = 0), 
# the sum will represent the number of cells within our buffer radius 
# that are developed.  
est.1000m <- focal(ccap_reclass, w = fm.1000m, fun = "sum")
#View a plot of it: 
plot(est.1000m)

# Now change the sum of the other cells into a proportion of "estuarine" cells: 
est.prop.1000m <- est.1000m / sum(fm.1000m)*100
#Check the plot - this will look exactly the same as the last plot, but 
# note that the scale is now a proportion instead of a number of cells. 
plot(est.prop.1000m)
plot(ccap_reclass)

# You're done! This raster represents the proportion of land cover area
# within a 1000 m radius that is "estuarine marsh". 

#save/export the raster file - then you can take it into GIS or whatever other
# program you are using for adjusting the resolution/extent of your 
# raster files. 

f <- file.path("C:/Users/Rachel Anderson//Desktop/MarshGIS/ccap10m", "est.prop1000.tif")
writeRaster(est.prop.1000m, f)



#----------------------------------------------------Prop palustrine marsh----------------------

#------------------------------------------Prop Est Marsh--------------------------------------------
# create a vector specifying how you want to reclassify each cover type. 
# this vector should be the same length as your above vector  
ccap_pal <- c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)
#ccap_class<- c(0,2,8,11,12,13,14,15,17,18,19,20,21,22,23)

#Combine vectors into a matrix 
reclass.matp <- cbind(ccap_class, ccap_pal)

#View it to make sure that you assigned the reclassification categories correctly
# double check the classification legend. 
View(reclass.matp)
#looks good

#Now, reclassify your raster file based on this reclassification matrix
ccap_reclassp <- classify(ccap.ext, reclass.matp)

#plot 
plot(ccap_reclassp)

#Make sure the values are correct
unique(ccap_reclassp)

# Now that we have the reclassified raster for our landcover type of 
# interest, we want to do the moving window analysis for each 10 m raster cell. 
# An important thing to consider is how wide of a buffer you are interested 
# in / how wide of a buffer will impact your species of interest. 
# Some people choose a few different buffer radii and look for the ones 
# that are most closely correlated with their species (you would use a 
# pairplot like the onle Carlos sent you).  
# I might recommend trying 200 m, 500 m, and 1000 m at least. 
# Or base it on something known about the species already. 

# set buffer radius of interest:
buffer.radius <- (1000)

# We are going to create a matrix of focal weights to use as the 
# moving window.  I am using a circle here, so it will be a buffer
# that is a circle with a 100m radius.  
fm.1000m.p <- focalMat(ccap_reclassp, buffer.radius, type = "circle")
fm.1000m.p 

#Convert nonzero values to 1 - so the matrix becomes binary
fm.1000m.p <- ifelse(fm.1000m.p > 0, 1, 0)
fm.1000m.p # Now you can see they are all 0s and 1s


#### Conduct the moving window analysis  
# w is your matrix of weights, fun = "sum" means that the focal 
# function will add up all of the values within that moving window. 
# since we made the raster binary (developed = 1, nondeveloped = 0), 
# the sum will represent the number of cells within our buffer radius 
# that are developed.  
pal.1000m <- focal(ccap_reclassp, w = fm.1000m.p, fun = "sum")
#View a plot of it: 
plot(pal.1000m)

# Now change the sum of the other cells into a proportion of "estuarine" cells: 
pal.prop.1000m <- pal.1000m / sum(fm.1000m.p)*1000
#Check the plot - this will look exactly the same as the last plot, but 
# note that the scale is now a proportion instead of a number of cells. 
plot(pal.prop.1000m)
plot(ccap_reclassp)

# You're done! This raster represents the proportion of land cover area
# within a 1000 m radius that is "estuarine marsh". 

#save/export the raster file - then you can take it into GIS or whatever other
# program you are using for adjusting the resolution/extent of your 
# raster files. 

f <- file.path("C:/Users/Rachel Anderson//Desktop/MarshGIS/ccap10m", "pal.prop1000.tif")
writeRaster(pal.prop.1000m, f)

#------------------------------------------------------------------------------
###Shannon's diversity index
setwd("C:/Users/Rachel Anderson/Desktop/MarshGIS/ccap10m")
library(landscapemetrics)
library(terra)

#Load data
ccap.r <- ("C:/Users/Rachel Anderson/Desktop/MarshGIS/ccap10m/ccap_crop.tif")
ccap.r <- rast(ccap.r)
plot(ccap.r)

##Add points to map  
pots <- readOGR(".","stupid")

#check projections match  
ccap.r
pots  #Extents don't match       

ext(pots)
ext(ccap.r)

ccap.r.e <- as.vector(ext(ccap.r))
ext.points <- crop(pots, ccap.r.e)

ccap.r.e <- as.vector(ext(pots))
ext.points <- crop(ccap.r, ccap.r.e)

repro <- project(pots, crs(ccap.r))

pots$bbox

################################################################################
###1. Edge density in 200 meter buffer 
##Utilize landscapemetrics package 
#look at class metrics

lsm.landscape <- list_lsm(level="landscape")   #landscape-level metrics
View(lsm.landscape) #shdi

#calculate ed
shdi.200 <- sample_lsm(ccap.r, y = pots, plot_id = pots$points, shape="circle", size=200, 
                     what =c("lsm_l_shdi"), all_classes=TRUE, 
                     return_raster=FALSE, verbose=TRUE, progress=TRUE)

shdi

shdi200.df<-as.data.frame(shdi.200)
View(shdi200.df)
write.csv(shdi200.df, "C:/Users/Rachel Anderson/Desktop/MarshGIS/ccap10m/shdi200.csv")


# set buffer radius of interest:
buffer.radius <- (200)

# We are going to create a matrix of focal weights to use as the 
# moving window.  I am using a circle here, so it will be a buffer
# that is a circle with a 100m radius.  
fm.200m.w <- focalMat(ccap.r, buffer.radius, type = "circle")
fm.200m.w 

#Convert nonzero values to 1 - so the matrix becomes binary
fm.200m.w <- ifelse(fm.200m.w > 0, 1, 0)
fm.200m.w # Now you can see they are all 0s and 1s

shdi200.window <- window_lsm(ccap.r, window=fm.200m.w, what=c("lsm_l_shdi"), level="landscape", progress=TRUE)

# Now change the sum of forest cells into a proportion of developed cells: 
shdi.prop.200m <- shdi200.window / sum(fm.200m.w)*100
#Check the plot - this will look exactly the same as the last plot, but 
# note that the scale is now a proportion instead of a number of cells. 
plot(shdi.prop.200m)

#save/export the raster file - then you can take it into GIS or whatever other
# program you are using for adjusting the resolution/extent of your 
# raster files.

f <- file.path("C:/Users/lab/Desktop/Rachel/MarshGIS/", "shdi.prop200.tif")

writeRaster(shdi.prop.200m, f) 




###=====================================================
#small crop of study area for landon

#import ccap reclass file 
ccap <- ("C:/Users/Rachel Anderson//Desktop/MarshGIS/ccap10m/ccap_reclass.tif")

#Import landcover data file as a raster file 
ccap.r <- rast(ccap)
ccap.r #NAD83

#Read shape file
shape <- st_read("C:/Users/Rachel Anderson//Desktop/MarshGIS/din/din.shp")

#project shape file crs of ccap
repro <- project(shape, crs(ccap.r))

#Check to see if this worked
repro

#crop study area to shape file extent 
#set ccap extent to match shape file 
extent <- as.vector(ext(repro))
ccap.ext <- crop(ccap.r, extent)

#Check this too 
ext(repro)
ext(ccap.ext)
#This doesn't match exactly but close...Idk if that's a problem

#now plot again 
plot(ccap.ext)

#export 
f <- file.path("C:/Users/Rachel Anderson//Desktop/MarshGIS/ccap10m", "din_crop.tif")
writeRaster(ccap.ext, f) 


#import ccap reclass file 
din <- ("C:/Users/Rachel Anderson//Desktop/MarshGIS/ccap10m/din_crop.tif")
#Import landcover data file as a raster file 
din.r <- rast(din)
din.r #NAD83

plot(din.r)
