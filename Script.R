###SHRINKING WATER LEVELS AT THE THEEWATERSKLOOF RESERVOIR, CAPE TOWN, SA ###



####################### Background and goals ####################### 

# The region of Cape town, SA, suffered from a severe drought from 2015 to 2017.
# Those were three of the top five driest years since 1977.
# The average precipitation was only 46% of the median, in 2017 even only 30%.
# All drinking water for the city of cape town comes from a few huge reservoirs in the surroundings.
# Due to the lack of precipitation, the volume of water in the reservoirs fell down up to 10%
# The last 10% of a reservoir can not be used as drinking water.
# This led to severe Water restrictions and rationing by the government for the city of Cape Town.

# The largest reservoir is the Theewaterskloff reservoir in the mountains east of Cape town.
# It suffered from massive dehydration and the water level dropped a lot, leaving huge wastelands behind.
# In this exercise, we simply want to find out:
# how much the waterbody was shrinking and how much dried up area has been generated 
# We do this only using Landsat 8 satellite imagery



############################################################################
################## Preparation before processing the data ################## 
############################################################################


####################### Data download ######################


####################### Install and load the needed packages####################### 

#install the following packages:
install.packages("raster")
install.packages("rgdal")
install.packages("RStoolbox")
install.packages("ggplot2")

#load the packages to use them in this script:
library(raster)
library(rgdal)  
library(RStoolbox)
library(ggplot2)    

####################### set the working directory####################### 

#check your current working directory:
getwd()

#set the directory to the folder you want to work in and get & store your data:
setwd("C:/Users/Jakob/OneDrive/EAGLE M.Sc/Term 1 (Winter 2018 - 2019)/Macroecology/Data/") 
#this is important, when you don't want to change your directory for every dataset you want to load in your R workspace



############################################################################
################### Landsat Data import and preprocessing ################## 
############################################################################


#for this part of the script the data is not provided due to file size
#you can download Landsat data for example @ www.earthexplorer.usgs.com for free
#you can do this code with every landsat 8 dataset

####################### preparing the landsat scenes ####################### 

#load in the original downloaded landsat 8 scenes via adressing the metadatafile (*MTL.txt):
meta2014 <- readMeta("raw/LC08_L1TP_175083_20140103_20170427_01_T1/LC08_L1TP_175083_20140103_20170427_01_T1_MTL.txt")
meta2018 <- readMeta("raw/LC08_L1TP_175083_20180114_20180120_01_T1/LC08_L1TP_175083_20180114_20180120_01_T1_MTL.txt")

#look at the structure and content of the data:
meta2014

#We do't want to adress every single band/file while processing for every step, so we
#create a layer stack from the bands to make processing faster:

stack2014 <- stackMeta(meta2014)
stack2018 <- stackMeta(meta2018)

#look at the data again, the structure has changed:
stack2014

####################### Calibration ####################### 

#The raster values are stored in digital numbers to reduce file size. 
#To recalculate the actual radiation measured at the sensor, we need to apply sensor- and band-specific parameters
#Sensor Calibration is used for converting the digital numbers to meaningful units (reflectance, radiation, temperature)
#in this case, we are converting the dns of the multispectral bands to apparent reflectance ("apref") and the thermal bands to surface temperature

#calibration:
stack2014_cal <- radCor(stack2014, metaData=meta2014, method = "apref") 
stack2018_cal <- radCor(stack2018, metaData=meta2014, method = "apref") 

#check that the data type of the bands has changed from integer to float:
dataType(stack2014)
dataType(stack2014_cal)

#new range of values also can be seen when comparing:
stack2014
stack2014_cal

#export the calibrated layerstack to the disk as a *.tif file:
writeRaster(stack2014_cal, filename="stack2014_cal", format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(stack2018_cal, filename="stack2018_cal", format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))



############################################################################
################### Rezising Data to the study area ################## 
############################################################################


####################### Loading study area ####################### 

#our area of interest is the Theewaterskloof reservoir, but the landsat scene is 170 km x 185 km
#to save disk space and computation time, we crop the data by a predefined shapefile of the study area 

#load the study area (stored as a shapefile):
Theewaterskloof_reservoir <- readOGR("Theewaterskloof_reservoir.shp")

#reproject the coordinate system of the shapefile
#it's required for further analysis that the shapefile has the same coordinate system as the layerstack:
Theewaterskloof_reservoir <- spTransform(Theewaterskloof_reservoir, CRS(proj4string(stack2014_cal)))

#plot the study area on top of the Landsat image to see if the reprojection has worked:
plotRGB(stack2014_cal, r="B4_tre", g="B3_tre", b="B2_tre", stretch="lin")
plot(Theewaterskloof_reservoir, col="yellow", add=TRUE)


####################### Cropping ####################### 

#resize the data using the crop function:
Theewaterskloof_2014 <- crop(stack2014_cal, Theewaterskloof_reservoir)
Theewaterskloof_2018 <- crop(stack2018_cal, Theewaterskloof_reservoir)

#export the data as a multilayer *.tif file to disk:
writeRaster(Theewaterskloof_2014, filename="Theewaterskloof_2014", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(Theewaterskloof_2018, filename="Theewaterskloof_2018", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))

#reload the resized data as a raster brick (raster brick save ram while processing):
Theewaterskloof_2014 <- brick("Theewaterskloof_2014.tif")
Theewaterskloof_2018 <- brick("Theewaterskloof_2018.tif")

#have a look at the data:
Theewaterskloof_2014
plot(Theewaterskloof_2014) #plots all bands seperately
ggRGB(Theewaterskloof_2014, stretch="lin") #plots a rgb image from selected bands


####################### Renaming bands ####################### 

#as you can see, the bands are named not very suitable for further analysis:
names(Theewaterskloof_2014)

#That's why we rename them to their real names using a list with the names (band names can be easily found out using Google):
Landsat8_band_names <- c("Coastal Aerosol", "Blue", "Green", "Red", "NIR", "SWIR 1", "SWIR 2", "Cirrus", "TIR1", "TIR2")
names(Theewaterskloof_2014) <- Landsat8_band_names
names(Theewaterskloof_2018) <- Landsat8_band_names

#let's have a look at the data again to see the changes:
names(Theewaterskloof_2014)

plot(Theewaterskloof_2014) #plots all bands
plot(Theewaterskloof_2018)


####################### Get a Overview of the study area ####################### 

#compare the landsat imagery from both aquisition dates:
par(mfrow=c(1,2)) #split the data viewer in 2 columns

plotRGB(Theewaterskloof_2014, r="Red", g="Green", b="Blue", stretch="lin")
legend("top", legend = NA, title = "2014", bty = "n", cex = 2) #Adds title

plotRGB(Theewaterskloof_2018, r="Red", g="Green", b="Blue", stretch="lin")
legend("top", legend = NA, title = "2018", bty = "n", cex = 2)
scalebar(5000, xy=click(), type='bar', divs=4, below="Meters", cex= 0.5) #adds scalebar

#we can clearly see change of the water body



############################################################################
############################### NDWI analysis ##############################
############################################################################


#to detect areas covered with water from remote sensing data, the NDWI indexcan be used
#NDWI = Normalized Difference Water Index

################### Computation of NDWI and MNDWI ################### 

#there are two versions of the NDWI:
#the "normal" version, which is defined as (Green-NIR/Green+NIR)
#the "modified" version, which is defined as (Green-SWIR1/Green+SWIR1)
#we have to evaluate, which version fits our area the best

#compute the NDWI:

NDWI_2014 <- (Theewaterskloof_2014[["Green"]]-Theewaterskloof_2014[["NIR"]])/(Theewaterskloof_2014[["Green"]]+Theewaterskloof_2014[["NIR"]])
NDWI_2018 <- (Theewaterskloof_2018[["Green"]]-Theewaterskloof_2018[["NIR"]])/(Theewaterskloof_2018[["Green"]]+Theewaterskloof_2018[["NIR"]])

#compute the MDWI:

MNDWI_2014 <- (Theewaterskloof_2014[["Green"]]-Theewaterskloof_2014[["SWIR.1"]])/(Theewaterskloof_2014[["Green"]]+Theewaterskloof_2014[["SWIR.1"]])
MNDWI_2018 <- (Theewaterskloof_2018[["Green"]]-Theewaterskloof_2018[["SWIR.1"]])/(Theewaterskloof_2018[["Green"]]+Theewaterskloof_2018[["SWIR.1"]])


#compare the indices by plotting them all together:
indices <- stack(NDWI_2014, NDWI_2018, MNDWI_2014, MNDWI_2018)
names(indices) <- c("NDWI 2014", "NDWI 2018", "MNDWI 2014", "MNDWI 2018")

plot(indices)

#--> the MNDWI values show the waterbodies better (the background (=soil) has more differentiable values)

#export all indices to disk, but as seperate *.tif files (named after band):
writeRaster(indices, ".tif", bylayer = T, overwrite=TRUE)



################### Calculate waterbody polygons ################### 

#at the end, we want to have shapefiles of the waterbody of 2014 and 2018 as well as of the dried up area
#to calculate the extent of the waterbodies in 2014 and 2018 we need to define thresholds: hat values are representing water? --> >=0 

#therefore we let all other values be NA (so they will not be recognised, when we transfer the raster into a shapefile):
MNDWI_2014[MNDWI_2014 < 0] <- NA
MNDWI_2018[MNDWI_2018 < 0] <- NA

#all other values need to be the same, otherwise every different value will be transferred into a single polygon:
MNDWI_2014[MNDWI_2014 >= 0] <- 1
MNDWI_2018[MNDWI_2018 >= 0] <- 1

#have a look at the new values:
dev.off()
plot(MNDWI_2014, col="blue")

#generate shapefiles (dissolve=TRUE for joining all similar values to one polygon):
Waterbody_2014 <- rasterToPolygons(MNDWI_2014, fun=function(MNDWI_2014){MNDWI_2014>0}, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
Waterbody_2018 <- rasterToPolygons(MNDWI_2018, fun=function(MNDWI_2018){MNDWI_2018>0}, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)

#export the shapefiles to disk:
writeOGR(Waterbody_2014, ".", "Waterbody_2014", driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR(Waterbody_2018, ".", "Waterbody_2018", driver="ESRI Shapefile", overwrite_layer=TRUE)

#have a look at the generated shapefiles:
par(mfrow=c(1,2))
plot(Waterbody_2014, col="blue", border=FALSE, main="Waterbody 2014")
plot(Waterbody_2018, col="blue", border=FALSE, main="Waterbody 2018")

ggRGB(Theewaterskloof_2014, r="Red", g="Green", b="Blue", stretch="lin")+
  geom_polygon(data = Waterbody_2014, aes(x = long, y = lat, group = group), fill = "blue", size = 0.25)+
  ggtitle("Waterbody 2014")

ggRGB(Theewaterskloof_2018, r="Red", g="Green", b="Blue", stretch="lin")+
  geom_polygon(data = Waterbody_2018, aes(x = long, y = lat, group = group), fill = "blue", size = 0.25)+
  ggtitle("Waterbody 2018")


#have a look the attributes of the shapefile:
Waterbody_2018
summary(Waterbody_2014)

#as you can see, we have the polygon, but no information on the area covered


################### Calculate water area ################### 

#calculate the area covered with water:

area(Waterbody_2014) #in m?
area(Waterbody_2014)/1000000 #in km?

area(Waterbody_2018)
area(Waterbody_2018)/1000000

#add the area information as a column to the shapefile table:
Waterbody_2014$area_sqm <- area(Waterbody_2014)
Waterbody_2014$area_sqkm <- area(Waterbody_2014)/1000000

Waterbody_2018$area_sqm <- area(Waterbody_2018)
Waterbody_2018$area_sqkm <- area(Waterbody_2018)/1000000

Waterbody_2014

################### Dried up Area Detection ################### 

#generate the dried up area:
Dried_up_area <- erase(Waterbody_2014, Waterbody_2018)

writeOGR(Dried_up_area, ".", "Dried_up_area", driver="ESRI Shapefile", overwrite_layer=TRUE)

dev.off()

#plot the dried up area on top of the landsat imagery:
plotRGB(Theewaterskloof_2018, r="Red", g="Green", b="Blue", stretch="lin")
plot(Dried_up_area, col="orange", border=FALSE, add=TRUE)

#calculate the dried up area:

area(Dried_up_area) #in m²
area(Dried_up_area)/1000000 #in km²

#add it to the shapefile as a column:
Dried_up_area$area_sqm <- area(Dried_up_area)
Dried_up_area$area_sqkm <- area(Dried_up_area)/1000000

#THE END