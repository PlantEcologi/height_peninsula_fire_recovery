##############################################################################
######## Code to calculate veld age for rasters from fire polygon data
##############################################################################
######## Compiled by Jasper Slingsby 2019
######## Last edited: 15 April 2026
##############################################################################
# Data required: fire layers, a raster grid, dates of interest to calculate age at
##############################################################################

library(terra)
library(doParallel)
library(reshape2)
library(ggplot2)
library(dplyr)

### Get and subset fire data
fires <- vect("data/TMNP_fire/tmnp_fire_1962_2021.shp")

### Get a raster to use the grid (MODIS in this case)
NDVI <- rast("data/NDVI_stack_2001_2026.tif")
NAflag(NDVI) <- -9999
crs(NDVI) <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m"
NDVI <- project(NDVI, y = crs(fires)) #Make sure it's the same projection

# Load dates
NDVIdates <- read.csv("data/NDVI_dates_2001_2026.csv")

# Assign time
time(NDVI) <- as.Date(NDVIdates$date)

# Select only NDVI within the fire observation period
NDVI <- NDVI[[time(NDVI) < as.Date("2023-12-31")]]

###A fix for my specific dataset
#Set fires with no day or month to be on the 1 Jan each year
fires$STARTDATE <- paste0(substr(fires$STARTDATE,1,4), gsub("0000", "0101", substr(fires$STARTDATE,5,8)))
fires$STARTDATE <- paste0(substr(fires$STARTDATE,1,6), gsub("00", "01", substr(fires$STARTDATE,7,8)))
fires$STARTDATE <- gsub("19939315", "19930315", fires$STARTDATE)

### Crop fires to NDVI extent
fires <- crop(fires, NDVI)

#Fudge to make all "day numbers" (Days since 1970-01-01) positive
#fires$STARTDATE <- as.numeric(as.character(fires$STARTDATE))
fires[["STARTDATE"]] <- as.numeric(as.Date(as.character(fires$STARTDATE), format = "%Y%m%d")) 
#fires$STARTDATE <- type.convert(fires$STARTDATE)

###Set dates
date <- as.Date(time(NDVI), origin = "1970-01-01")
dates <- as.numeric(date)

### Run basic for loop
rfi <- list()
for (i in 1:length(dates)){ #(y=dates,.combine=c,.packages="terra") %dopar% { 
  # Loop through dates making a raster of image date - fire dates
  rfi[[i]] <- dates[i] - rasterize(fires[which(fires$STARTDATE<dates[i]),], NDVI, field="STARTDATE", fun="max", background=NA, na.rm = T) #Make background 1 Jan 1970
  # Return the individual rasterStack of "days since last fire" for each image date
} # End parallelized code

rfi <- rast(rfi)

### Parallel processing for age calculation - not working...
# # Set up cluster
# no_cores <- detectCores() - 1
# cl <- makeCluster(no_cores, type = "FORK")
# registerDoParallel(cl)

# # Start parallel processing
# rfi <- foreach(y=dates,.combine=c,.packages="terra") %dopar% { 
#   # Loop through dates making a raster of image date - fire dates
#   td <- y - terra::rasterize(fires[which(fires$STARTDATE<y),], NDVI, field="STARTDATE", fun="min", background=NA, na.rm = T) #Make background 1 Jan 1970
#   # Return the individual rasterStack of "days since last fire" for each image date
#   return(td)
# } # End parallelized code
#
# stopImplicitCluster()

names(rfi) <- date
# rfi <- setZ(rfi, date)

writeRaster(rfi, "data/veldage2023.tif", overwrite = TRUE)

### Code to process output raster stack - in this case I select only sites that burnt on record
fires$IDs <- 1:nrow(fires) # Add an unique ID for each fire
firecount <- rasterizeGeom(fires, NDVI, fun="count") # A count of #fires per pixel
firecount <- rasterize(fires, NDVI, fun="count") # A count of #fires per pixel
firecount <- firecount>0
rfi <- mask(rfi, firecount>0, maskvalue = 0)

# ### Reshape data to dataframe (if required) - messy, need to clean up!!!
# adat <- as.points(rfi)
# adatXY <- adat[,1:2]
# adat <- as.data.frame(t(adat[,3:ncol(adat)])) #transpose data.frame
# adat$Date <- as.Date(rownames(adat), format = "X%Y.%m.%d")
# adat <- melt(adat, id=c("Date"))
# adatXY$variable <- paste0("V", rownames(adatXY))
# adat <- merge(adat, adatXY)
# names(adat)[which(names(adat)=="value")] <- "Age"
# adat <- na.omit(adat)
# adat$Age <- adat$Age/365.25 #Age in years
