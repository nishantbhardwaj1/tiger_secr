install.packages(c("camtrapR", "tidyr","raster","secr"))
library(camtrapR)
library(raster)
library(secr)
library(tidyr)
#exiftool is required for reading data from camera trap images 
# for more info read https://jniedballa.github.io/camtrapR/ 
exiftool_dir <- "E:/Software/EXIF"


## Read and plot the shapefile 


#Reading the shapefile for the wildlife sanctuary
NWLS<-shapefile('nwls_mask_2024.shp')
plot(NWLS)



## Reading the Camera operation file to create a matrix
NWLS_camop<-read.csv("GPS_2024_utm44.csv")
head(NWLS_camop)
cameraOperation

#creating on-off matrix/trap matrix for camera trap stations to understand effort 
NWLS_camop_table<-cameraOperation(NWLS_camop,
                                  stationCol = "station",
                                  setupCol = "Start_Date",
                                  retrievalCol = "end_date",
                                  dateFormat = "%d-%m-%Y",
                                  writecsv= TRUE,
                                  outDir= "H:/NWLS_data_2024/Tiger2024/Output")

head(NWLS_camop_table)

# this command is used for creating a capture matrix for the tigers
# It is not used here as directory 'experiment' which has tiger images cannot 
# be shared without permission
{tiger_record_table<-recordTableIndividual(inDir = "experiment",
                                           hasStationFolders = TRUE,           
                                           IDfrom = "directory",
                                           camerasIndependent = TRUE,
                                           minDeltaTime = 01,
                                           deltaTimeComparedTo = "lastIndependentRecord",
                                           timeZone = "Asia/Calcutta",
                                           stationCol = "station",
                                           writecsv = TRUE,
   
                                                                                   outDir = "H:/NWLS_data_2024/Tiger2024/Output")}
# output for the above command is provided to allow analysis i.e capture matrix
tiger_record_table<-read.csv("Output/record_table_individuals1min_deltaT_2024-10-15.csv")
#

# creating a capthist object required for secr
# need trap matrix and capture matrix
tiger_capthist1 <-spatialDetectionHistory(tiger_record_table,
                                           species = "experiment",
                                           output = "binary",
                                           camOp = NWLS_camop_table,
                                           CTtable = NWLS_camop,
                                           stationCol = "station",
                                           Xcol = "Longitude",
                                           Ycol= "Latitude",
                                           individualCol = "Individual",
                                           recordDateTimeCol     = "DateTimeOriginal",
                                           recordDateTimeFormat  = "%d-%m-%Y %H:%M",
                                           occasionLength        = 1,
                                           day1 = "station",
                                           includeEffort         = TRUE,
                                           timeZone              = "Asia/Calcutta"
)
write.csv(tiger_capthist1, "Output/nwls_tiger_capthist.csv")

summary(tiger_capthist1)
plot(tiger_capthist1, tracks = TRUE)

# command to understand how much buffer is required 
suggest.buffer(tiger_capthist1)
hist(unlist(moves(tiger_capthist1)))

#creating a habitat Mask
nwls_mask<-make.mask(traps(tiger_capthist1), buffer = 10000, spacing = 1000, type = "trapbuffer", poly = NWLS, poly.habitat = TRUE)
plot(nwls_mask)


##Null Model-
tiger_secr<-secr.fit(tiger_capthist1, mask = nwls_mask)
##detection probablity + heterogenity-
tiger_secr_go_het<-secr.fit(tiger_capthist1, model = g0~h2, mask = nwls_mask)
##movement parameters + heterogenity-
tiger_secr_sig_het<-secr.fit(tiger_capthist1, model = sigma~h2, mask = nwls_mask)
##Combination
tiger_secr_go_sig_het<-secr.fit(tiger_capthist1, model = list(g0~h2, sigma~h2), mask = nwls_mask)

AIC(tiger_secr, tiger_secr_go_het, tiger_secr_sig_het, tiger_secr_go_sig_het)

tiger_secr
tiger_secr_go_het
tiger_secr_sig_het
tiger_secr_go_sig_het

#summary for best fit model
summary(tiger_secr_sig_het)
tiger_secr_sig_het

region.N(tiger_secr_sig_het)


#addtional code
plot(NWLS)
plot(tiger_capthist1, add = T)
plotMCP(tiger_capthist1, add = T)



