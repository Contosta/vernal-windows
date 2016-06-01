#######################################################################################################
#
#

# One of several scripts used to identify the thresholds that mark the onset of 
# spring within weather, snow, terrestrial, and aquatic systems.

# This script is written to calculate the threshold date of the end of snow melt,
# defined here as the date where SWE drops to zero and does not exceed zero again.
# Also applied to snow depth in place of swe.

# Developed as part of the analysis of the sensor datasets generated as part of 
# research under the New Hampshire ESPCoR Ecosystems and Society grant.  
# Code by A. Adolph with portion including import of netcdf data by M. Green. 
# See companion paper by Contosta et al (under review, 2016).

# Input data required for script:
#   1. NetCDF file from NOHRSC containing daily snow water equivalent
#   2. Latitude and Longitude coordinates of site points
#   OR (with data input modifications)
#   1. Daily time series of snow depth and snow density from the CoCoRAHS dataset
#      collected through the EPSCoR Program and available for download at the
#      data discovery center online through the Univeristy of New Hampshire

# Output:
#   1. A .csv file listing thresholds for each site 
#	
#
#######################################################################################################


library(RNetCDF)
library(rgdal)
library(raster)
library(foreign)
library(zoo)
library(segmented)
library(matrixStats)
library(hydroTSM)

Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# Opening up netcdf file and extracting relevant data
nc = open.nc("G02158_86562.sub.nc")
file.inq.nc(nc)
print.nc(nc)
#lat
lat1 = var.inq.nc(nc, 10)
lat = var.get.nc(nc, lat1$name, unpack=TRUE)
#long
long1 = var.inq.nc(nc, 11)
long = var.get.nc(nc, long1$name, unpack=TRUE)
#date
date1 = var.inq.nc(nc,1)
date = var.get.nc(nc, date1$name, unpack=TRUE)
dates = as.POSIXct(date*86400,origin ="1601-01-01")
#SWE
swe1 = var.inq.nc(nc, 2)
swe = var.get.nc(nc, swe1$name, unpack=TRUE)
swe[swe == -9999] <- NA



#-----------------------------------
sitedat <- read.table("AllSites_Latitude_Longitude.csv", head = TRUE, sep = ",")
# sites 1-94 are LoVoTECS, sites 95-109 are USGS, sites 110-130 are CoCoRAHS, sites 131-136 are Intensive sites
latlook <- sitedat$Latitude
longlook <- sitedat$Longitude
sitename <- sitedat$Site.Name
latc <- sapply(latlook, function(i) which.min(abs(i - lat)))
longc <- sapply(longlook, function(i) which.min(abs(i - long)))

num_sites <- length(sitename)

# Initialize matricies
final <- matrix(data=NA, nrow=num_sites, ncol=3)
name <-character(length = num_sites)
end.day2 <-matrix(data=NA,nrow=1, ncol=1)
#------------------------------------------
# Looping through the sites
for( i in 1:num_sites ) {
  swesite = swe[longc[i],latc[i],]
  swesitename <- sitename[i]
  # Looping thtough the years
  for (year in 2012:2014) {
    start = paste (toString(year-1), "-11-01", sep = "", collapse = NULL)
    end = paste (toString(year), "-06-30", sep = "", collapse = NULL)
    start2 <- as.POSIXct(start)
    end2 <- as.POSIXct(end)
    z <- zoo(swesite,dates[1:1331])
    sweyear <- window( z, start = start2, end=end2 )
    
      # Finding the first day when the SWE dropped to zero and does not return above zero to set as end
      temp.series <- coredata(sweyear)
      date.series <- index(sweyear)
      lengz <- length(sweyear)
      for (day in 1:lengz)
      {
        if (temp.series[(lengz-day)]>0)
        {end.day <- date.series[(lengz-day)+1]}
        if (temp.series[(lengz-day)]>0)
        {break}
      }
      
    end.day2[1] <- as.POSIXct(end.day, origin="1970-01-01" )
    base<- (year-2012)+1
    final[i,base]<- round((((end.day2)/(3600*24)/365.25)-(year-1970))*365.25,digits=0)
    
    
    #resultsout[i, ]<- results
    name[i] <- swesitename 
  }
}

output_final <- data.frame(cbind(lapply(sitename, as.character),final))
my.df <- data.frame(lapply(output_final, as.character), stringsAsFactors=)
write.csv(my.df, file = "Thresholds_NOHRSC_SnowFree.csv")
dev.off()
