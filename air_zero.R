#######################################################################################################
#
#

# One of several scripts used to identify the thresholds that mark the onset of 
# spring within weather, snow, terrestrial, and aquatic systems.

# This script is written to analyze air temperature data and determine the date at which
# the smoothed temperature profile crosses from below 0 degrees C to above 0 degrees C
# for the first time.


# Developed as part of the analysis of the sensor datasets generated as part of 
# research under the New Hampshire ESPCoR Ecosystems and Society grant.  
# Code by A. Adolph with portion including import of netcdf data by Mark Green. 
# See companion paper by Contosta et al (under review, 2016).

# Use a Monte Carlo approach:

#  1. Across 1000 iterations, vary the parameters that define the dataset for analysis: 
#   - smooth the data using rollmedian, varying the size of the window for the median calculation 
#   - vary the start date used to subset the data for analysis: January 15 plus or minus 15 days
#   - vary the end date: May 15 plus or minus 15 days 

#  2. For each iteration, calculate the day of year of the threshold

#  3. The final threshold = the mode of the calculated thresholds for all 1000 iterations 

#  4. Confidence intervals = the 2.5 and 97.5 exceedance values of the 1000 iterations

# Data used in this sample script:
#   1. instream conductivity data from the New Hampshire LoVoTECS network, available
#   at the New Hampshire EPSCoR Data Discovery Center, downloaded to a local
#   folder in the working directory
#	2. temperature thresholds for the same set of site, previously calculated

# Output:
#	a csv file listing thresholds and confidence intervals for each site 

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
library(Imap)
Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# Set up pdf for graphing
pdf("NOHRSC_ZeroTemp_gridded_AllSites_2.pdf")
par(mfrow = c(2,2), mar = c(5, 5, 1, 1))

# User defined variables
g <- 200 # number of iterations to run the piecewise linear regression
break.est.var <- 15

files <- dir(path = "C:/Users/Alden Adolph/Dropbox/Research/Albedo Project/Sensors Team/R_thresholds/gridded_data/ncep_temp")
lats <- substr(files,2,6)
longs <- substr(files,8,13)
lat <- as.numeric(lats)/1000
long <- as.numeric(longs)/1000
num_points <- length(files)
#-----------------------------------
sitedat <- read.table("AllSites_Latitude_Longitude.csv", head = TRUE, sep = ",")
# sites 1-94 are LoVoTECS, sites 95-109 are USGS, sites 110-130 are CoCoRAHS, sites 131-136 are Intensive sites
latlook <- sitedat$Latitude
longlook <- sitedat$Longitude
sitename <- sitedat$Site.Name
num_sites <- length(sitename)
#latc <- sapply(latlook, function(i) which.min(abs(i - lat)))
#longc <- sapply(longlook, function(i) which.min(abs(i - long)))
dist <- matrix(data=NA, nrow=num_points, ncol=1)
filenum <- matrix(data=NA,nrow=num_sites, ncol=1) 


for (j in 1:num_sites){
  for (k in 1:num_points){
    dist[k] <- gdist(long[k], lat[k], longlook[j], latlook[j], units = "nm", a = 6378137.0, b = 6356752.3142, verbose = FALSE)
  }
  filenum[j]<-which.min(dist)
}


# sample distance calc - gdist(long.1, lat.1, lon.2, lat.2, units = "nm", a = 6378137.0, b = 6356752.3142, verbose = FALSE)



# Initialize matricies
output.doys.2012 <- matrix(data=NA, nrow=num_sites, ncol=g)
output.doys.2013 <- matrix(data=NA, nrow=num_sites, ncol=g)
output.doys.2014 <- matrix(data=NA, nrow=num_sites, ncol=g)
final <- matrix(data=NA, nrow=num_sites, ncol=12)
name <-character(length = num_sites)

#------------------------------------------
# Looping through the sites
for( i in 1:num_sites ) {
  rename <- paste0("ncep_temp/",files[filenum[i]], collapse=NULL)
  name[i] <- files[filenum[i]]
  data <- read.table(rename, head = TRUE, sep = ",")
  # Looping thtough the years
  for (year in 2013:2014) {
    start = paste (toString(year-1), "-11-01", sep = "", collapse = NULL)
    end = paste (toString(year), "-06-30", sep = "", collapse = NULL)
    start2 <- as.POSIXct(start)
    end2 <- as.POSIXct(end)
    temp_date<- as.POSIXct(data$date)
    z <- zoo(data$temp_c, temp_date)
    tempyear <- window( z, start = start2, end=end2 )
    
    # Looping through the iterations
    for( p in 1:g )	{	
      w1 <- sample( seq( from=5, to=80 ), size=1  )
      w <- floor( w1 / 2 ) * 2 + 1 # convert numbers to odd integers
      z.median <- rollmedian( tempyear[!is.na(tempyear)], k=w )     
      
      # Finding start day through randomly generated variable
      startseed1 <- paste (toString(year), "-01-15", sep = "", collapse = NULL)
      startseed <- as.POSIXct(startseed1)
      endseed1 <- paste (toString(year), "-05-15", sep = "", collapse = NULL)
      endseed <- as.POSIXct(endseed1)
      var <- break.est.var 
      start.day <- sample( seq( from = startseed - var, to = startseed + var, by = "day"), size = 1 ) # generate random value within variability of window
      end.day <- sample( seq( from = endseed - var, to = endseed + var, by = "day"), size = 1 ) # generate random value within variability of window
      
      if( end.day-start.day < 10 )   {next}
     
      # Running piecewise linear regression
      med.for.regression.temp <- window( z.median, start = start.day, end=end.day )
      temp.series <- coredata(med.for.regression.temp)
      date.series <- index(med.for.regression.temp)
      lengz <- length(med.for.regression.temp)
      # run through the dataset starting at the end and find the 0 threshold
      if (min(temp.series)>0) {next}
      
      for (day in 1:lengz)
      {
        if (temp.series[(lengz-day)]<0)
        {begin.temp <- date.series[(lengz-day)]}
        if (temp.series[(lengz-day)]<0)
        {break}
      }
      
      if (year == 2012)
        {output.doys.2012[i,p]<-begin.temp}
      if (year == 2013)
        {output.doys.2013[i,p]<-begin.temp}
      if (year == 2014)
        {output.doys.2014[i,p]<-begin.temp}
    }
    
    if (year == 2012)
    {outputs<-output.doys.2012[i,]}
    if (year == 2013)
    {outputs<-output.doys.2013[i,]}
    if (year == 2014)
    {outputs<-output.doys.2014[i,]}
    
    if (all(is.na(outputs))) {next} #keep going if no results are achieved for a site in a year
    
    title <- paste (sitename[i], year, sep = " ", collapse = NULL)
    hist(as.POSIXct(outputs, origin = "1970-01-01"),"days",freq = TRUE, xlab="Date",ylab="count",main = title)
    
    startplot = as.POSIXct(paste (toString(year-1), "-12-15", sep = "", collapse = NULL))
    endplot = as.POSIXct(paste (toString(year), "-05-30", sep = "", collapse = NULL))

    
    zplots <- window( tempyear, start = startplot, end=endplot )
    plot(zplots, xlab = "", ylab = "Air Temp",  main = title)
    results <- outputs
    
    abline(v=Mode(results, na.rm=TRUE), col="blue")
    abline(v=quantile(results, 0.025, na.rm = TRUE), col="red")
    abline(v=quantile(results, 0.975, na.rm = TRUE), col="red")
    
    base<- (year-2012)*4
    final[i,base+1]<- round((((Mode(results, na.rm=TRUE))/(3600*24)/365.25)-(year-1970))*365.25,digits=0)
    final[i,base+2]<- round((((quantile(results, 0.025, na.rm = TRUE))/(3600*24)/365.25)-(year-1970))*365.25,digits=0)
    final[i,base+3]<- round((((quantile(results, 0.975, na.rm = TRUE))/(3600*24)/365.25)-(year-1970))*365.25,digits=0)
    final[i,base+4]<- round((((median(results, na.rm=TRUE))/(3600*24)/365.25)-(year-1970))*365.25,digits=0)
    
    #resultsout[i, ]<- results
    name[i] <- sitename[i] 
  }
}

output_final <- data.frame(cbind(lapply(sitename, as.character),final))
my.df <- data.frame(lapply(output_final, as.character), stringsAsFactors=)
#output_results <- data.frame(cbind(name,resultsout))
#write.csv(output_results, file="FullResults_NOHRSC_SWE_Grid.csv")
write.csv(my.df, file = "Thresholds_NOHRSC_ZeroTemp_Allsites_1.csv")
dev.off()
