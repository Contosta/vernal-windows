#######################################################################################################
#
#

# One of several scripts used to identify the thresholds that mark the onset of 
# spring within weather, snow, terrestrial, and aquatic systems.

# This script is written to analyze conductivity data and identifies two thresholds:
# (1) peak and (2) trough closest to the temperature threshold

# Developed as part of the analysis of the sensor datasets generated as part of 
# research under the New Hampshire ESPCoR Ecosystems and Society grant. Code 
# developed by M Green and A Contosta. See companion paper by Contosta et al (under review, 2016).

# Use a Monte Carlo approach:

#  1. Across 1000 iterations, vary the parameters that define the dataset for analysis: 
#   - smooth the data using rollmedian, varying the size of the window for the median calculation 
#   - vary the start date used to subset the data for analysis: February 15 plus or minus 15 days
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
#	a csv file listing thresholds and confidence intervals for each site / data file
#	thresholds printed to standard output, as they are calculated

#	
#
#######################################################################################################

#------------------------------------------------------------------------------------------------
#initialize
#------------------------------------------------------------------------------------------------

# set directories
directory.working <- "C://work/rsch - epscor sensors" # working directory
directory.data <- "data_lovotecs/" # local directory with instream data
file.thresholds <- "output_thresholds/thresholds_20160110-1837.csv" # the file with the previously calculated temp thresholds

# user-specified parameters
n = 10 # how many iterations?

# start & end dates for smoothing
dates <- data.frame(
	start = c( as.POSIXlt( "2012-12-15" ), as.POSIXlt( "2013-12-15" ) ),
	end = c( as.POSIXlt( "2013-06-30" ), as.POSIXlt( "2014-06-30" ) )
)

# start and end dates for threshold analysis; will be modified +/- 15d for each Monte Carlo iteration
d.start <- "02-01"
d.end <- "05-15"

d <- c( "2012-2013", "2013-2014" ) # text of years in POR for labeling

begin.m = 1
begin.d = 1
begin.y = 2013

end.m = 12
end.d = 31
end.y = 2013


#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

setwd(directory.working)
files = list.files(directory.data)
sites = substr(files, 10, 12)


begin = as.integer(ISOdatetime(begin.y, begin.m, begin.d, 0, 0, 0, tz="UTC"))
end = as.integer(ISOdatetime(end.y, end.m, end.d, 0, 0, 0, tz="UTC"))


#------------------------------------------------------------------------------------------------
# pull temp thresholds & match with input files
#------------------------------------------------------------------------------------------------

temp.thresholds = read.csv(file.thresholds, header=T)

sites.matched = match(temp.thresholds[,1], sites)
sites.analysis = sites[!is.na(sites.matched)]
site.files = files[!is.na(sites.matched)]

# initialize matrix to store results
results = matrix(0, length(site.files), 2)


#------------------------------------------------------------------------------------------------
# loop through sites and years
#------------------------------------------------------------------------------------------------

for( i in 1:length(site.files)) {

	# read instream conductivity data
	setwd(directory.data)
	d <- read.csv(site.files[j], header=T)
		
	for( j in 1:length(dates[,1])) {
		
		# get previously calculated temp threshold date for this year at this site
		if(j==1) temp.threshold.col = 2 # column with threshold dates for year1
		if(j==2) temp.threshold.col = 5 # column with threshold dates for year2
		temp.thresholds.row = match(sites.analysis[i], temp.thresholds[,1])
		temp.day = temp.thresholds[ temp.thresholds.row, temp.threshold.col ]
		
		# convert temp threshold from day of year to date
		year = strftime( dates$end[j], "%Y" )
		temp.date = as.POSIXlt( temp.day*24*60*60, origin=paste(year,"-01-01", sep=""), tz="GMT" )
		
		
		#------------------------------------------------------------------------------------------------
		# loop through Monte Carlo iterations
		#------------------------------------------------------------------------------------------------
		
		sub.results = matrix(0, n, 2) # a place to hold Monte Carlo iterative results
		
		for(k in 1:n){
			
			# select the size of the window, in days, to look for the conductivity thresholds
			cond.window = runif(1, 15, 40)			
			
			# subset the data to the window
			d$TS = as.POSIXct(d$TS, origin="1970-01-01")
			d1 = d[ d$TS > temp.date-cond.window*86400  &  d$TS < temp.date+cond.window*86400, ]
						
			#------------------------------------------------------------------------------------------------
			#smoothing
			#------------------------------------------------------------------------------------------------
			
			x = data1.[,2] # put your independent variable here
			y = data1.[,4] # put your dependent variable here
			
			y = subset(y, y!=-9999)
			x = subset(x, y!=-9999)
			
			Kh = runif(1, 2, 5)*86400 # this is the window size
			step = 86400 # this is the step increment
			
			if ( length(y)>0){
				if ( max(x, na.rm=T)-min(x, na.rm=T) > Kh ){
					#------------------------------------------------------------------------------------------------
					lowest = min(x, na.rm=T)+(Kh/2)
					highest = max(x, na.rm=T)-(Kh/2)
					
					trend = matrix(0, length(seq(lowest, highest, step)), 3)
					index = 1
					
					for(i in seq(lowest, highest, step)) {
						minKh = i - (Kh/2)
						maxKh = i + (Kh/2)
						a = subset(y, x <= maxKh & x >= minKh)
						trend[index,1] = i
						trend[index,2] = length(a)
						trend[index,3] = median(a, na.rm=T)
						index = index+1
					}
					
					trend2 = subset(trend, trend[,1]!=0)
					
					#------------------------------------------------------------------------------------------------
					
					a1 = subset(trend2, trend2[,3]>quantile(trend2[,3], runif(1, 0.9, 0.98), na.rm=T)) # take the largest 3 percent of data
					a1. = abs(a1[,1]-as.numeric(wtdate2)) # calculate the time gap between the temp threshold and those data
					a1.1 = subset(a1, a1.==min(a1., na.rm=T)) # find the point that is closest date to the temp threshold
					a2 = subset(trend2, trend2[,3]<quantile(trend2[,3], runif(1, 0.02, 0.1), na.rm=T)) # take the smallest 3 percent of data
					a2. = abs(a2[,1]-as.numeric(wtdate2)) # calculate the time gap between the temp threshold and those data
					a2.1 = subset(a2, a2.==min(a2., na.rm=T)) # find the point that is closest date to the temp threshold
					if(length(a1.1[1]==1)){sub.results[k,1] = a1.1[1]}
					if(length(a2.1[1]==1)){sub.results[k,2] = a2.1[1]}
				}
			}
		}
		
		doys1 = strptime(as.POSIXct(sub.results[,1], origin="1970-01-01"), format="%Y-%m-%d %H:%M:%S")$yday
		ux. = unique(doys1)
		ux = ux.[!is.na(ux.)]
		doys2 = ux[which.max(tabulate(match(doys1, ux)))]
		doys3 = median(subset(sub.results[,1], doys1==doys2), na.rm=T)
		results[j,1] = doys2[1]
		
		doys1 = strptime(as.POSIXct(sub.results[,2], origin="1970-01-01"), format="%Y-%m-%d %H:%M:%S")$yday
		ux. = unique(doys1)
		ux = ux.[!is.na(ux.)]
		doys2 = ux[which.max(tabulate(match(doys1, ux)))]
		doys4 = median(subset(sub.results[,2], doys1==doys2), na.rm=T)
		results[j,2] = doys2[1]
		
		#----diagnostic----
		data1. = subset(data1, data1[,2]>as.numeric(wtdate2)-30*86400 & data1[,2]<as.numeric(wtdate2)+30*86400 & data1[,4]>0)
		date1 = as.POSIXct(data1.[,2], origin="1970-01-01")
		plot(date1, data1.[,4], type="l")
		abline(v=as.numeric(wtdate2), col="red")
		abline(v=median(doys3, na.rm=T), col="blue")
		abline(v=median(doys4, na.rm=T))
	}		
		
	}
}




aa = data.frame(sites1, results)
setwd("/Users/mbgreen/Desktop")
write.table(aa, "SC_LoVoTECS_thresh_2014.csv", sep=",", row.names=F, col.names=c("Site", "Peak", "Trough"))
