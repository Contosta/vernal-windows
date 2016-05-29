#######################################################################################################
#
#
# This script identifies the day of the year where the breakpoint occurs
# in the best-fit piecewise linear regression for two lines in the time series
# from mid-winter through mid-spring.

# Developed as part of the analysis of the research under the New Hampshire ESPCoR Ecosystems and Society grant. 
# Code developed by D Burchsted, A Contasta, A Adolph, and M Green. 
# See companion paper by Contasta et al 2016.
# Date of most recent script draft: January 10, 2016.

# This script was used for identification of the onset of warming in streams, rivers, and soils.
# Below is a simplified version of the full script that was used, demonstrating the algorithms and functions.
# The full scripts also looped through multiple sites and multiple years and created output files and plots.

# Uses a Monte Carlo approach:

#  1. Across 1000 iterations, vary the parameters that define the dataset for analysis: 
#   - smooth the data using rollmedian, varying the size of the window for the median calculation 
#   - vary the start date used to subset the data for analysis: February 15 plus or minus 15 days
#   - vary the end date: May 15 plus or minus 15 days 

#  2. For each iteration, calculate the day of year of the threshold
#	(i.e., they doy for the breakpoint of the best-fit linear regression 
#	 of two lines across the time series)

#  3. The reported threshold = the mode of the calculated thresholds for all 1000 iterations 

#  4. Confidence intervals = the 2.5 and 97.5 exceedance values of the 1000 iterations

# Data used in this sample script:
#   instream temperature data from the New Hampshire LoVoTECS network, available
#   at the New Hampshire EPSCoR Data Discovery Center, downloaded to a local
#   folder in the working directory

#	
#
#######################################################################################################


###
### set-up
###


## call libraries
library( zoo ) # for rolling averages
library( "segmented" ) # for breakpoint linear regression


## user parameters

num.iter <- 1000 # set the # of iterations for the monte carlo
break.est.var <- 15 # variability of initial breakpoint calculation, in days


# start & end dates for smoothing
dates <- data.frame(
	start = c( as.POSIXlt( "2012-12-15" ), as.POSIXlt( "2013-12-15" ) ),
	end = c( as.POSIXlt( "2013-06-30" ), as.POSIXlt( "2014-06-30" ) )
)

# start and end dates for threshold analysis; will be modified +/- 15d for each Monte Carlo iteration
d.start <- "02-01"
d.end <- "05-15"


## define function to determine the mode of a distribution
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


###
### main loop to read & analyze data 
###

	## get data
	
	# read in data for current site
	df <- read.csv(
		sites[i], # point to the file with the site data
		na.strings = "-9999", # define the values or strings that indicate no data
		# this shows the structure of the input data file that was used
		colClasses = c( "character", "numeric", "numeric", "numeric", "numeric" )
	)
	
	# break if no data
	if( length(df[,1]) == 0 ) next
	
	
	## calculate stats
	
		# set up
		year = as.numeric( strftime( dates$start[j], "%Y" ) ) # this defines the year of analysis; the index of j is based on a loop through all years that is not included in this script
		
		# create holder to store breakpoints for each Monte Carlo iteration
		breakpoints <- data.frame( 
			smoothing_window = NA,
			min_date = strptime(NA, format="%m/%d/%Y"),
			max_date = strptime(NA, format="%m/%d/%Y"),
			breakpoint = strptime(NA, format="%m/%d/%Y")
		)
		
		# subset data
		z.data <- zoo( df$TEMP_C, as.POSIXct( df$TS, origin="1970-01-01" ) ) # creates zoo object from the data to merge with z.reg
		z.data <- z.data[!duplicated(index(z.data))] # delete duplicate values
		
		# aggregate data to hourly
		index.hourly <- as.POSIXct( round( as.double(index(z.data)) / (60*60) ) * (60*60), 
									origin="1970-01-01" ) # drop the minutes & seconds from the time stamps of the data
		z.hour <- aggregate( z.data, index.hourly, mean, na.rm=TRUE, na.action=NULL)
		rm(z.data, index.hourly)
		
		# store temp observations in a regular zoo object with hourly intervals
		z.reg <- zoo( , seq( dates$start[j]-120*24*60*60, dates$end[j]+120*24*60*60, by=60*60 ) )
		z <- window( merge( z.reg, z.hour ), index( z.reg ) ) 
		rm(z.hour, z.reg)
				
		# break for data sets that are too small
		if( length(z[!is.na(z)]) == 0 ) next
		if( index( tail(z[!is.na(z)],1) ) < as.POSIXct( paste( year+1, "-01-01", sep="" ), origin="1970-01-01" ) ) next
		d1.test <- as.POSIXct( paste( year+1, "-01-01", sep="" ), origin="1970-01-01" )
		d2.test <- as.POSIXct( paste( year+1, "-07-15", sep="" ), origin="1970-01-01" )
		z.test <- window( z, start=d1.test, end=d2.test )
		if( length( z.test[!is.na(z.test)] )  <  
				0.5 * length( seq(from=d1.test, to=d2.test, by="hour") ) )  next
		rm( d1.test, d2.test, z.test )
		
		# initial estimate of window for piecewise regression: 
		d.min <- as.POSIXct( paste( year+1, "-", d.start, sep=""), origin="1970-01-01" )
		d.max <- as.POSIXct( paste( year+1, "-", d.end, sep=""), origin="1970-01-01" )
		
		
		## calculations of piecewise linear regression & breakpoint
		
		# ~~~ begin loop through Monte Carlo iteration ~~~
		for( k in 1:num.iter ) {
			
			## generate smoothed curve
			
			# randomly choose an integer between 5-120 to be the # days of the window for smoothing
			w <- sample( seq( from=5, to=120 ), size=1  ) 
			# break if smoothing window is > length of data
			if( w > length(z[!is.na(z)])/24 ) next
			# calculate median of the window: window is above w in hours plus one to have an odd # for median
			z.median <- rollmedian( z[!is.na(z)], k=w*24+1 )
			
			
			## select data for regression
			
			# set min & max dates for regression
			var <- break.est.var * 24 * 60 * 60 # variation of breakpoint estimate in seconds
			d1 <- sample( seq( d.min-var, d.min+var, by="day"), size=1 ) # generate random value within variability of window
			d2 <- sample( seq( d.max-var, d.max+var, by="day"), size=1 ) # generate random value within variability of window
			d1 <- strptime( d1, format="%Y-%m-%d") # extract day, drop time
			d2 <- strptime( d2, format="%Y-%m-%d" ) # extract day, drop time
			
			# write parameters for regression to breakpoints table
			breakpoints[k,1] <- w
			breakpoints$min_date[k] <- d1
			breakpoints$max_date[k] <- d2
			
			# break if data window is less than 30 days
			if( d2-d1 < 30 ) 	next
			
			# subset data for breakpoint regression
			z.for.regression <- window( z.median, index = seq(from=d1, to=d2, by="day") )
			if( length(z.for.regression) == 0 ) next
			if( index(tail(z.for.regression,1))-index(head(z.for.regression,1)) < 30 ) 	next
			
			## calc breakpoint
			
			# first linear regression
			x <- as.numeric( index(z.for.regression) ) 
			y <- coredata( z.for.regression )
			
			# piecewise regression with error handling
			break.estimate <- ( as.numeric(d2) - as.numeric(d1) ) / 2
			psi <- head(x,1) + break.estimate
			rm(break.estimate)
			segmented.mod <- try( 
				segmented( lm(y~x),
						   seg.Z=~x,
						   psi=list(x=c(psi))
				), silent = TRUE
			) # try function allows for error handling
			
			# determine date of breakpoint & store in dataframe
			if( class(segmented.mod) == "try-error" ) { next
			} else {
				brk.pt <- as.POSIXct( segmented.mod$psi[2], origin="1970-01-01" )
				breakpoints[k, 4] <- format( brk.pt, "%Y-%m-%d" )
			}
						
			# ~~~ end loop through monte carlo iterations ~~~
		}
				
		
		## calc final threshold & confidence intervals
		
		# identify the two most common breakpoints
		breaks <- breakpoints$breakpoint[ !is.na(breakpoints$breakpoint) ]
		mode1 <- Mode( breaks )
		mode2 <- Mode( breaks[ breaks!=mode1 ] )
		
		# identify the confidence interval
		ci.lo <- quantile(breaks, 0.025)
		ci.hi <- quantile(breaks, 0.975)
		
		cat( paste("Threshold day of year = ", strftime( mode1, "%j" ), " Confidence interval = ", strftime( ci.lo, "%j" ), "-", strftime( ci.hi, "%j" ), "\n") )
				
