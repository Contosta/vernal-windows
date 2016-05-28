gages = c("01052500", "01064801", "01064500", "01073500", "01092000", "01137500", "01152500", "01161000", "01072800", "01076500", "01082000", "01086000", "01089500", "01094000", "01091000")

doys = matrix(0, length(gages), 4)

for(j in 1:length(gages)){
gage=gages[j]
url=paste("http://waterdata.usgs.gov/nh/nwis/dv?cb_00060=on&format=rdb&site_no=",gage,"&referred_module=sw&period=&begin_date=2012-01-01&end_date=2014-12-31",sep="")
data=read.table(url, sep="\t", header=T)
data1=data[2:length(data[,1]),]
q=as.numeric(as.character(data1[,4]))
dt=strptime(data1[,3], format="%Y-%m-%d")

q. = subset(q, dt$year==112)
dt. = subset(dt, dt$year==112)
dt.. = as.numeric(dt.)

plot(dt., q., log="y", type="l")

#------randomly find the threshold----

#------test window size------
n = 1000
results = matrix(0, n, 3)
for(z in 1:n){

sz = runif(1, 10, 60)
qn = runif(1, 0.1, 0.4)
#---- user inputs------

x = dt.. # put your independent variable here
y = q. # put your dependent variable here

Kh = sz*86400 # this is the window size
step = 86400 # this is the step increment

#-------------------------------------
lowest = min(x)+(Kh/2)
highest = max(x)-(Kh/2)

trend = matrix(0, length(y), 3)
index = 1

for(i in seq(lowest, highest, step)) {
minKh = i - (Kh/2)
maxKh = i + (Kh/2)
a = subset(y, x <= maxKh & x >= minKh)
trend[index,1] = i
trend[index,2] = length(a)
trend[index,3] = quantile(a, qn, na.rm=T)
index = index+1
}

trend2 = subset(trend, trend[,1]!=0)
trend2. = trend[1:(length(trend2[,1])/2.5),]
date2 = as.POSIXct(trend2.[,1], origin="1970-01-01")
date3 = strptime(date2, format="%Y-%m-%d %H:%M:%S")



kk = subset(date3$yday, trend2.[,3]==max(trend2.[,3], na.rm=T))

results[z,1] = sz
results[z,2] = qn
results[z,3] = kk[1]
}

step1 = subset(dt., dt.>as.POSIXct(quantile(results[,3], 0.025, na.rm=T)*86400, origin="2012-01-01") & dt.<as.POSIXct(quantile(results[,3], 0.975, na.rm=T)*86400, origin="2012-01-01"))

step2 = subset(q., dt.>as.POSIXct(quantile(results[,3], 0.025, na.rm=T)*86400, origin="2012-01-01") & dt.<as.POSIXct(quantile(results[,3], 0.975, na.rm=T)*86400, origin="2012-01-01"))

step3 = subset(step1, step2==max(step2))
step4 = subset(step2, step2==max(step2))

doys[j,1] = quantile(results[,3], 0.025, na.rm=T)
doys[j,2] = quantile(results[,3], 0.5, na.rm=T)
doys[j,3] = quantile(results[,3], 0.975, na.rm=T)
doys[j,4] = step3$yday

}

setwd("/Users/mbgreen/Google Drive/SnowSensorData/Discharge")
write.table(data.frame(gages, doys), "doys_2012.csv", sep=",", row.names=F, col.names=c("gage", "Pct2.5", "Pct50", "Pct97.5", "peak"))