
library(dplyr)
data2015<-read.csv("Beijing_2015_HourlyPM25_created20160201.csv",header=T,skip=3)
data2015<-data2015[data2015$Value>0,]
gp2015<-group_by(data2015,Year,Month,Day)
day_data2015<-summarize(gp2015,pm2.5=mean(Value))

data2014<-read.csv("Beijing_2014_HourlyPM25_created20150203.csv",header=T,skip=3)
data2014<-data2014[data2014$Value>0,]
gp2014<-group_by(data2014,Year,Month,Day)
day_data2014<-summarize(gp2014,pm2.5=mean(Value))

data2013<-read.csv("Beijing_2013_HourlyPM2.5_created20140325.csv",header=T,skip=3)
data2013<-data2013[data2013$Value>0,]
gp2013<-group_by(data2013,Year,Month,Day)
day_data2013<-summarize(gp2013,pm2.5=mean(Value))

data2012<-read.csv("Beijing_2012_HourlyPM2.5_created20140325.csv",header=T,skip=3)
data2012<-data2012[data2012$Value>0,]
gp2012<-group_by(data2012,Year,Month,Day)
day_data2012<-summarize(gp2012,pm2.5=mean(Value))

data2011<-read.csv("Beijing_2011_HourlyPM25_created20140709.csv",header=T,skip=3)
data2011<-data2011[data2011$Value>0,]
gp2011<-group_by(data2011,Year,Month,Day)
day_data2011<-summarize(gp2011,pm2.5=mean(Value))

data2010<-read.csv("Beijing_2010_HourlyPM25_created20140709.csv",header=T,skip=3)
data2010<-data2010[data2010$Value>0,]
gp2010<-group_by(data2010,Year,Month,Day)
day_data2010<-summarize(gp2010,pm2.5=mean(Value))

day_data2010_2015<-rbind(day_data2010,day_data2011,day_data2012,day_data2013,day_data2014,day_data2015)

write.table(day_data2010_2015,file="day_data_2010_2015.csv",row.names=F,col.names=T,sep=",")

day_data2010_2015$Time<-paste(day_data2010_2015$Year,day_data2010_2015$Month,day_data2010_2015$Day,sep='/')
day_data2010_2015$Time<-as.Date(day_data2010_2015$Time)

write.table(day_data2010_2015,file="day_data_2010_2015.csv",row.names=F,col.names=T,sep=",")

data<-read.csv("day_data_2010_2015.csv",header=T)
data$Time<-as.Date(data$Time)

temp<-read.csv("temp.csv",header=T)
temp=temp[,c(6,7,12)]
temp$TMAX[temp$TMAX==-9999]<-NA
temp$TMIN[temp$TMIN==-9999]<-NA
temp<-temp[1:1643,]
temp$Time<-as.Date(as.character(temp$DATE),format = "%Y%m%d")

write.table(temp,file="temp.csv",row.names=F,col.names=T,sep=",")

temp_pm<-left_join(data,temp)

temp_pm2<-temp_pm[,c(5,4,7,8)]
temp_pm2$TD<-temp_pm2$TMAX-temp_pm2$TMIN

write.table(temp_pm2,file="temp_pm25.csv",row.names=F,col.names=T,sep=",")

temp_an<-temp_pm2[1:1060,]
write.table(temp_an,file="temp_pm25_2010_2012.csv",row.names=F,col.names=T,sep=",")
