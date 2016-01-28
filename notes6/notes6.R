## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)

## ----data_file-----------------------------------------------------------
system("head huron_depth.csv",intern=TRUE)

## ----read_data-----------------------------------------------------------
dat <- read.table(file="huron_depth.csv",sep=",",header=TRUE)
dat$Date <- strptime(dat$Date,"%m/%d/%Y")
dat$year <- as.numeric(format(dat$Date, format="%Y"))
dat$month <- as.numeric(format(dat$Date, format="%m"))
head(dat)
huron_depth <- dat$Average
time <- dat$year + dat$month/12 # Note: we treat December 2011 as time 2012.0, etc
plot(huron_depth~time,type="l")

## ----sarima--------------------------------------------------------------
huron_sarma11x10 <- arima(huron_depth,
   order=c(1,0,1),
   seasonal=list(order=c(1,0,0),period=12)
)
huron_sarma11x10

## ----residuals-----------------------------------------------------------
acf(resid(huron_sarma11x10))

