## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)

## ----data_unadj----------------------------------------------------------
system("head unadjusted_unemployment.csv",intern=TRUE)
U1 <- read.table(file="unadjusted_unemployment.csv",sep=",",header=TRUE)
head(U1)

## ----reshape-------------------------------------------------------------
u1 <- t(as.matrix(U1[2:13]))
dim(u1) <- NULL
date <- seq(from=1948,length=length(u1),by=1/12)
plot(date,u1,type="l",ylab="Percent unemployment (unadjusted)")

## ----test----------------------------------------------------------------
plot(u1)

## ----data_adj------------------------------------------------------------
U2 <- read.table(file="adjusted_unemployment.csv",sep=",",header=TRUE)
u2 <- t(as.matrix(U2[2:13]))
dim(u2) <- NULL
plot(date,u1,type="l",ylab="percent",col="red")
lines(date,u2,type="l")
title("Unemployment. Raw (black) and seasonally adjusted (red)")

## ----adjustment_spectrum-------------------------------------------------
u1_ts <- ts(u1,start=1948,frequency=12)
u2_ts <- ts(u2,start=1948,frequency=12)
spectrum(ts.union(u1_ts,u2_ts),spans=c(3,5,3),main="Unemployment. Raw (black) and seasonally adjusted (red)")

## ----plot.ts-------------------------------------------------------------
plot(u1_ts)

## ----bls_filter----------------------------------------------------------
s <- spectrum(ts.union(u1_ts,u2_ts),plot=FALSE)

## ----s_names-------------------------------------------------------------
names(s)

## ----s_spec_dim----------------------------------------------------------
dim(s$spec)

## ----s_transfer----------------------------------------------------------
plot(s$freq,s$spec[,2]/s$spec[,1],type="l",log="y",
  ylab="frequency ratio", xlab="frequency",  
  main="frequency response (dashed lines at 0.9 and 1.1)")
abline(h=c(0.9,1.1),lty="dashed",col="red")

## ----loess---------------------------------------------------------------
u1_loess <- loess(u1~date,span=0.5)
plot(date,u1,type="l",col="red")
lines(u1_loess$x,u1_loess$fitted,type="l")

## ----loess_transfer------------------------------------------------------
s2 <- spectrum(ts.union(
  u1_ts,ts(u1_loess$fitted,start=1948,frequency=12)),
  plot=FALSE)
plot(s2$freq,s2$spec[,2]/s$spec[,1],type="l",log="y",
  ylab="frequency ratio", xlab="frequency", xlim=c(0,1.5),
  main="frequency response (dashed line at 1.0)")
abline(h=1,lty="dashed",col="red")

## ----cycles,fig.height=6-------------------------------------------------
u_low <- ts(loess(u1~date,span=0.5)$fitted,start=1948,frequency=12)
u_hi <- ts(u1 - loess(u1~date,span=0.1)$fitted,start=1948,frequency=12)
u_cycles <- u1 - u_hi - u_low
plot(ts.union(u1, u_low,u_hi,u_cycles),
  main="Decomposition of unemployment as trend + noise + cycles")

## ----freq_response-------------------------------------------------------
spec_cycle <- spectrum(ts.union(u1_ts,u_cycles),
  spans=c(3,3),
  plot=FALSE)
freq_response_cycle <- spec_cycle$spec[,2]/spec_cycle$spec[,1]
plot(spec_cycle$freq,freq_response_cycle,
  type="l",log="y",
  ylab="frequency ratio", xlab="frequency", xlim=c(0,1.2), ylim=c(5e-6,1.1),
  main="frequency response (dashed line at 1.0)")
abline(h=1,lty="dashed",col="red")  


## ----show range----------------------------------------------------------
cut_fraction <- 0.5
plot(spec_cycle$freq,freq_response_cycle,
  type="l",log="y",
  ylab="frequency ratio", xlab="frequency", xlim=c(0,0.9), ylim=c(1e-4,1.1),
  main=paste("frequency response, showing region for ratio >", cut_fraction))
abline(h=1,lty="dashed",col="red")  
freq_cycles <- range(spec_cycle$freq[freq_response_cycle>cut_fraction]) 
abline(v=freq_cycles,lty="dashed",col="blue") 
abline(h=cut_fraction,lty="dashed",col="blue")

