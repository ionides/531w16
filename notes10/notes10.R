## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)

## ----read_data_e---------------------------------------------------------
e_data <- read.table(file="life_expectancy_usa.txt",header=TRUE)
head(e_data)

## ----read_data_u---------------------------------------------------------
u_data <- read.table(file="unadjusted_unemployment.csv",sep=",",header=TRUE)
head(u_data)

## ----clean_data----------------------------------------------------------
t <- intersect(e_data$Yr,u_data$Year)
e <- e_data$e0[e_data$Yr %in% t]
u <- apply(u_data[u_data$Year %in% t, 2:13],1,mean)

## ----data_plots,fig.height=5---------------------------------------------
plot(ts(cbind(e,u),start=1948),main="Percent unemployment (u) and life expectancy (e) for USA",xlab="Year")

## ----hp------------------------------------------------------------------
require(mFilter)
e_hp <- hpfilter(e, freq=100,type="lambda",drift=F)$cycle
u_hp <- hpfilter(u, freq=100,type="lambda",drift=F)$cycle

## ----hpplots, fig.cap="Figure 2. Detrended unemployment (black; left axis) and detrended life expectancy at birth (red; right axis)."----
plot(t,u_hp,type="l",xlab="Year",ylab="")
par(new=TRUE)
plot(t,e_hp,col="red",type="l",axes=FALSE,xlab="",ylab="")
axis(side=4, col="red")

## ----hp_b----------------------------------------------------------------
arima(e_hp,xreg=u_hp,order=c(1,0,0))

## ----lrt-----------------------------------------------------------------
log_lik_ratio <- as.numeric(
   logLik(arima(e_hp,xreg=u_hp,order=c(1,0,0))) -
   logLik(arima(e_hp,order=c(1,0,0)))
)
LRT_pval <- 1-pchisq(2*log_lik_ratio,df=1)

## ----hp_early_fit--------------------------------------------------------
tt <- 1:45
arima(e_hp[tt],xreg=u_hp[tt],order=c(1,0,0))

## ----aic_table-----------------------------------------------------------
aic_table <- function(data,P,Q,xreg=NULL){
  table <- matrix(NA,(P+1),(Q+1))
  for(p in 0:P) {
    for(q in 0:Q) {
       table[p+1,q+1] <- arima(data,order=c(p,0,q),xreg=xreg)$aic
    }
  }
  dimnames(table) <- list(paste("<b> AR",0:P, "</b>", sep=""),paste("MA",0:Q,sep=""))
  table
}
e_aic_table <- aic_table(e_hp,4,5,xreg=u_hp)
require(knitr)
kable(e_aic_table,digits=2)

## ----arma31--------------------------------------------------------------
arima(e_hp,xreg=u_hp,order=c(3,0,1))

## ----resid---------------------------------------------------------------
r <- resid(arima(e_hp,xreg=u_hp,order=c(1,0,0)))
plot(r)

## ----acf-----------------------------------------------------------------
acf(r)

## ----clean_data_again----------------------------------------------------
delta_e <- e - e_data$e0[e_data$Yr %in% (t-1)]

## ----plots, fig.cap="unemployment (black; left axis) and differenced life expectancy (red; right axis)."----
plot(t,u,type="l",xlab="Year",ylab="")
par(new=TRUE)
plot(t,delta_e,col="red",type="l",axes=FALSE,xlab="",ylab="")
axis(side=4,col="red")

## ----arma----------------------------------------------------------------
arima(delta_e,xreg=u,order=c(1,0,1))

