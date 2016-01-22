## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)

## ----stationarity_sim, echo=FALSE----------------------------------------
N <- 500
times <- 1:N
T1 <- 120
T2 <- 37
set.seed(73413)
y <- sin(2*pi*(times/T1 + runif(1))) +   sin(2*pi*(times/T2 + runif(1))) + rnorm(N)
x <- y[1:50]
oldpars <- par(mfrow=c(1,2))
plot(x,ty="l",xlab="")
plot(y,ty="l",xlab="")
par(oldpars)

## ----sinusoidal,echo=FALSE-----------------------------------------------
np <- 500
U <- seq(from=0,to=1,length=np)
epsilon1 <- sin(2*pi*U)
epsilon2 <- sin(2*pi*2*U)
epsilon3 <- sin(2*pi*3*U)
matplot(U,cbind(epsilon1,epsilon2,epsilon3),col=c("black","red","blue"),lty=c(1,2,4),ylab="",ty="l",xlab="U")
abline(h=0,lty="dotted")
abline(v=c(1/4,1/2,3/4),lty="dotted")


## ----ar_arima_sim,fig.width=4--------------------------------------------
set.seed(123456789)
ar1 <- arima.sim(list(ar=0.6),n=100,sd=1)
plot(ar1,type="l")

## ----ar_sim,fig.width=4--------------------------------------------------
set.seed(123456789)
N <- 100
X <- numeric(N)
X[1] <- rnorm(1,sd=1.56)
for(n in 2:N) X[n] <- 0.6 * X[n-1] + rnorm(1)
plot(X,type="l")

## ----ma_sim--------------------------------------------------------------
N <- 100
set.seed(123456789)
X1 <- arima.sim(list(ma=c(1.5,1)),n=N,sd=1)
set.seed(123456789)
epsilon <- rnorm(N+2)
X2 <- numeric(N)
for(n in 1:N) X2[n] <- epsilon[n+2]+1.5*epsilon[n+1]+epsilon[n]
oldpars <- par(mfrow=c(1,2))
plot(X1,type="l")
plot(X2,type="l")
par(oldpars)

## ----check---------------------------------------------------------------
all(X1==X2)

## ----noise_sim,fig.width=4-----------------------------------------------
N <- 100
set.seed(123456789)
epsilon <- rnorm(N)
plot(epsilon,type="l")

## ----sp500---------------------------------------------------------------
dat <- read.table("sp500.csv",sep=",",header=TRUE)
N <- nrow(dat)
sp500 <- dat$Close[N:1] # data are in reverse order in sp500.csv
par(mfrow=c(1,2))
plot(sp500,type="l")
plot(log(sp500),type="l")

## ----sp500params---------------------------------------------------------
mu <- mean(diff(log(sp500)))
sigma <- sd(diff(log(sp500)))
set.seed(95483123)
X1 <- log(sp500[1])+cumsum(c(0,rnorm(N-1,mean=mu,sd=sigma)))
set.seed(324324587)
X2 <- log(sp500[1])+cumsum(c(0,rnorm(N-1,mean=mu,sd=sigma)))
par(mfrow=c(1,2))
plot(X1,type="l")
plot(X2,type="l")

## ----sp500_acf-----------------------------------------------------------
z <- diff(log(sp500))
acf(z)

## ----sp500_abs_return_acf------------------------------------------------
acf(abs(z-mean(z)),lag.max=200)

