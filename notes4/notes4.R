## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)

## ----root----------------------------------------------------------------
roots <- polyroot(c(1,2,2))
roots

## ----abs_roots-----------------------------------------------------------
abs(roots)

## ----reducibility--------------------------------------------------------
list(AR_roots=polyroot(c(1,-5/6,1/6)),MA_roots=polyroot(c(1,-1,1/4)))

## ----quasi_periodic------------------------------------------------------
omega <- 0.1
ar_coefs <- c(2/(1+omega^2), - 1/(1+omega^2))
set.seed(8395200)
X <- arima.sim(list(ar=ar_coefs),n=500,sd=1)
par(mfrow=c(1,2))
plot(X)
plot(ARMAacf(ar=ar_coefs,lag.max=500),type="l",ylab="ACF of X")

