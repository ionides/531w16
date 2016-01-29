## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)

## ----eigen---------------------------------------------------------------
N <- 100
phi <- 0.8
sigma <- 1
V <- matrix(NA,N,N)
for(m in 1:N) for(n in 1:N) V[m,n] <- sigma^2 * phi^abs(m-n) / (1-phi^2)
V_eigen <- eigen(V,symmetric=TRUE)
oldpars <- par(mfrow=c(1,2))
matplot(V_eigen$vectors[,1:5],type="l")
matplot(V_eigen$vectors[,6:9],type="l")
par(oldpars)

## ----evals---------------------------------------------------------------
round(V_eigen$values[1:9],2)

