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

## ----weather_data_file---------------------------------------------------
system("head ann_arbor_weather.csv",intern=TRUE)

## ----weather_data--------------------------------------------------------
x <- read.table(file="ann_arbor_weather.csv",header=TRUE)
head(x)
low <- x$Low

## ----replace_na----------------------------------------------------------
low[is.na(low)] <- mean(low, na.rm=TRUE)

## ----periodogram---------------------------------------------------------
spectrum(low, main="Unsmoothed periodogram")

## ----smoothed_periodogram------------------------------------------------
spectrum(low,spans=c(3,5,3), main="Smoothed periodogram",ylim=c(15,100))

## ----ar_periodogram------------------------------------------------------
spectrum(low,method="ar", main="Spectrum estimated via AR model picked by AIC")

## ----poly_fit------------------------------------------------------------
lm0 <- lm(Low~1,data=x)
lm1 <- lm(Low~Year,data=x)
lm2 <- lm(Low~Year+I(Year^2),data=x)
lm3 <- lm(Low~Year+I(Year^2)+I(Year^3),data=x)
poly_aic <- matrix( c(AIC(lm0),AIC(lm1),AIC(lm2),AIC(lm3)), nrow=1,
   dimnames=list("<b>AIC</b>", paste("order",0:3)))
require(knitr)
kable(poly_aic,digits=1)

## ----plot_jan_temp,fig.width=5-------------------------------------------
plot(Low~Year,data=x,type="l")

## ----read_glob_temp------------------------------------------------------
y <- read.table("Global_Temperature.txt",header=TRUE)
global_temp <- y$Annual[y$Year %in% x$Year]
lm_global <- lm(Low~global_temp,data=x)
AIC(lm_global)

## ----glob_temp_fit-------------------------------------------------------
summary(lm_global)

