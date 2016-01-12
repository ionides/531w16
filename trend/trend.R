## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)

## ----data----------------------------------------------------------------
global_temp <- read.table("Global_Temperature.txt",header=TRUE)
str(global_temp)

## ----data_plot-----------------------------------------------------------
plot(Annual~Year,data=global_temp,ty="l")

## ----glob_temp_lm--------------------------------------------------------
lm_fit <- lm(Annual~Year+I(Year^2),data=global_temp)
summary(lm_fit)

## ----glob_temp_lm_plot---------------------------------------------------
yr <- 1880:2026
Z <- cbind(1,yr,yr^2)
beta <- coef(lm_fit)
prediction <- Z%*%beta
plot(Annual~Year,data=global_temp,ty="l",xlim=range(yr),ylim=range(c(global_temp$Annual,prediction),na.rm=TRUE),lty="dashed")
lines(x=yr,y=prediction,col="red")

## ----acf_global_temp-----------------------------------------------------
acf(resid(lm_fit))

