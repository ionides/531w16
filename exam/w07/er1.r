um=read.table("um-ertimeseries.txt",header=T)
colnames(um)=c("date","hr","patients","beds","frac",
paste("X",1:17,sep=""))

attach(um)

pdf(file="ER-um-aug-sep05.pdf")
plot(er,ylab="ER occupancy fraction",xlab="Time, in days since the start of August 2005")
dev.off()

pdf(file="ER-um-05.pdf",wi=6.5,he=4)
plot(ts(frac,freq=24),ylab="ER occupancy fraction",xlab="Time, in days since the start of July, 2005")
dev.off()

x=ts(frac,freq=24)
acf(x,lag.max=15*24)

acf(x-lag(x,24),lag.max=15*24)

spectrum(x-lag(x,24),spans=c(5,7,9))
tmp=spectrum(x-lag(x,24),spans=c(7,9,13))

spectrum(ts(frac,freq=24),spans=c(5,7,9))
spectrum(ts(frac,freq=24),spans=c(9,11,13),xlim=c(0,5))
spectrum(x,spans=c(5,7,9),xlim=c(0,5))
spectrum(ts(frac,freq=24),spans=c(7,9,13))#,xlim=c(0,3))
#abline(v=1/7,lty=3)
lines(x=c(1,1)/7,y=c(1e-06,5e-02),lty=3)
spectrum(ts(frac,freq=24),spans=c(7,9,13),xlim=c(0,3))


pdf(file="ER-um-spec-closeup.pdf")
tmp=spectrum(ts(frac,freq=24),spans=c(7,9,13),plot=F)
plot.freq=(tmp$freq<3)
tmp$spec=tmp$spec[plot.freq]
tmp$freq=tmp$freq[plot.freq]
plot.spec(tmp,main="")
text(2,2e-01,expression(bold(A)),cex=2)
dev.off()



pdf(file="ER-um-acf.pdf")
acf(ts(frac,freq=24),main="",lag.max=15*24)
text(10,0.8,expression(bold(B)),cex=2)
dev.off()


pdf(file="ER-um-spec.pdf")
spectrum(ts(frac,freq=24),spans=c(5,7,9),main="")
dev.off()

spectrum(ts(frac,freq=24),method="ar")

spectrum(ts(frac,freq=24))
,method="ar")

spectrum(er,spans=c(3,3,3))

m1=arima(x1,order=c(1,0,1),seasonal=c(1,0,1))

tsdiag(m1)

m2=arima(x1,order=c(1,0,1),seasonal=c(2,0,1))

m3=arima(x1,order=c(1,0,1),seasonal=c(1,0,1))

m4=arima(x1,order=c(1,0,1),seasonal=c(0,1,1))
tsdiag(m4)
acf(resid(m4))
plot(resid(m4))
cpgram(resid(m4))



m5=arima(er,order=c(2,0,1),seasonal=c(0,1,1))




spectrum(ts(frac,freq=24)

acf(ts(frac,freq=24),lag.max=10*24)

#mu=lowess(x,span=0.5)

time=1:length(frac)
mu=ts(loess(frac~time,data=um,span=0.5)$fitted,freq=24)

lo=loess(frac~time,data=um,span=0.5)
lo.p=predict(lo,time=seq(from=1,to=length(frac),length=200),se=T)

tt=seq(from=1,to=length(frac),length=200)
lo.p=predict(lo,
  data.frame(time=tt,se=T)


pdf(file="ER-um-smo.pdf",wi=6.5,he=4)
matplot(y=cbind(lo.p$fit-2*lo.p$se,lo.p$fit,lo.p$fit+2*lo.p$se),x=tt/24,ty="l",lty=c(3,1,3),col=1,ylab="mu.hat",xlab="time (days)")
dev.off()


plot(mu)

pdf(file="ER-um-trend-spec.pdf")
tmp=spectrum(mu,spans=c(7,9,13),plot=F)
plot.freq=(tmp$freq<1.5)
tmp$spec=tmp$spec[plot.freq]
tmp$freq=tmp$freq[plot.freq]
plot.spec(tmp,main="")
#text(2,2e-01,expression(bold(A)),cex=2)
dev.off()

spectrum(mu)

tsmu$fitted)

plot(ts(mu$y,freq=24))

aic.table<-function(dat,I=0,Nar=3,Nma=3,...){
 aic<-matrix(NA,Nar+1,Nma+1)
 for(i in 0:Nar){
  for(j in 0:Nma){ 
    if(i+j != 0){
      aic[i+1,j+1]<-(arima(dat,
        order=c(i,I,j),...))$aic}
  }
 }
 rownames(aic)=paste("AR",c(0:Nar))
 colnames(aic)=paste("MA",c(0:Nma))
 round(aic,1)
}

int1=31*24+1:(61*24) #august + sept 2005
x1= ts( (frac- mu)[int1],freq=24)
u1= ts( (frac)[int1],freq=24)
#note 7/1/2005 was a friday
#     8/1/2005 was a monday

int2=5833:(5833+61*24-1) # march + april 2006
x2= ts( (frac- mu)[int2],freq=24)
u2= ts( (frac)[int2],freq=24)
plot(x2)

plot(cbind(u1,u2))

pdf(file="ER-2case.pdf")
plot(cbind(AugSep05=x1,MarApr06=x2),main="")
dev.off()


spectrum(rbind(AugSep05=x1,MarApr06=x2))

spectrum(x1)

pdf(file="ER-2case-spec.pdf")
spectrum(cbind(x1,x2),span=c(3,5,7),col=1,main="")

dev.off()

par(mfrow=c(2,2))
spectrum(x1)
acf(x1)
spectrum(x2)
acf(x2)
par(mfrow=c(1,1))

aic.table(x1)

a1=aic.table(x1,Nar=4,Nma=4)
(as1=aic.table(x1,Nar=2,Nma=2,seasonal=c(0,1,1)))
(as2=aic.table(x2,Nar=2,Nma=2,seasonal=c(0,1,1)))

(as1.3=aic.table(x1,Nar=3,Nma=3,seasonal=c(0,1,1)))
(as2.3=aic.table(x2,Nar=3,Nma=3,seasonal=c(0,1,1)))
> (as1.3=aic.table(x1,Nar=3,Nma=3,seasonal=c(0,1,1)))
        MA 0    MA 1    MA 2    MA 3
AR 0      NA -1612.9 -2258.7 -2528.5
AR 1 -3060.0 -3058.8 -3057.2 -3055.2
AR 2 -3058.8 -3057.1 -3055.2 -3053.2
AR 3 -3057.2 -3054.8 -3059.1 -3054.5
> (as2.3=aic.table(x2,Nar=3,Nma=3,seasonal=c(0,1,1)))
        MA 0    MA 1    MA 2    MA 3
AR 0      NA -1168.5 -1844.2 -2273.7
AR 1 -2944.9 -2944.4 -2943.2 -2941.7
AR 2 -2944.5 -2943.1 -2941.3 -2940.4
AR 3 -2943.2 -2941.3 -2939.8 -2938.1

f1=arima(x1,order=c(1,0,1),seasonal=c(1,0,1))

f2=arima(x2,order=c(1,0,1),seasonal=c(1,0,1))

f1d=arima(x1,order=c(1,0,1),seasonal=c(0,1,1))
f2d=arima(x2,order=c(1,0,1),seasonal=c(0,1,1))
f1u=arima(u1,order=c(1,0,1),seasonal=c(0,1,1))
f2u=arima(u2,order=c(1,0,1),seasonal=c(0,1,1))

f11=arima(x1,order=c(1,0,0),seasonal=c(0,1,1))
f12=arima(x2,order=c(1,0,0),seasonal=c(0,1,1))

pdf(file="ER-resid-acf.pdf",height=4,width=4)
acf(resid(f1),main="",lag.max=2.5*24)
text(1,0.8,expression(bold(A)),cex=1.5)
dev.off()


pdf(file="ER-resid.pdf",he=4,wi=4)
plot(resid(f1))
dev.off()

pdf(file="ER-qqnorm.pdf",he=4,wi=4)
qqnorm(resid(f1),main="")
text(-2,0.2,expression(bold(B)),cex=1.5)
dev.off()

