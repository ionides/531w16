# Reported monthly cases of mumps, 1928-1972 New York City

x<-ts(scan("mumps.dat"),start=1928,freq=12)
pdf(file="mumps.pdf",width=8,height=4)
plot(x,xlab="") #, ylab="Monthly cases of mumps")
dev.off()

pdf(file="mumps-spec.pdf",width=8,height=4)
#spectrum(x) #,main="Monthly mumps in NYC")
spectrum(x,spans=c(3,5,7),main="",sub="",ylim=c(0.8e2,1e6))
text(0,1e5,"(1)")
text(0.35,1.3e5,"(2)")
text(1,5.5e5,"(3)")
text(1.35,3.5e4,"(4)")
text(2,2.5e4,"(5)")
#,main="monthly mumps in NYC")
#dev.copy2eps(file="mumps-spec.eps")
dev.off()

#spectrum(x,method="ar",main="monthly mumps in NYC")

#acf(x)

#arima

aic.table<-function(dat,I=0,Nar=3,Nma=3,...){
 aic<-matrix(NA,Nar+1,Nma+1)
 for(i in 0:Nar){
  for(j in 0:Nma){ 
    if(i+j != 0){
      aic[i+1,j+1]<-(arima(dat,
        order=c(i,I,j),...))$aic}
  }
 }
 aic
}
#aic.table(x,seasonal=c(1,1,1),Nar=4,Nma=4) 

# (aic<-aic.table(log(x),seasonal=c(0,1,1),Nar=4,Nma=4) )
#          [,1]      [,2]      [,3]       [,4]      [,5]
#[1,]        NA  312.7628   92.7453  -42.91403 -131.8598
#[2,] -213.9453 -211.9458 -224.4315 -227.35447 -225.5215
#[3,] -211.9459 -212.5350 -236.8260 -223.70594 -224.7305
#[4,] -224.9618 -237.9834 -236.1537 -234.41349 -232.4061
#[5,] -229.8224 -221.1222 -236.4941 -235.21621 -239.7320

M1=arima(x,order=c(3,0,0),seasonal=c(0,1,1))
#spectrum(resid(M1),method="ar")
#cpgram(resid(M1))

pdf(file="m1-resid.pdf",width=5,height=5)
acf(resid(M1))
text(0.5,0.85,"(a1)",cex=2)
dev.off()

pdf(file="m1-qq.pdf",width=5,height=5)
qqnorm(resid(M1))
text(-1.8,500,"(a2)",cex=2)
dev.off()

pdf(file="m1-resid2.pdf",width=5,height=5)
acf(abs(resid(M1)))
text(0.5,0.85,"(a3)",cex=2)
dev.off()

M2<-arima(log(x),order=c(3,0,0),seasonal=c(0,1,1))
#spectrum(resid(M2),method="ar")
#cpgram(resid(M2))

pdf(file="m2-resid.pdf",width=5,height=5)
acf(resid(M2))
text(0.5,0.85,"(b1)",cex=2)
dev.off()

pdf(file="m2-qq.pdf",width=5,height=5)
qqnorm(resid(M2))
text(-1.8,0.35,"(b2)",cex=2)
dev.off()

pdf(file="m2-resid2.pdf",width=5,height=5)
acf(abs(resid(M2)))
text(0.5,0.85,"(b3)",cex=2)
dev.off()


(m3<-arima(log(x),order=c(3,0,0),seasonal=c(1,1,1)))

matplot(matrix(x[1:528],nr=12),ty="l",lty=1)


