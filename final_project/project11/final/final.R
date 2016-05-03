set.seed(100000)
library(ggplot2)
library(plyr)
library(reshape2)
library(pomp)
library(tseries)
library(doParallel)
library(foreach)
library(doMC)

data = read.csv("fx.csv")
data$Rate = as.numeric(as.character(data$Rate))
# Remove 2790, which is closing market
dt = na.omit(data)
plot(dt$Rate,type = 'l')
plot(dt$Rate[2601:2900],type = 'l')
# Take the data subset with 400 samples
dt2 = dt[2601:2900,]
dt2$Date=1:300
# pomp
fx = pomp(dt2,times="Date",t0=0)
plot(fx)
# delta t = 1
stochStep <- Csnippet("
                      e = rnorm(0,sigma);
                      N = N*exp((mu-delta*delta/2)/300+delta/sqrt(300)*e);
                      ")
pomp(fx,rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),
     paramnames=c("mu","delta","sigma"),statenames=c("N","e")) -> fx
set.seed(100)
sim <- simulate(fx,params=c(N.0=1.5,e.0=0,mu=0.1,delta=0.2,sigma=0.3),as.data.frame=TRUE,states=TRUE)
plot(dt2$Rate,type = 'l')
lines(N~time,data=sim,col = "red",type = 'l')
set.seed(1000000)
sim <- simulate(fx,params=c(N.0=1.5,e.0=0,mu=0.1,delta=0.3,sigma=0.1),as.data.frame=TRUE,states=TRUE,nsim=3, include.data = TRUE)
sim$N[1:300] = sim$Rate[1:300]
ggplot(data=sim,mapping=aes(x=time,y=N,colour = sim))+geom_line()

# logged, t = 1, no difference
dt3 = log(dt2$Rate)
plot(dt3, type = 'l')

# add measure, from sd normal
dmeas <- Csnippet("lik = dnorm(Rate,0,N,give_log);")
rmeas <- Csnippet("Rate = rnorm(0,N);")
Ne_initializer <- "
 N=rpois(1.5);
 e=rpois(1);
"
stopifnot(packageVersion("pomp")>="0.75-1")
pomp(data = dt2,
     times="Date",
     t0=0,
     rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),
     rmeasure = rmeas,
     dmeasure=dmeas, 
     obsnames = "Rate",
     paramnames=c("mu","delta","sigma"),
     statenames=c("N","e"),
     initializer=Csnippet(Ne_initializer)
     ) -> fx



#test = c(N.0=1.5,e.0=0,mu=1,delta=2,sigma=2)
# likelihood
run_level = 1
level_Np = c(100,1000,2000)
level_Nmif = c(10,100,200)
level_Nreps_eval = c(4,10,20)
level_Nreps_local = c(10,20,20)
level_Nreps_global = c(10,20,100)
registerDoParallel()

# slice
sliceDesign(
  c(mu=0.1,delta=0.2,sigma=0.4),
  mu=rep(seq(from=-10,to=10,length=40),each=3),
  delta=rep(seq(from=0.1,to=3,length=40),each=3),
  sigma=rep(seq(from=0.1,to=3,length=40),each=3)
  ) -> p

registerDoMC(cores=5)
set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)

foreach (theta=iter(p,"row"),.combine=rbind,
         .inorder=FALSE,.options.multicore=mcopts) %dopar% 
         {
           pfilter(fx,params=unlist(theta),Np=5000) -> pf
           pf
           theta$loglik <- logLik(pf)
           theta
         } -> p

foreach (v=c("mu","delta","sigma")) %do% 
{
  x <- subset(p,slice==v)
  plot(x[[v]],x$loglik,xlab=v,ylab="loglik")
}

test = c(N.0=1.5,e.0=0,mu=-0.26,delta=0.62,sigma=1.2)
fx_box <- rbind(
  mu = c(-3,3),
  delta = c(0.5,0.9),
  sigma = c(1,1.5)
)



# filter on simulated data
registerDoParallel()
stew(file=sprintf("pf1.rda",run_level),{
  t.pf1 <- system.time(
    pf1 <- foreach(i=1:level_Nreps_eval[run_level],.packages='pomp',
                   .options.multicore=list(set.seed=TRUE)) %dopar% try(
                     pfilter(fx,params=test,
                             Np=level_Np[run_level])
                   )
  )
},seed=493536993,kind="L'Ecuyer")
logmeanexp(sapply(pf1,logLik),se=TRUE)

# fitting model to data
fx.sd_rp <- 0.002
fx.sd_ivp <- 0.1
fx_cooling.fraction.50 <- 0.1

stew("mif1.rda",{
  t.if1 <- system.time({
    if1 <- foreach(i=1:level_Nreps_local[run_level],
                   .packages='pomp', .combine=c,
                   .options.multicore=list(set.seed=TRUE)) %dopar% try(
                     mif2(fx,
                          start=test,
                          Np=level_Np[run_level],
                          Nmif=level_Nmif[run_level],
                          cooling.type="geometric",
                          cooling.fraction.50=fx_cooling.fraction.50,
                          transform=TRUE,
                          rw.sd = rw.sd(
                            mu = fx.sd_rp,
                            delta = fx.sd_rp,
                            sigma = fx.sd_rp
                          )
                     )
                   )
    
    L.if1 <- foreach(i=1:level_Nreps_local[run_level],.packages='pomp',
                     .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% 
                     {
                       logmeanexp(
                         replicate(level_Nreps_eval[run_level],
                                   logLik(pfilter(fx,params=coef(if1[[i]]),Np=level_Np[run_level]))
                         ),
                         se=TRUE)
                     }
  })
},seed=318817883,kind="L'Ecuyer")

r.if1 <- data.frame(logLik=L.if1[,1],logLik_se=L.if1[,2],t(sapply(if1,coef)))
if (run_level>1) 
  write.table(r.if1,file="fx_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)
summary(r.if1$logLik,digits=5)

pairs(~logLik+mu+delta+sigma,data=subset(r.if1,logLik>max(logLik)-20))


#  Likelihood maximization using randomized starting values

stew(file="box_eval.rda",{
  t.box <- system.time({
    if.box <- foreach(i=1:level_Nreps_global[run_level],.packages='pomp',.combine=c,
                      .options.multicore=list(set.seed=TRUE)) %dopar%  
      mif2(
        if1[[1]],
        start=apply(fx_box,1,function(x)runif(1,x))
      )
    
    L.box <- foreach(i=1:level_Nreps_global[run_level],.packages='pomp',.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932+i)
                       logmeanexp(
                         replicate(level_Nreps_eval[run_level],
                                   logLik(pfilter(fx,params=coef(if.box[[i]]),Np=level_Np[run_level]))
                         ), 
                         se=TRUE)
                     }
  })
},seed=290860873,kind="L'Ecuyer")


r.box <- data.frame(logLik=L.box[,1],logLik_se=L.box[,2],t(sapply(if.box,coef)))
if(run_level>1) write.table(r.box,file="sp500_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)
summary(r.box$logLik,digits=5)

pairs(~logLik+mu+delta+sigma,data=subset(r.box,logLik>max(logLik)-10))

plot(if.box)

# GARCH, compare likelihood, a little better
fxg = garch(dt2$Rate,order = c(0,1),coef = NULL, itmax = 200,eps=NULL, grad = c("analytic"))
logLik(fxg)
plot(fxg)

