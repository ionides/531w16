## ----opts,include=FALSE--------------------------------------------------
library(knitr)
prefix <- "sp500"
opts_chunk$set(
  progress=TRUE,
  prompt=FALSE,tidy=FALSE,highlight=TRUE,
  strip.white=TRUE,
  warning=FALSE,
  message=FALSE,
  error=FALSE,
  echo=TRUE,
  cache=TRUE,
#  cache=FALSE,
  cache.extra=list(rand_seed,R.version.string),
  results='markup',
  fig.show='asis',
  size='small',
  fig.lp="fig:",
  fig.path=paste0("figure/",prefix,"-"),
  cache.path=paste0("cache/",prefix,"-"),
  fig.pos="h!",
  fig.align='center',
  fig.height=5,fig.width=6.83,
  dpi=300,
  dev='png',
  dev.args=list(bg='transparent')
  )

options(cores=20)

library(ggplot2)
theme_set(theme_bw())

## ----sp500---------------------------------------------------------------
dat <- read.table("sp500.csv",sep=",",header=TRUE)
N <- nrow(dat)
sp500 <- dat$Close[N:1] # data are in reverse order in sp500.csv
par(mfrow=c(1,2))
plot(sp500,type="l")
plot(log(sp500),type="l")

## ----garch---------------------------------------------------------------
require(tseries)
fit.garch <- garch(sp500,grad = "numerical", trace = FALSE)
L.garch <- logLik(fit.garch)

## ----data----------------------------------------------------------------
load(file="sp500-2002-2012.rda")
plot(sp500.ret.demeaned,xlab="business day (2002-2012)",ylab="demeaned S&P 500 return", type="l")

## ----load_pomp-----------------------------------------------------------
require(pomp)

## ----names---------------------------------------------------------------
sp500_statenames <- c("H","G","Y_state")
sp500_rp_names <- c("sigma_nu","mu_h","phi","sigma_eta")
sp500_ivp_names <- c("G_0","H_0")
sp500_paramnames <- c(sp500_rp_names,sp500_ivp_names)
sp500_covarnames <- "covaryt"

## ----rproc---------------------------------------------------------------
rproc1 <- "
  double beta,omega,nu;
  omega = rnorm(0,sigma_eta * sqrt( 1- phi*phi ) * sqrt(1-tanh(G)*tanh(G)));
  nu = rnorm(0, sigma_nu);
  G += nu;
  beta = Y_state * sigma_eta * sqrt( 1- phi*phi );
  H = mu_h*(1 - phi) + phi*H + beta * tanh( G ) * exp(-H/2) + omega;
"
rproc2.sim <- "
  Y_state = rnorm( 0,exp(H/2) );
 "

rproc2.filt <- "
  Y_state = covaryt;
 "
sp500_rproc.sim <- paste(rproc1,rproc2.sim)
sp500_rproc.filt <- paste(rproc1,rproc2.filt)

## ----initializer---------------------------------------------------------
sp500_initializer <- "
  G = G_0;
  H = H_0;
  Y_state = rnorm( 0,exp(H/2) );
"

## ----measure-------------------------------------------------------------
sp500_rmeasure <- "
   y=Y_state;
"

sp500_dmeasure <- "
   lik=dnorm(y,0,exp(H/2),give_log);
"

## ----transforms----------------------------------------------------------
sp500_toEstimationScale <- "
  Tsigma_eta = log(sigma_eta);
  Tsigma_nu = log(sigma_nu);
  Tphi = logit(phi);
"

sp500_fromEstimationScale <- "
  Tsigma_eta = exp(sigma_eta);
  Tsigma_nu = exp(sigma_nu);
  Tphi = expit(phi);
"

## ----sp_pomp-------------------------------------------------------------
sp500.filt <- pomp(data=data.frame(y=sp500.ret.demeaned,
                     time=1:length(sp500.ret.demeaned)),
              statenames=sp500_statenames,
              paramnames=sp500_paramnames,
              covarnames=sp500_covarnames,
              times="time",
              t0=0,
              covar=data.frame(covaryt=c(0,sp500.ret.demeaned),
                     time=0:length(sp500.ret.demeaned)),
              tcovar="time",
              rmeasure=Csnippet(sp500_rmeasure),
              dmeasure=Csnippet(sp500_dmeasure),
              rprocess=discrete.time.sim(step.fun=Csnippet(sp500_rproc.filt),delta.t=1),
              initializer=Csnippet(sp500_initializer),
              toEstimationScale=Csnippet(sp500_toEstimationScale), 
              fromEstimationScale=Csnippet(sp500_fromEstimationScale)
)

## ----sim_pomp------------------------------------------------------------
expit<-function(real){1/(1+exp(-real))}
logit<-function(p.arg){log(p.arg/(1-p.arg))}
params_test <- c(
     sigma_nu = exp(-4.5),  
     mu_h = -0.25,  	 
     phi = expit(4),	 
     sigma_eta = exp(-0.07),
     G_0 = 0,
     H_0=0
  )

sim1.sim <- pomp(sp500.filt, 
               statenames=sp500_statenames,
               paramnames=sp500_paramnames,
               covarnames=sp500_covarnames,
               rprocess=discrete.time.sim(step.fun=Csnippet(sp500_rproc.sim),delta.t=1)
)

sim1.sim <- simulate(sim1.sim,seed=1,params=params_test)

## ----build_sim1.filt-----------------------------------------------------

sim1.filt <- pomp(sim1.sim, 
  covar=data.frame(
    covaryt=c(obs(sim1.sim),NA),
    time=c(timezero(sim1.sim),time(sim1.sim))),
  tcovar="time",
  statenames=sp500_statenames,
  paramnames=sp500_paramnames,
  covarnames=sp500_covarnames,
  rprocess=discrete.time.sim(step.fun=Csnippet(sp500_rproc.filt),delta.t=1)
)



## ----run_level-----------------------------------------------------------
run_level <- 3 
sp500_Np <-          c(100,1e3,2e3)
sp500_Nmif <-        c(10, 100,200)
sp500_Nreps_eval <-  c(4,  10,  20)
sp500_Nreps_local <- c(10, 20, 20)
sp500_Nreps_global <-c(10, 20, 100)

## ----parallel-setup,cache=FALSE------------------------------------------
require(doParallel)
registerDoParallel()

## ----pf1-----------------------------------------------------------------
stew(file=sprintf("pf1.rda",run_level),{
  t.pf1 <- system.time(
    pf1 <- foreach(i=1:sp500_Nreps_eval[run_level],.packages='pomp',
                   .options.multicore=list(set.seed=TRUE)) %dopar% try(
                     pfilter(sim1.filt,Np=sp500_Np[run_level])
                   )
  )
},seed=493536993,kind="L'Ecuyer")
(L.pf1 <- logmeanexp(sapply(pf1,logLik),se=TRUE))

## ----mif-----------------------------------------------------------------
sp500_rw.sd_rp <- 0.02
sp500_rw.sd_ivp <- 0.1
sp500_cooling.fraction.50 <- 0.5

stew("mif1.rda",{
   t.if1 <- system.time({
   if1 <- foreach(i=1:sp500_Nreps_local[run_level],
                  .packages='pomp', .combine=c,
                  .options.multicore=list(set.seed=TRUE)) %dopar% try(
                    mif2(sp500.filt,
                         start=params_test,
                         Np=sp500_Np[run_level],
                         Nmif=sp500_Nmif[run_level],
                         cooling.type="geometric",
                         cooling.fraction.50=sp500_cooling.fraction.50,
                         transform=TRUE,
                         rw.sd = rw.sd(
                            sigma_nu  = sp500_rw.sd_rp,
                            mu_h      = sp500_rw.sd_rp,
                            phi       = sp500_rw.sd_rp,
                            sigma_eta = sp500_rw.sd_rp,
                            G_0       = ivp(sp500_rw.sd_ivp),
                            H_0       = ivp(sp500_rw.sd_ivp)
                         )
                    )
                  )
    
    L.if1 <- foreach(i=1:sp500_Nreps_local[run_level],.packages='pomp',
                      .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% 
                      {
                        logmeanexp(
                          replicate(sp500_Nreps_eval[run_level],
                                    logLik(pfilter(sp500.filt,params=coef(if1[[i]]),Np=sp500_Np[run_level]))
                          ),
                          se=TRUE)
                      }
  })
},seed=318817883,kind="L'Ecuyer")

r.if1 <- data.frame(logLik=L.if1[,1],logLik_se=L.if1[,2],t(sapply(if1,coef)))
if (run_level>1) 
  write.table(r.if1,file="sp500_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)
summary(r.if1$logLik,digits=5)

## ----pairs---------------------------------------------------------------
pairs(~logLik+sigma_nu+mu_h+phi+sigma_eta,data=subset(r.if1,logLik>max(logLik)-20))

## ----box-----------------------------------------------------------------
sp500_box <- rbind(
 sigma_nu=c(0.005,0.05),
 mu_h    =c(-1,0),
 phi = c(0.95,0.99),
 sigma_eta = c(0.5,1),
 G_0 = c(-2,2),
 H_0 = c(-1,1)
)

## ----box_eval------------------------------------------------------------
stew(file="box_eval.rda",{
  t.box <- system.time({
    if.box <- foreach(i=1:sp500_Nreps_global[run_level],.packages='pomp',.combine=c,
                  .options.multicore=list(set.seed=TRUE)) %dopar%  
      mif2(
        if1[[1]],
        start=apply(sp500_box,1,function(x)runif(1,x))
      )
    
    L.box <- foreach(i=1:sp500_Nreps_global[run_level],.packages='pomp',.combine=rbind,
                      .options.multicore=list(set.seed=TRUE)) %dopar% {
                        set.seed(87932+i)
                        logmeanexp(
                          replicate(sp500_Nreps_eval[run_level],
                                    logLik(pfilter(sp500.filt,params=coef(if.box[[i]]),Np=sp500_Np[run_level]))
                          ), 
                          se=TRUE)
                      }
  })
},seed=290860873,kind="L'Ecuyer")


r.box <- data.frame(logLik=L.box[,1],logLik_se=L.box[,2],t(sapply(if.box,coef)))
if(run_level>1) write.table(r.box,file="sp500_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)
summary(r.box$logLik,digits=5)

## ----pairs_global--------------------------------------------------------
pairs(~logLik+log(sigma_nu)+mu_h+phi+sigma_eta+H_0,data=subset(r.box,logLik>max(logLik)-10))

## ----garch_benchmark-----------------------------------------------------
require(tseries)
fit.garch.benchmark <- garch(sp500.ret.demeaned,grad = "numerical", trace = FALSE)
L.garch.benchmark <- logLik(fit.garch.benchmark)

