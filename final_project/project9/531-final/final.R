
dow=read.csv(file = "dow.csv",header = T)

dow_log=diff(log(dow$Adj_Close))

dow_hp=dow_log-mean(dow_log)

require(pomp)

dow_statenames <- c("V","Y_state")
dow_rp_names <- c("sigma_omega","phi","theta")
dow_ivp_names <- c("V_0")
dow_paramnames <- c(dow_rp_names,dow_ivp_names)
dow_covarnames <- "covaryt"

rproc1 <- "
double omega;
omega = rnorm(0,sigma_omega);
V = theta*(1 - phi) + phi*sqrt(V) + sqrt(V)*omega;
"
rproc2.sim <- "
Y_state = rnorm( 0,sqrt(V) );
"
rproc2.filt <- "
Y_state = covaryt;
"
dow_rproc.sim <- paste(rproc1,rproc2.sim)
dow_rproc.filt <- paste(rproc1,rproc2.filt)



dow_initializer <- "
V = V_0;
Y_state = rnorm( 0,sqrt(V) );
"
dow_rmeasure <- "
y=Y_state;
"
dow_dmeasure <- "
lik=dnorm(y,0,sqrt(V),give_log);
"

dow_toEstimationScale <- "
  Tsigma_omega = log(sigma_omega);
  Ttheta = log(theta);
  Tphi = logit(phi);
"
dow_fromEstimationScale <- "
  Tsigma_omega = exp(sigma_omega);
  Ttheta = exp(theta);
  Tphi = expit(phi);
"

dow.filt <- pomp(data=data.frame(y=dow_hp,
                                 time=1:length(dow_hp)),
                 statenames=dow_statenames,
                 paramnames=dow_paramnames,
                 covarnames=dow_covarnames,
                 times="time",
                 t0=0,
                 covar=data.frame(covaryt=c(0,dow_hp),
                                  time=0:length(dow_hp)),
                 tcovar="time",
                 rmeasure=Csnippet(dow_rmeasure),
                 dmeasure=Csnippet(dow_dmeasure),
                 rprocess=discrete.time.sim(step.fun=Csnippet(dow_rproc.filt),delta.t=1),
                 initializer=Csnippet(dow_initializer),
                 toEstimationScale=Csnippet(dow_toEstimationScale), 
                 fromEstimationScale=Csnippet(dow_fromEstimationScale)
)

expit<-function(real){1/(1+exp(-real))}
logit<-function(p.arg){log(p.arg/(1-p.arg))}

params_test <- c(
  sigma_omega = 0.001,  
  phi = 0.001,
  theta = 1,
  V_0= 0.02
)

run_level <- 3 
dow_Np <-          c(100,2e3,5e3)
dow_Nmif <-        c(10, 100,200)
dow_Nreps_eval <-  c(4,  10,  20)
dow_Nreps_local <- c(10, 20, 20)
dow_Nreps_global <-c(10, 20, 40)

require(doParallel)
registerDoParallel(15)

dow_rw.sd_rp <- 0.001
dow_rw.sd_ivp <- 0.001
dow_cooling.fraction.50 <- 0.5

##----------------------------------------------------------------------
stew(file=sprintf("3mif1-%d.rda",run_level),{
  t.if1 <- system.time({
    if1 <- foreach(i=1:dow_Nreps_global[run_level],
                   .packages='pomp', .combine=c,
                   .options.multicore=list(set.seed=TRUE)) %dopar% try(
                     mif2(dow.filt,
                          start=params_test,
                          Np=dow_Np[run_level],
                          Nmif=dow_Nmif[run_level],
                          cooling.type="geometric",
                          cooling.fraction.50=dow_cooling.fraction.50,
                          transform=TRUE,
                          rw.sd = rw.sd(
                            sigma_omega  = dow_rw.sd_rp,
                            theta      = dow_rw.sd_rp,
                            phi       = dow_rw.sd_rp,
                            V_0       = ivp(dow_rw.sd_ivp)
                          )
                     )
                   )
    L.if1 <- foreach(i=1:dow_Nreps_global[run_level],.packages='pomp',
                     .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% 
                     {
                       logmeanexp(
                         replicate(dow_Nreps_eval[run_level],
                                   logLik(pfilter(dow.filt,params=coef(if1[[i]]),Np=dow_Np[run_level]))
                         ),
                         se=TRUE)
                     }
    
  })
},seed=318817883,kind="L'Ecuyer")

##---------------------------------------------------------
dow_box <- rbind(
  sigma_omega=c(0.001,0.05),
  theta    = c(0.1,1),
  phi = c(0.1,1),
  V_0 = c(0.02,0.03)
)

stew(file=sprintf("3box_eval-%d.rda",run_level),{
  t.box <- system.time({
    if.box <- foreach(i=1:dow_Nreps_global[run_level],.packages='pomp',.combine=c,
                      .options.multicore=list(set.seed=TRUE)) %dopar%  
     try( mif2(
        if1[[1]],
        start=apply(dow_box,1,function(x)runif(1,x[1],x[2]))
      ))
    
    L.box <- foreach(i=1:dow_Nreps_global[run_level],.packages='pomp',.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932+i)
                       logmeanexp(
                         replicate(dow_Nreps_eval[run_level],
                                   logLik(pfilter(dow.filt,params=coef(if.box[[i]]),Np=dow_Np[run_level]))
                         ), 
                         se=TRUE)
                     }
  })
},seed=290860873,kind="L'Ecuyer")

plot(if.box)







