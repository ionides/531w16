## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)
require(pomp)
require(ggplot2)

## ----load_data-----------------------------------------------------------
S=read.table("SP500table.csv",header=TRUE,sep=',')
Date=S$Date
S=rev(S$Close)
plot(S,type="l")

## ----building_BS---------------------------------------------------------
## BS

SP_statenames <- c("H","sp_state")
SP_rp_names_BS <- c("alpha","sigma","eta")
SP_ivp_names <- c("H_0")
SP_paramnames_BS <- c(SP_rp_names_BS,SP_ivp_names)
SP_covarnames <- "Co_sp"

SP_rproc1_BS <- "
  double Z;
  Z = rnorm(0,1);
  H = log(sp_state)+alpha-sigma*sigma/2+sigma*Z;
"
SP_rproc2.sim <- "
  sp_state = exp(H)*exp(rnorm(0,eta));
"

SP_rproc2.filt <- "
  sp_state = Co_sp;
"
SP_rproc_BS.sim <- paste(SP_rproc1_BS,SP_rproc2.sim)
SP_rproc_BS.filt <- paste(SP_rproc1_BS,SP_rproc2.filt)

SP_initializer <-"
  H = H_0;
  sp_state = exp(H)*exp(rnorm(0,eta));
"
SP_rmeasure <-"
  sp=sp_state;
"
SP_dmeasure <-"
  lik=dnorm(log(sp),H,eta,give_log);
"
SP_toEst_BS <-"
  Tsigma=log(sigma);
  Teta=log(eta);
"
SP_fromEst_BS <-"
  Tsigma=exp(sigma);
  Teta=exp(eta);
"

SP_BS.filt <- pomp(data=data.frame(sp=S[-1],time=1:(length(S)-1)),
                statenames=SP_statenames,
                paramnames=SP_paramnames_BS,
                covarnames=SP_covarnames,
                times="time",
                t0=0,
                covar=data.frame(Co_sp=S,time=0:(length(S)-1)),
                tcovar="time",
                rmeasure=Csnippet(SP_rmeasure),
                dmeasure=Csnippet(SP_dmeasure),
                rprocess=discrete.time.sim(step.fun=Csnippet(SP_rproc_BS.filt),
                                           delta.t=1),
                initializer=Csnippet(SP_initializer),
                toEstimationScale=Csnippet(SP_toEst_BS),
                fromEstimationScale=Csnippet(SP_fromEst_BS)
)


## ----building_MJD--------------------------------------------------------
SP_rp_names_MJD <- c("alpha","sigma","lambda","mu","delta","eta")
SP_paramnames_MJD <- c(SP_rp_names_MJD,SP_ivp_names)

SP_rproc1_MJD <- "
  double N1,Z,sumY,k;
  N1 = rpois(lambda);
  Z = rnorm(0,1);
  sumY=0;
  if (N1>0){
    sumY = rnorm(N1*mu,delta*sqrt(N1));
  };
  k = exp(mu+1/2*delta*delta)-1;
  H = log(sp_state)+alpha-sigma*sigma/2-lambda*k+sigma*Z+sumY;
"

SP_rproc_MJD.sim <- paste(SP_rproc1_MJD,SP_rproc2.sim)
SP_rproc_MJD.filt <- paste(SP_rproc1_MJD,SP_rproc2.filt)

SP_toEst_MJD <-"
  Tsigma=log(sigma);
  Tdelta=log(delta);
  Tlambda=log(lambda);
  Teta=log(eta);
"
SP_fromEst_MJD <-"
  Tsigma=exp(sigma);
  Tdelta=exp(delta);
  Tlambda=exp(lambda);
  Teta=exp(eta);
"

SP_MJD.filt<- pomp(data=data.frame(sp=S[-1],time=1:(length(S)-1)),
                statenames=SP_statenames,
                paramnames=SP_paramnames_MJD,
                covarnames=SP_covarnames,
                times="time",
                t0=0,
                covar=data.frame(Co_sp=S,time=0:(length(S)-1)),
                tcovar="time",
                rmeasure=Csnippet(SP_rmeasure),
                dmeasure=Csnippet(SP_dmeasure),
                rprocess=discrete.time.sim(step.fun=Csnippet(SP_rproc_MJD.filt),
                                           delta.t=1),
                initializer=Csnippet(SP_initializer),
                toEstimationScale=Csnippet(SP_toEst_MJD),
                fromEstimationScale=Csnippet(SP_fromEst_MJD)
)

## ----build_HMJD----------------------------------------------------------
SP_statenames_HMJD <- c("H","sp_state","I")
SP_rp_names_HMJD <- c("alpha","sigma","lambda","mu","delta","eta","p")
SP_ivp_names_HMJD <- c("H_0","I_0")
SP_paramnames_HMJD <- c(SP_rp_names_HMJD,SP_ivp_names_HMJD)

SP_rproc1_HMJD <- "
  double N1,Z,sumY,k;
  I = rbinom(1,p);
  N1 = rpois(lambda);
  Z = rnorm(0,1);
  sumY=0;
  if (N1>0){
    sumY = rnorm(N1*mu,delta*sqrt(N1));
  };
  k = exp(mu+1/2*delta*delta)-1;
  H = log(sp_state)+alpha-sigma*sigma/2+sigma*Z-I*(lambda*k-sumY);
"
SP_rproc_HMJD.sim <- paste(SP_rproc1_HMJD,SP_rproc2.sim)
SP_rproc_HMJD.filt <- paste(SP_rproc1_HMJD,SP_rproc2.filt)

SP_initializer_HMJD <-"
  H = H_0;
  I = 0;
  if (I_0>0.5)  I = 1;
  sp_state = exp(H)*exp(rnorm(0,eta));
"
SP_toEst_HMJD <-"
  Tsigma=log(sigma);
  Tdelta=log(delta);
  Tlambda=log(lambda);
  Teta=log(eta);
  TI_0=logit(I_0);
  Tp=logit(p);
"
SP_fromEst_HMJD <-"
  Tsigma=exp(sigma);
  Tdelta=exp(delta);
  Tlambda=exp(lambda);
  Teta=exp(eta);
  TI_0=expit(I_0);
  Tp=expit(p);
"

SP_HMJD.filt <- pomp(data=data.frame(sp=S[-1],time=1:(length(S)-1)),
                statenames=SP_statenames_HMJD,
                paramnames=SP_paramnames_HMJD,
                covarnames=SP_covarnames,
                times="time",
                t0=0,
                covar=data.frame(Co_sp=S,time=0:(length(S)-1)),
                tcovar="time",
                rmeasure=Csnippet(SP_rmeasure),
                dmeasure=Csnippet(SP_dmeasure),
                rprocess=discrete.time.sim(step.fun=Csnippet(SP_rproc_HMJD.filt),
                                           delta.t=1),
                initializer=Csnippet(SP_initializer_HMJD),
                toEstimationScale=Csnippet(SP_toEst_HMJD),
                fromEstimationScale=Csnippet(SP_fromEst_HMJD)
)

## ----simulation_HJMD-----------------------------------------------------
params_test <-c(
  sigma=0.01,
  eta=0.005,
  H_0=log(1280),
  alpha=-0.0001,
  mu=0.5,
  lambda=0.002,
  delta=0.005,
  I_0=0.5,
  p=0.5
  )

sim1_HMJD.sim <-pomp(SP_HMJD.filt,
                statenames=SP_statenames_HMJD,
                paramnames=SP_paramnames_HMJD,
                covarnames=SP_covarnames,
                rprocess=discrete.time.sim(step.fun=Csnippet(SP_rproc_HMJD.sim),delta.t=1)
                )
sim <- simulate(sim1_HMJD.sim,seed=1,params=params_test,nsim=20,as=TRUE,include=TRUE)
ggplot(sim,mapping=aes(x=time,y=sp,group=sim,color=sim=='data'))+
  geom_line()+guides(color=FALSE)


## ----runlevel------------------------------------------------------------
run_level <-3
SP_Np = c(100,1e3,5e3)
SP_Nmif=c(10,100,200)
SP_Neval=c(10,20,30)
SP_Nglobal=c(10,10,100)
SP_Nlocal=c(10,10,20)

## ----Parallel------------------------------------------------------------
require(doParallel)
cores <-20
registerDoParallel(cores)
mcopts = list(set.seed=TRUE)
set.seed(931129,kind = "L'Ecuyer")


## ----pfilter-------------------------------------------------------------

##pfilter for BS
stew(file=sprintf("pf_BS-%d.rda",run_level),{
  t_pf_BS <- system.time(
    pf_BS <- foreach(i=1:SP_Neval[run_level],.packages='pomp',
                  .options.multicore=mcopts) %dopar% try(
                    pfilter(SP_BS.filt,params=params_test[1:4],Np=SP_Np[run_level],filter.mean=TRUE)
                  )
  )
  
},seed=1320290398,kind="L'Ecuyer")

L_pf_BS <- logmeanexp(sapply(pf_BS,logLik),se=TRUE)

##pfilter for MJD
stew(file=sprintf("pf_MJD-%d.rda",run_level),{
  t_pf_MJD <- system.time(
    pf_MJD <- foreach(i=1:SP_Neval[run_level],.packages='pomp',
                  .options.multicore=mcopts) %dopar% try(
                    pfilter(SP_MJD.filt,params=params_test[1:7],Np=SP_Np[run_level],filter.mean=TRUE)
                  )
  )
  
},seed=1320290398,kind="L'Ecuyer")

L_pf_MJD <- logmeanexp(sapply(pf_MJD,logLik),se=TRUE)

## ----result_pfilter------------------------------------------------------
L_pf_BS
L_pf_MJD


## ----local_MLE-----------------------------------------------------------
SP_rw.sd_rp <- 0.02
SP_rw.sd_ivp <- 0.1
SP_cooling.fraction.50 <- 0.5

#BS
stew(file=sprintf("mif1_BS-%d.rda",run_level),{
  t.if1_BS <- system.time({
    if1_BS <- foreach(i=1:SP_Nlocal[run_level],
                   .packages='pomp', .combine=c,
                   .options.multicore=list(set.seed=TRUE)) %dopar% try(
                     mif2(SP_BS.filt,
                          start=params_test[1:4],
                          Np=SP_Np[run_level],
                          Nmif=SP_Nmif[run_level],
                          cooling.type="geometric",
                          cooling.fraction.50=SP_cooling.fraction.50,
                          transform=TRUE,
                          rw.sd = rw.sd(
                            alpha  = SP_rw.sd_rp,
                            sigma  = SP_rw.sd_rp,
                            eta    = SP_rw.sd_rp,
                            H_0    = ivp(SP_rw.sd_ivp)
                          )
                     )
                   )
    
    L.if1_BS <- foreach(i=1:SP_Nlocal[run_level],.packages='pomp',
                     .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% 
{
  logmeanexp(
    replicate(SP_Neval[run_level],
              logLik(pfilter(SP_BS.filt,params=coef(if1_BS[[i]]),Np=SP_Np[run_level]))
    ),
    se=TRUE)
}
  })
},seed=318817883,kind="L'Ecuyer")

r.if1_BS <- data.frame(logLik=L.if1_BS[,1],logLik_se=L.if1_BS[,2],t(sapply(if1_BS,coef)))
if (run_level>1) 
  write.table(r.if1_BS,file="sp_params_BS.csv",append=TRUE,col.names=FALSE,row.names=FALSE,sep=',')


#MJD
stew(file=sprintf("mif1_MJD-%d.rda",run_level),{
  t.if1_MJD <- system.time({
    if1_MJD <- foreach(i=1:SP_Nlocal[run_level],
                   .packages='pomp', .combine=c,
                   .options.multicore=list(set.seed=TRUE)) %dopar% try(
                     mif2(SP_MJD.filt,
                          start=params_test[1:7],
                          Np=SP_Np[run_level],
                          Nmif=SP_Nmif[run_level],
                          cooling.type="geometric",
                          cooling.fraction.50=SP_cooling.fraction.50,
                          transform=TRUE,
                          rw.sd = rw.sd(
                            alpha  = SP_rw.sd_rp,
                            mu     = SP_rw.sd_rp,
                            sigma  = SP_rw.sd_rp,
                            eta    = SP_rw.sd_rp,
                            lambda = SP_rw.sd_rp,
                            delta  = SP_rw.sd_rp,
                            H_0    = ivp(SP_rw.sd_ivp)
                          )
                     )
                   )
    
    L.if1_MJD <- foreach(i=1:SP_Nlocal[run_level],.packages='pomp',
                     .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% 
{
  logmeanexp(
    replicate(SP_Neval[run_level],
              logLik(pfilter(SP_MJD.filt,params=coef(if1_MJD[[i]]),Np=SP_Np[run_level]))
    ),
    se=TRUE)
}
  })
},seed=318817883,kind="L'Ecuyer")

r.if1_MJD <- data.frame(logLik=L.if1_MJD[,1],logLik_se=L.if1_MJD[,2],t(sapply(if1_MJD,coef)))
if (run_level>1) 
  write.table(r.if1_MJD,file="sp_params_MJD.csv",append=TRUE,col.names=FALSE,row.names=FALSE,sep=',')

#HMJD
stew(file=sprintf("mif1_HMJD-%d.rda",run_level),{
  t.if1_HMJD <- system.time({
    if1_HMJD <- foreach(i=1:SP_Nlocal[run_level],
                   .packages='pomp', .combine=c,
                   .options.multicore=list(set.seed=TRUE)) %dopar% try(
                     mif2(SP_HMJD.filt,
                          start=params_test,
                          Np=SP_Np[run_level],
                          Nmif=SP_Nmif[run_level],
                          cooling.type="geometric",
                          cooling.fraction.50=SP_cooling.fraction.50,
                          transform=TRUE,
                          rw.sd = rw.sd(
                            alpha  = SP_rw.sd_rp,
                            mu     = SP_rw.sd_rp,
                            sigma  = SP_rw.sd_rp,
                            eta    = SP_rw.sd_rp,
                            lambda = SP_rw.sd_rp,
                            delta  = SP_rw.sd_rp,
                            p      = SP_rw.sd_rp,
                            H_0    = ivp(SP_rw.sd_ivp),
                            I_0    = ivp(SP_rw.sd_ivp)
                          )
                     )
                   )
    
    L.if1_HMJD <- foreach(i=1:SP_Nlocal[run_level],.packages='pomp',
                     .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% 
{
  logmeanexp(
    replicate(SP_Neval[run_level],
              logLik(pfilter(SP_HMJD.filt,params=coef(if1_HMJD[[i]]),Np=SP_Np[run_level]))
    ),
    se=TRUE)
}
  })
},seed=318817883,kind="L'Ecuyer")

r.if1_HMJD <- data.frame(logLik=L.if1_HMJD[,1],logLik_se=L.if1_HMJD[,2],t(sapply(if1_HMJD,coef)))
if (run_level>1) 
  write.table(r.if1_HMJD,file="sp_params_HMJD.csv",append=TRUE,col.names=FALSE,row.names=FALSE,sep=',')

## ----results_local-------------------------------------------------------
summary(r.if1_BS$logLik,digits=5)
summary(r.if1_MJD$logLik,digits=5)
summary(r.if1_HMJD$logLik,digits=5)

## ----local_geometry------------------------------------------------------
pairs(~logLik+log(sigma)+mu+eta+alpha+log(delta)+log(lambda)+H_0+p+I_0,data=subset(r.if1_HMJD,logLik>max(logLik)-400))


## ----parameter_box-------------------------------------------------------
SP_box <- rbind(
  sigma=c(0.005,0.1),
  alpha=c(-0.02,0.02),
  mu = c(-4,2),
  lambda=c(0.001,0.02),
  eta = c(0.01,0.1),
  delta=c(0.001,0.1),
  p=c(0.01,0.6),
  H_0 = c(6.5,8),
  I_0=c(0.2,0.8)
)

## ----global_MLE----------------------------------------------------------
###BS
stew(file=sprintf("box_eval_BS-%d.rda",run_level),{
  t.box_BS <- system.time({
    if.box_BS <- foreach(i=1:SP_Nglobal[run_level],.packages='pomp',.combine=c,
                      .options.multicore=list(set.seed=TRUE)) %dopar%  
      mif2(
            SP_BS.filt,
            start=apply(SP_box[c(1,2,5,8),],1,function(x){
              ans=runif(1,x[1],x[2])
              if (min(x)>0) ans=exp(runif(1,log(x[1]),log(x[2])))
              return(ans)}),
            Np=SP_Np[run_level],
            Nmif=SP_Nmif[run_level],
            cooling.type="geometric",
            cooling.fraction.50=SP_cooling.fraction.50,
            transform=TRUE,
            rw.sd = rw.sd(
                    alpha  = SP_rw.sd_rp,
                    sigma  = SP_rw.sd_rp,
                    eta    = SP_rw.sd_rp,
                    H_0    = ivp(SP_rw.sd_ivp)
                    )
      )
    
    L.box_BS <- foreach(i=1:SP_Nglobal[run_level],.packages='pomp',.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932+i)
                       logmeanexp(
                         replicate(SP_Neval[run_level],
                                   logLik(pfilter(SP_BS.filt,params=coef(if.box_BS[[i]]),Np=SP_Np[run_level]))
                         ), 
                         se=TRUE)
                     }
  })
},seed=290860873,kind="L'Ecuyer")


r.box_BS <- data.frame(logLik=L.box_BS[,1],logLik_se=L.box_BS[,2],t(sapply(if.box_BS,coef)))
if(run_level>1) write.table(r.box_BS,file="sp500_params_box_BS.csv",append=TRUE,col.names=FALSE,row.names=FALSE)




###MJD
stew(file=sprintf("box_eval_MJD-%d.rda",run_level),{
  t.box_MJD <- system.time({
    if.box_MJD <- foreach(i=1:SP_Nglobal[run_level],.packages='pomp',.combine=c,
                      .options.multicore=list(set.seed=TRUE)) %dopar%  
      mif2(
            SP_MJD.filt,
            start=apply(SP_box[c(1:6,8),],1,function(x){
              ans=runif(1,x[1],x[2])
              if (min(x)>0) ans=exp(runif(1,log(x[1]),log(x[2])))
              return(ans)}),
            Np=SP_Np[run_level],
            Nmif=SP_Nmif[run_level],
            cooling.type="geometric",
            cooling.fraction.50=SP_cooling.fraction.50,
            transform=TRUE,
            rw.sd = rw.sd(
                            alpha  = SP_rw.sd_rp,
                            mu     = SP_rw.sd_rp,
                            sigma  = SP_rw.sd_rp,
                            eta    = SP_rw.sd_rp,
                            lambda = SP_rw.sd_rp,
                            delta  = SP_rw.sd_rp,
                            H_0    = ivp(SP_rw.sd_ivp)
                    )
      )
    
    L.box_MJD <- foreach(i=1:SP_Nglobal[run_level],.packages='pomp',.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932+i)
                       logmeanexp(
                         replicate(SP_Neval[run_level],
                                   logLik(pfilter(SP_MJD.filt,params=coef(if.box_MJD[[i]]),Np=SP_Np[run_level]))
                         ), 
                         se=TRUE)
                     }
  })
},seed=290860873,kind="L'Ecuyer")


r.box_MJD <- data.frame(logLik=L.box_MJD[,1],logLik_se=L.box_MJD[,2],t(sapply(if.box_MJD,coef)))
if(run_level>1) write.table(r.box_MJD,file="sp500_params_box_MJD.csv",append=TRUE,col.names=FALSE,row.names=FALSE)


###HMJD
stew(file=sprintf("box_eval_HMJD-%d.rda",run_level),{
  t.box_HMJD <- system.time({
    if.box_HMJD <- foreach(i=1:SP_Nglobal[run_level],.packages='pomp',.combine=c,
                      .options.multicore=list(set.seed=TRUE)) %dopar%  
      mif2(
            SP_HMJD.filt,
            start=apply(SP_box,1,function(x){
              ans=runif(1,x[1],x[2])
              if (min(x)>0) ans=exp(runif(1,log(x[1]),log(x[2])))
              return(ans)}),
            Np=SP_Np[run_level],
            Nmif=SP_Nmif[run_level],
            cooling.type="geometric",
            cooling.fraction.50=SP_cooling.fraction.50,
            transform=TRUE,
            rw.sd = rw.sd(
                            alpha  = SP_rw.sd_rp,
                            mu     = SP_rw.sd_rp,
                            sigma  = SP_rw.sd_rp,
                            eta    = SP_rw.sd_rp,
                            lambda = SP_rw.sd_rp,
                            delta  = SP_rw.sd_rp,
                            p      = SP_rw.sd_rp,
                            H_0    = ivp(SP_rw.sd_ivp),
                            I_0    = ivp(SP_rw.sd_ivp)
                    )
      )
    
    L.box_HMJD <- foreach(i=1:SP_Nglobal[run_level],.packages='pomp',.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932+i)
                       logmeanexp(
                         replicate(SP_Neval[run_level],
                                   logLik(pfilter(SP_HMJD.filt,params=coef(if.box_HMJD[[i]]),Np=SP_Np[run_level]))
                         ), 
                         se=TRUE)
                     }
  })
},seed=290860873,kind="L'Ecuyer")


r.box_HMJD <- data.frame(logLik=L.box_HMJD[,1],logLik_se=L.box_HMJD[,2],t(sapply(if.box_HMJD,coef)))
if(run_level>1) write.table(r.box_HMJD,file="sp500_params_box_HMJD.csv",append=TRUE,col.names=FALSE,row.names=FALSE)

## ----results_global------------------------------------------------------
summary(r.box_BS$logLik,digits=5)
summary(r.box_MJD$logLik,digits=5)
summary(r.box_HMJD$logLik,digits=5)

## ----global_geometry-----------------------------------------------------
pairs(~logLik+log(sigma)+mu+eta+alpha+log(delta)+log(lambda)+H_0+p+I_0,data=subset(r.box_HMJD,logLik>max(logLik)-400))

## ----diagnostic----------------------------------------------------------
plot(if.box_MJD)
plot(if.box_HMJD)

## ----likelihood_surface--------------------------------------------------
para=r.box_HMJD[which.max(r.box_HMJD$logLik),3:11]
s=50
expand.grid(alpha=seq(from=-0.005,to=0.01,length=s),
            mu=seq(from=-0.5,to=0.5,length=s)) -> p
p_grid=cbind(p,(sapply(para[c(1,4:9)],rep,s*s)))


stew(file=sprintf("lik_grid23_HJMD_1_%d.rda",run_level),
     {p <- foreach(i=1:(s*s),.packages='pomp',.combine=rbind,
                      .options.multicore=list(set.seed=TRUE)) %dopar% try(
 {
   pfilter(SP_HMJD.filt, params=unlist(p_grid[i,]),Np=5000) -> pf
   theta=p_grid[i,]
   theta$loglik <- logLik(pf)
   theta
 })
 },seed=290860873,kind="L'Ecuyer")

require(dplyr)
pp <- mutate(p,loglik=ifelse(loglik>max(loglik)-600,loglik,NA))
ggplot(data=pp,mapping=aes(x=alpha,y=mu,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  geom_contour(color='black',binwidth=100)+
  geom_point(aes(x=0,y=0),color="red")+
  scale_fill_gradient()+
  labs(x=expression(alpha),y=expression(mu))


## ----profile_alpha-------------------------------------------------------

## construct_parameter_box
It=30
nprof=40
profile.box <- profileDesign(  
  alpha=seq(-0.01,0.01,length.out=It),
  lower=SP_box[-2,1],
  upper=SP_box[-2,2],
  nprof=nprof
)

## mif1
stew(file=sprintf("profile alpha-1-%d.rda",It),{
  
  t_global.prf1 <- system.time({
      prof.llh<- foreach(i=1:(It*nprof),.packages='pomp', .combine=rbind, .options.multicore=mcopts) %dopar%{
        # Find MLE
        mif2(
            SP_HMJD.filt,
            start=unlist(profile.box[i,]),
            Np=3000,
            Nmif=50,
            cooling.type="geometric",
            cooling.fraction.50=SP_cooling.fraction.50,
            transform=TRUE,
            rw.sd = rw.sd(
                            mu     = SP_rw.sd_rp,
                            sigma  = SP_rw.sd_rp,
                            eta    = SP_rw.sd_rp,
                            lambda = SP_rw.sd_rp,
                            delta  = SP_rw.sd_rp,
                            p      = SP_rw.sd_rp,
                            H_0    = ivp(SP_rw.sd_ivp),
                            I_0    = ivp(SP_rw.sd_ivp)
                    )
        )->mif_prf
        # evaluate llh
        evals = replicate(10, logLik(pfilter(mif_prf,Np=3000)))
        ll=logmeanexp(evals, se=TRUE)        
        
        data.frame(as.list(coef(mif_prf)),
                   loglik = ll[1],
                   loglik.se = ll[2])
      }
  })
},seed=931129,kind="L'Ecuyer")

t_global.prf1

## filiter again on the maxima

prof.llh %>%   
  mutate(alpha=signif(alpha,digits=6))%>%
  ddply(~alpha,subset,rank(-loglik)<=10) %>%
  subset(select=SP_paramnames_HMJD) -> pars


## mif2 again
stew(file=sprintf("profile alpha-2-%d.rda",It),{
  
  t_global.prf2 <- system.time({
    prof.llh<- foreach(i=1:(nrow(pars)),.packages='pomp', .combine=rbind, .options.multicore=mcopts) %dopar%{
      # Find MLE
      mif2(
            SP_HMJD.filt,
            start=unlist(pars[i,]),
            Np=5000,
            Nmif=100,
            cooling.type="geometric",
            cooling.fraction.50=SP_cooling.fraction.50,
            transform=TRUE,
            rw.sd = rw.sd(
                            mu     = SP_rw.sd_rp,
                            sigma  = SP_rw.sd_rp,
                            eta    = SP_rw.sd_rp,
                            lambda = SP_rw.sd_rp,
                            delta  = SP_rw.sd_rp,
                            p      = SP_rw.sd_rp,
                            H_0    = ivp(SP_rw.sd_ivp),
                            I_0    = ivp(SP_rw.sd_ivp)
                    )
      )->mif_prf2
      # evaluate llh 
      pf= replicate(10,pfilter(mif_prf2,Np=5000))
      evals=sapply(pf,logLik)
      ll=logmeanexp(evals, se=TRUE)  
      nfail=sapply(pf,getElement,"nfail")
      
      data.frame(as.list(coef(mif_prf2)),
                 loglik = ll[1],
                 loglik.se = ll[2],
                 nfail.max=max(nfail))
    }
  })
},seed=931129,kind="L'Ecuyer")
t_global.prf2

## plot_profile
prof.llh %<>%
  mutate(alpha=signif(alpha,digits=6))%>%
  ddply(~alpha,subset,rank(-loglik)<=2)

a=prof.llh$loglik[which(rank(-prof.llh$loglik)==5)]
b=a-1.92
CI=which(prof.llh$loglik>=b)


prof.llh %>%
  ggplot(aes(x=alpha,y=loglik))+
  geom_point()+
  geom_smooth(method="loess")+
  geom_hline(aes(yintercept=a),linetype="dashed",color="red")


