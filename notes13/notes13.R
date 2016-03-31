## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8"
  )

## ----prelims,echo=F,cache=F----------------------------------------------
set.seed(594709947L)
require(ggplot2)
theme_set(theme_bw())
require(plyr)
require(reshape2)
require(foreach)
require(doMC)
require(pomp)
stopifnot(packageVersion("pomp")>="0.69-1")

## ----sir-construct-------------------------------------------------------
bsflu <- read.table("bsflu_data.txt")

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_init <- Csnippet("
  S = nearbyint(N)-1;
  I = 1;
  R = 0;
  H = 0;
")

dmeas <- Csnippet("lik = dbinom(B,H,rho,give_log);")
rmeas <- Csnippet("B = rbinom(H,rho);")

pomp(bsflu,times="day",t0=0,
     rprocess=euler.sim(sir_step,delta.t=1/5),
     initializer=sir_init,rmeasure=rmeas,dmeasure=dmeas,
     zeronames="H",statenames=c("H","S","I","R"),
     paramnames=c("Beta","gamma","rho","N")) -> sir

## ----bbs-mc-like-2,results='markup'--------------------------------------
simulate(sir,params=c(Beta=2,gamma=1,rho=0.5,N=2600),
         nsim=10000,states=TRUE) -> x
matplot(time(sir),t(x["H",1:50,]),type='l',lty=1,
        xlab="time",ylab="H",bty='l',col='blue')
lines(time(sir),obs(sir,"B"),lwd=2,col='black')

## ----bbs-mc-like-3,results='markup',cache=T------------------------------
ell <- dmeasure(sir,y=obs(sir),x=x,times=time(sir),log=TRUE,
                params=c(Beta=2,gamma=1,rho=0.5,N=2600))
dim(ell)

## ----bbs-mc-like-4,results='markup'--------------------------------------
ell <- apply(ell,1,sum); summary(exp(ell)); logmeanexp(ell,se=TRUE)

## ----sir-sim1------------------------------------------------------------
sims <- simulate(sir,params=c(Beta=2,gamma=1,rho=0.8,N=2600),nsim=20,
                 as.data.frame=TRUE,include.data=TRUE)

ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)

## ----sir-pfilter-1,results='markup',cache=T------------------------------
pf <- pfilter(sir,Np=5000,params=c(Beta=2,gamma=1,rho=0.8,N=2600))
logLik(pf)

## ----sir-pfilter-2,results='markup',cache=T------------------------------
pf <- replicate(10,pfilter(sir,Np=5000,params=c(Beta=2,gamma=1,rho=0.8,N=2600)))
ll <- sapply(pf,logLik); ll
logmeanexp(ll,se=TRUE)

## ----sir-like-slice,cache=TRUE,results='hide'----------------------------
sliceDesign(
  c(Beta=2,gamma=1,rho=0.8,N=2600),
  Beta=rep(seq(from=0.5,to=4,length=40),each=3),
  gamma=rep(seq(from=0.5,to=2,length=40),each=3)) -> p

require(foreach)
require(doMC)

registerDoMC(cores=5)        
## number of cores 
## usually the number of cores on your machine, or slightly smaller

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)

foreach (theta=iter(p,"row"),.combine=rbind,
         .inorder=FALSE,.options.multicore=mcopts) %dopar% 
 {
   pfilter(sir,params=unlist(theta),Np=5000) -> pf
   theta$loglik <- logLik(pf)
   theta
 } -> p

## ----sir-like-slice-plot,cache=F,results="hide"--------------------------
foreach (v=c("Beta","gamma")) %do% 
{
  x <- subset(p,slice==v)
  plot(x[[v]],x$loglik,xlab=v,ylab="loglik")
}

## ----sir-grid1-----------------------------------------------------------
expand.grid(Beta=seq(from=1,to=4,length=50),
            gamma=seq(from=0.7,to=3,length=50),
            rho=0.8,
            N=2600) -> p

foreach (theta=iter(p,"row"),.combine=rbind,
         .inorder=FALSE,.options.multicore=mcopts) %dopar% 
 {
   pfilter(sir,params=unlist(theta),Np=5000) -> pf
   theta$loglik <- logLik(pf)
   theta
 } -> p


## ----sir-grid1-plot------------------------------------------------------
pp <- mutate(p,loglik=ifelse(loglik>max(loglik)-100,loglik,NA))
ggplot(data=pp,mapping=aes(x=Beta,y=gamma,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  geom_contour(color='black',binwidth=3)+
  scale_fill_gradient()+
  labs(x=expression(beta),y=expression(gamma))

## ----load_bbs------------------------------------------------------------
bsflu_data <- read.table("bsflu_data.txt")

## ----bsflu_names---------------------------------------------------------
bsflu_statenames <- c("S","I","R1","R2")
bsflu_paramnames <- c("Beta","mu_I","rho","mu_R1","mu_R2")

## ----bsflu_obsnames------------------------------------------------------
(bsflu_obsnames <- colnames(bsflu_data)[1:2])

## ----csnippets_bsflu-----------------------------------------------------
bsflu_dmeasure <- "
  lik = dpois(B,rho*R1+1e-6,give_log);
"

bsflu_rmeasure <- "
  B = rpois(rho*R1+1e-6);
  C = rpois(rho*R2);
"

bsflu_rprocess <- "
  double t1 = rbinom(S,1-exp(-Beta*I*dt));
  double t2 = rbinom(I,1-exp(-dt*mu_I));
  double t3 = rbinom(R1,1-exp(-dt*mu_R1));
  double t4 = rbinom(R2,1-exp(-dt*mu_R2));
  S -= t1;
  I += t1 - t2;
  R1 += t2 - t3;
  R2 += t3 - t4;
"

bsflu_fromEstimationScale <- "
 TBeta = exp(Beta);
 Tmu_I = exp(mu_I);
 Trho = expit(rho);
"

bsflu_toEstimationScale <- "
 TBeta = log(Beta);
 Tmu_I = log(mu_I);
 Trho = logit(rho);
"

bsflu_initializer <- "
 S=762;
 I=1;
 R1=0;
 R2=0;
"

## ----pomp_bsflu----------------------------------------------------------
require(pomp)
stopifnot(packageVersion("pomp")>="0.75-1")
bsflu2 <- pomp(
  data=bsflu_data,
  times="day",
  t0=0,
  rprocess=euler.sim(
    step.fun=Csnippet(bsflu_rprocess),
    delta.t=1/12
  ),
  rmeasure=Csnippet(bsflu_rmeasure),
  dmeasure=Csnippet(bsflu_dmeasure),
  fromEstimationScale=Csnippet(bsflu_fromEstimationScale),
  toEstimationScale=Csnippet(bsflu_toEstimationScale),
  obsnames = bsflu_obsnames,
  statenames=bsflu_statenames,
  paramnames=bsflu_paramnames,
  initializer=Csnippet(bsflu_initializer)
)
plot(bsflu2)

## ----run_level-----------------------------------------------------------
run_level <- 3
switch(run_level,
       {bsflu_Np=100; bsflu_Nmif=10; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10}, 
       {bsflu_Np=20000; bsflu_Nmif=100; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10}, 
       {bsflu_Np=60000; bsflu_Nmif=300; bsflu_Neval=10; bsflu_Nglobal=100; bsflu_Nlocal=20}
)

## ----bsflu_params--------------------------------------------------------
bsflu_params <- data.matrix(read.table("mif_bsflu_params.csv",row.names=NULL,header=TRUE))
bsflu_mle <- bsflu_params[which.max(bsflu_params[,"logLik"]),][bsflu_paramnames]

## ----fixed_params--------------------------------------------------------
bsflu_fixed_params <- c(mu_R1=1/(sum(bsflu_data$B)/512),mu_R2=1/(sum(bsflu_data$C)/512))

## ----parallel-setup,cache=FALSE------------------------------------------
require(doParallel)
cores <- 20  # The number of cores on this machine 
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)

set.seed(396658101,kind="L'Ecuyer")

## ----alternative_parallel, eval=FALSE------------------------------------
## require(doMC)
## registerDoMC(cores=20)

## ----pf------------------------------------------------------------------
stew(file=sprintf("pf-%d.rda",run_level),{
  t_pf <- system.time(
    pf <- foreach(i=1:20,.packages='pomp',
                  .options.multicore=mcopts) %dopar% try(
                    pfilter(bsflu2,params=bsflu_mle,Np=bsflu_Np)
                  )
  )
  
},seed=1320290398,kind="L'Ecuyer")

(L_pf <- logmeanexp(sapply(pf,logLik),se=TRUE))

## ----set_cache, eval=FALSE-----------------------------------------------
## opts_chunk$set(
##   cache=TRUE,
##   )

## ----box_search_local----------------------------------------------------
bsflu_rw.sd <- 0.02
bsflu_cooling.fraction.50 <- 0.5

stew(file=sprintf("local_search-%d.rda",run_level),{
  
  t_local <- system.time({
    mifs_local <- foreach(i=1:bsflu_Nlocal,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  {
      mif2(
        bsflu2,
        start=bsflu_mle,
        Np=bsflu_Np,
        Nmif=bsflu_Nmif,
        cooling.type="geometric",
        cooling.fraction.50=bsflu_cooling.fraction.50,
        transform=TRUE,
        rw.sd=rw.sd(
          Beta=bsflu_rw.sd,
          mu_I=bsflu_rw.sd,
          rho=bsflu_rw.sd
        )
      )
      
    }
  })
  
},seed=900242057,kind="L'Ecuyer")


## ----lik_local_eval------------------------------------------------------
stew(file=sprintf("lik_local-%d.rda",run_level),{
    t_local_eval <- system.time({
    liks_local <- foreach(i=1:bsflu_Nlocal,.packages='pomp',.combine=rbind) %dopar% {
      evals <- replicate(bsflu_Neval, logLik(pfilter(bsflu2,params=coef(mifs_local[[i]]),Np=bsflu_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")

results_local <- data.frame(logLik=liks_local[,1],logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))
summary(results_local$logLik,digits=5)

## ----pairs_local---------------------------------------------------------
pairs(~logLik+Beta+mu_I+rho,data=subset(results_local,logLik>max(logLik)-50))

## ----box-----------------------------------------------------------------
bsflu_box <- rbind(
  Beta=c(0.001,0.01),
  mu_I=c(0.5,2),
  rho = c(0.5,1)
)

## ----box_eval------------------------------------------------------------
stew(file=sprintf("box_eval-%d.rda",run_level),{
  
  t_global <- system.time({
    mifs_global <- foreach(i=1:bsflu_Nglobal,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  mif2(
      mifs_local[[1]],
      start=c(apply(bsflu_box,1,function(x)runif(1,x[1],x[2])),bsflu_fixed_params)
    )
  })
},seed=1270401374,kind="L'Ecuyer")

## ----lik_global_eval-----------------------------------------------------
stew(file=sprintf("lik_global_eval-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:bsflu_Nglobal,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(bsflu_Neval, logLik(pfilter(bsflu2,params=coef(mifs_global[[i]]),Np=bsflu_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

results_global <- data.frame(logLik=liks_global[,1],logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))
summary(results_global$logLik,digits=5)

## ----save_params,eval=FALSE----------------------------------------------
## if (run_level>2)
##   write.table(rbind(results_local,results_global),
##               file="mif_bsflu_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)

## ----pairs_global--------------------------------------------------------
pairs(~logLik+Beta+mu_I+rho,data=subset(results_global,logLik>max(logLik)-250))

