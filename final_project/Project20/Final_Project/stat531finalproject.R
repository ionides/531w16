#setwd("C:\\Users\\User\\Google Drive\\Winter 2016\\STAT 531 - Time Series Analysis\\Final Project\\R Code")

#load libraries
require(ggplot2)
require(plyr)
require(reshape2)
require(magrittr)
require(pomp)
require(foreach)

set.seed(123456)

delhi <- read.csv("delhi.csv",head=T)

##---plot data----------------------------------------------------
month <- delhi$month
cases <- delhi$cases
plot(month,cases,type='o',
     main="Monthly Reported TB Cases in Delhi (2007-2012)",
     xlab="Months",ylab="No. of Reported Cases")

##---r process----------------------------------------------------
sir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-lambda*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-gamma*I*dt));
  double dN_IQ = rbinom(I,1-exp(-alpha*Q*dt));
  double dN_QR = rbinom(Q,1-exp(-eta*R*dt));
  double dN_SV = rbinom(S,1-exp(-sigma*V*dt));
  double dN_VS = rbinom(V,1-exp(-rho*S*dt));
  double dN_RS = rbinom(R,1-exp(-epsilon*S*dt));
  double dN_IR = rbinom(I,1-exp(-phi*R*dt));
  S -= dN_SE - dN_SV + dN_VS + dN_RS;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IQ - dN_IR;
  Q += dN_IQ - dN_QR;
  R += dN_QR - dN_RS + dN_IR;
  V += dN_SV - dN_VS;
  H += dN_EI;
")

##---initialize variables in r process----------------------------
sir_init <- Csnippet("
  S = N-1;
  E = 1;
  I = 0;
  Q = 0;
  R = 0;
  V = 0;
  H = 0;
")

##---create pomp object
seirqv <- pomp(delhi,times = "month", t0=1,
               rprocess=euler.sim(sir_step,delta.t=1/2),
               initializer=sir_init,
               paramnames=c("lambda","gamma","alpha","eta",
                            "sigma","rho","epsilon","phi","N"),
               statenames=c("S","E","I","Q","R","V","H"))

##---add dmeasure and r measure into pomp object------------------
dmeas <- Csnippet("lik = dbinom(cases,H,psi,give_log);")
rmeas <- Csnippet("cases = rbinom(H,psi);")
seirqv <- pomp(seirqv,rmeasure=rmeas,dmeasure=dmeas,
               statenames="H",zeronames="H",paramnames="psi")

##---simulate from pomp model we created--------------------------
sims <- simulate(seirqv,
                 params=c(N=64000,psi=0.60,
                          lambda=38400,gamma=0.3,alpha=0.1,eta=0.4,
                          sigma=0.3,rho=0.6,epsilon=0.3,phi=0.68),
                 nsim=1000,as=TRUE,include=TRUE)

ggplot(sims,aes(x=time,y=cases,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)

##---evaluate log likelihood--------------------------------------
p <- c(N=64000,psi=0.60,
       lambda=38400,gamma=0.3,alpha=0.1,eta=0.4,
       sigma=0.3,rho=0.6,epsilon=0.3,phi=0.68)
t <- c(1:72)
t0 <- 0
x0 <- init.state(seirqv,params=p,t0=t0)
x <- rprocess(seirqv,xstart=x0,times=c(t0,t),params=p,offset=1)
y <- rmeasure(seirqv,params=p,x=x,times=t)
ll <- dmeasure(seirqv,y=y,x=x,times=t,params=p,log=TRUE)
ell <- apply(ll,1,sum)
summary(exp(ell))
logmeanexp(ell,se=TRUE)

##---simulate and plot the simulated parameters-------------------
simulate(seirqv,params=c(N=64000,psi=0.90,
                         lambda=1,gamma=0.5,alpha=0.3,eta=1,
                         sigma=1,rho=1,epsilon=1,phi=0.5),
         nsim=10000,states=TRUE) -> x
matplot(time(seirqv),t(x["H",1:50,]),type='l',lty=1,
        xlab="time",ylab="H",bty='l',col='blue')
lines(time(seirqv),obs(seirqv,"cases"),lwd=2,col='black')

matplot(time(seirqv),t(x["psi",1:50,]),type='l',lty=1,
        xlab="time",ylab="psi",bty='l',col='blue')
lines(time(seirqv),obs(seirqv,"cases"),lwd=2,col='black')

##---estimate parameters of pomp model----------------------------
cores <- 2  # The number of cores on this machine 
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)
delhi_box <- rbind(
  psi=c(0.1,0.90),
  lambda=c(0.5,3),
  gamma=c(0.1,0.90),
  alpha=c(0.1,0.90),
  phi=c(0.1,0.90),
  epsilon = c(0.1,0.90),
  eta=c(0.1,0.90)
)
delhi_rw.sd <- 0.02

run_level <- 2
switch(run_level,
       {delhi_Np=100; delhi_Nmif=10; delhi_Neval=10; delhi_Nglobal=10; delhi_Nlocal=10}, 
       {delhi_Np=20000; delhi_Nmif=100; delhi_Neval=10; delhi_Nglobal=10; delhi_Nlocal=10}, 
       {delhi_Np=60000; delhi_Nmif=300; delhi_Neval=10; delhi_Nglobal=100; delhi_Nlocal=20},
       {delhi_Np=50000; delhi_Nmif=280; delhi_Neval=10; delhi_Nglobal=100; delhi_Nlocal=20}
)
delhi_cooling.fraction.50 <- 0.6
delhi_fixed_params <- c(N=64000,sigma=1,rho=1)

stew(file=sprintf("box_eval-%d.rda",run_level),{
  t_global <- system.time({
    mifs_global <- foreach(i=1:delhi_Nglobal,.packages='pomp', 
                           .combine=c, .options.multicore=mcopts) %dopar% 
      mif2(     
        seirqv,
        start=c(apply(delhi_box,1,function(x)runif(1,x[1],x[2])),delhi_fixed_params),
        Np=delhi_Np,
        Nmif=delhi_Nmif,
        cooling.type="geometric",
        cooling.fraction.50=delhi_cooling.fraction.50,
        transform=TRUE,
        rw.sd=rw.sd(
          psi=delhi_rw.sd,
          lambda=delhi_rw.sd,
          gamma=delhi_rw.sd,
          alpha=delhi_rw.sd,
          phi=delhi_rw.sd,
          epsilon = delhi_rw.sd,
          eta=delhi_rw.sd
        )
      )
  })
},seed=1270401374,kind="L'Ecuyer")

stew(file=sprintf("lik_global_eval-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:delhi_Nglobal,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(delhi_Neval, logLik(pfilter(seirqv,params=coef(mifs_global[[i]]),Np=delhi_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")





