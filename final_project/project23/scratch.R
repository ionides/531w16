require(ggplot2)
require(plyr)
require(reshape2)
require(magrittr)
require(pomp)
require(foreach)
delhi <- read.csv("delhi.csv",head=T)

##---r process----------------------------------------------------
sir_step <- Csnippet("
                     double dN_SL = rbinom(S,1-exp(-(1-k)*dt));
                     double dN_SA = rbinom(S,1-exp(-(k)*dt));
                     double dN_LA = rbinom(L,1-exp(-(x2)*dt));
                     double dN_RA = rbinom(R,1-exp(-(x7+c*lambda)*dt));
                     double dN_AO = rbinom(A,1-exp(-(x3-p*x3)*dt));
                     double dN_AB = rbinom(A,1-exp(-(p*x3)*dt));
                     double dN_EB = rbinom(E,1-exp(-(x5)*dt));
                     double dN_BE = rbinom(B,1-exp(-(1-phi)*dt));
                     double dN_BC = rbinom(B,1-exp(-(phi)*dt));
                     double dN_CR = rbinom(C,1-exp(-(alpha)*dt));
                     double dN_CE = rbinom(C,1-exp(-(1-alpha)*dt));
                     S -= dN_SL-dN_SA-mu;
                     L += dN_SL-dN_LA-mu;
                     A += dN_LA+dN_SA+dN_RA-dN_AO-dN_AB-mu-mud;
                     O += dN_AO-mu-mud;
                     B += dN_AB+dN_EB-dN_BE-dN_BC-mu-mud;
                     C += dN_BC-dN_CR-dN_CE-mu-mud;
                     R += dN_CR-dN_RA-mu;
                     E += dN_BE+dN_CE-dN_EB-mu-mud;
                     lambda += beta*(A+O+B+C+E);
                     H += dN_LA;
                     ")

##---initialize variables in r process----------------------------
sir_init <- Csnippet("
  S = N-1;
  L = 1;
  A = 0;
  O = 0;
  B = 0;
  C = 0;
  R = 0;
  E = 0;
  H = 0;
  lambda = 0;
")

model <-  pomp(delhi,times = "month", t0=1,
               rprocess=euler.sim(sir_step,delta.t=1/12),
               initializer=sir_init,
               paramnames=c("N","k","x2","x7","c","x3","p",
                            "x5","phi","alpha","mu","mud","beta"),
               statenames=c("S","L","A","O","B","C","R","E","H","lambda"))

##---add dmeasure and r measure into pomp object------------------
dmeas <- Csnippet("lik = dbinom(cases,H,psi,give_log);")
rmeas <- Csnippet("cases = rbinom(H,psi);")
model <- pomp(model,rmeasure=rmeas,dmeasure=dmeas,
              statenames="H",zeronames="H",paramnames="psi")

##---simulate from pomp model we created--------------------------
sims <- simulate(model,
                 params=c(N=64000,k=0.1,x2=0.01/12,x7=0.0033/12,c=0.01,
                          x3=12/12,p=0.9,x5=12/12,phi=0.13,alpha=0.58,
                          mu=0.0151/12,mud=0.24/12,beta=11.03,psi=0.54),
                 nsim=1000,as=TRUE,include=TRUE,states=TRUE)

ggplot(sims,aes(x=time,y=cases,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)

##---evaluate log likelihood--------------------------------------
p <- c(N=64000,k=0.1,x2=0.01/12,x7=0.0033/12,c=0.01,
       x3=12/12,p=0.9,x5=12/12,phi=0.13,alpha=0.58,
       mu=0.0151/12,mud=0.24/12,beta=11.03,psi=0.54)
t <- c(1:72)
t0 <- 0
x0 <- init.state(model,params=p,t0=t0)
x <- rprocess(model,xstart=x0,times=c(t0,t),params=p,offset=1)
y <- rmeasure(model,params=p,x=x,times=t)
ll <- dmeasure(model,y=y,x=x,times=t,params=p,log=TRUE)
ell <- apply(ll,1,sum); summary(exp(ll)); logmeanexp(ll,se=TRUE)





