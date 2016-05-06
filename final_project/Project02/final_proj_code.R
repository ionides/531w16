##--  Set Working Directory --##
## setwd("C:/Users/zkeller/Google Drive/stats_531/final_project")
library(pomp)
library(dplyr)
library(magrittr)
library(plyr)
library(ggplot2)
library(reshape2)
##--  Read In Data --##
mca <- read.csv("measles_ca.csv", skip=2, header=T, na.strings="-")
names(mca) <- c("year", "week", "cases")

##--  Format Data --##
##--  Remove 0s that occur between spikes --##
require(dplyr)
mca <- mca %>%
  filter(year >= 1989 & year <= 1991) %>%
  mutate(time = year +  week/52) %>%
  select(time,cases)

orig <- mca

mca.over.lim <- mca$cases >= 100
mca$cases[c(mca.over.lim[2:nrow(mca)], F) & c(F, mca.over.lim[1:(nrow(mca)-1)]) & (mca$cases <= 1)] <- NA
mca <- mca[complete.cases(mca),]
##--  Plot Data --##
plot(orig$time, orig$cases, type="l", xaxt = "n", ylab="Cases", xlab="Time", main="California Measles Outbreak \n Removed Data in Black")
axis(1, at=c(1989,1990,1991,1992))
lines(mca$time, mca$cases, col="red")

##--  Get Covariates --##

##
## Population
##

read.census.data.pre91 <- function(dat.url, yr){
  pop.data <- read.table(dat.url, skip=29, header=F)
  names(pop.data) <- c("FIPS", "State","Age", "Pop", "Male", "Female")
  pop.data <- pop.data %>% 
    select(State, Pop, Age) %>%
    filter(State == "CA", Age<=9) %>%
    group_by(State) %>%
    summarize(pop = sum(Pop)) %>%
    mutate(year = yr)
  
  return(pop.data)
}

# pop87 <- read.census.data.pre91("www.census.gov/popest/data/state/asrh/1980s/tables/stiag787.txt", 1987)
# pop88 <- read.census.data.pre91("www.census.gov/popest/data/state/asrh/1980s/tables/stiag788.txt",1988)
# pop89 <- read.census.data.pre91("www.census.gov/popest/data/state/asrh/1980s/tables/stiag789.txt",1989)
# pop90 <- read.census.data.pre91("www.census.gov/popest/data/state/asrh/1980s/tables/stiag490.txt", 1990)

pop87 <- read.table("pop87.txt")
pop88 <- read.table("pop88.txt")
pop89 <- read.table("pop89.txt")
pop90 <- read.table("pop90.txt")

## next data from "https://www.census.gov/popest/data/state/asrh/1990s/tables/ST-99-08.txt"

pop91 <- c(2664214+2296748, 1991)
pop92 <- c(2752513+2314709, 1992)
pop93 <- c(2807471+2344063, 1993)
tot.pop <- rbind(pop87, pop88, pop89, pop90, pop91, pop92, pop93)



##
## Births
##

## Data from: http://www.cdph.ca.gov/data/statistics/Documents/VSC-2005-0201.pdf

births <- as.data.frame(list("births"=c(503376, 532708, 569308, 611666, 609228, 600838, 584483),
                             "year"=c(1987, 1988, 1989, 1990, 1991, 1992, 1993)))

## smoothing
demog <- inner_join(tot.pop, births, by=c("year"))

demog %>% 
  summarize(
    time=seq(from=min(year),to=max(year),by=1/52),
    pop=predict(smooth.spline(x=year,y=pop),x=time)$y,
    birthrate=predict(smooth.spline(x=year,y=births),x=time)$y # Note +0.5 from Aaron King's Code
      ) -> covar

birthrate_lag1 <- covar %>% select(time, birthrate) %>% mutate(time = time+1) %>% dplyr::rename(birthrate_lag1=birthrate) 
birthrate_lag05 <- covar %>% select(time, birthrate) %>% mutate(time = time+0.5) %>% dplyr::rename(birthrate_lag05=birthrate)

covar <- inner_join(covar, birthrate_lag1, by=c("time")) %>% inner_join(birthrate_lag05, by=c("time")) %>% filter(time >= 1989 & time <= 1992)

 
plot(pop~time, data=covar, type="l")
points(pop~year, data=demog)

plot(birthrate~time, data=covar, type="l")
points(births~year, data=demog)
  

##--  Our Compartment Model  --##
library(DiagrammeR)
DiagrammeR("digraph SEIR {
           graph [rankdir=TD, overlap=false, fontsize = 10]
           node[shape=egg, label='B'] b;
           subgraph {
           rank=same;
           node[shape=oval, label='S'] S;
           node[shape=oval, label='E'] E;
           node[shape=oval, label='I'] I;
           node[shape=oval, label='R'] R;
           S->E E->I I->R
           }
           node[shape=diamond, label='dead'] d;
           b->S
           {S E I R}->d
           }",type="grViz",engine="dot",height=300,width=800)


##-- Our Process Model  --##
## From http://kingaa.github.io/sbied/measles/measles.html ##

## Should we remove cohort effect and seasonality here?
## Vaccination Coverage
## http://apps.who.int/gho/data/node.main.A826 Vaccination Coverage
rproc <- Csnippet("
                  double beta, br, seas, foi, dw, births, vac;
                  double rate[6], trans[6];
                  
                  // cohort effect
                  // if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt) 
                  // br = cohort*birthrate_lag05/dt + (1-cohort)*birthrate_lag05;
                  // else 
                  // br = (1.0-cohort)*birthrate_lag05;
                  
                  // term-time seasonality
                  t = (t-floor(t))*365.25;
                  if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
                  seas = 1.0+amplitude*0.2411/0.7589;
                  else
                  seas = 1.0-amplitude;
                  
                  // transmission rate
                  beta = R0*(gamma+mu)*seas;
                  // expected force of infection
                  foi = beta*pow(I+iota,alpha)/pop;
                  // white noise (extrademographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  rate[0] = foi*dw/dt;  // stochastic force of infection
                  rate[1] = mu;			    // natural S death
                  rate[2] = sigma;		  // rate of ending of latent stage
                  rate[3] = mu;			    // natural E death
                  rate[4] = gamma;		  // recovery
                  rate[5] = mu;			    // natural I death
                  
                  // Poisson births
                  births = rpois(birthrate_lag05*dt);

                  // Vaccination
                  vac = nearbyint(vr*birthrate_lag1*.95*dt);
                  // printf(\" vr is: %f, birthrate_lag is: %f, dt is: %f , vac is: %f \", vr, birthrate_lag1, dt, vac);

                  // transitions between classes
                  reulermultinom(2,S,&rate[0],dt,&trans[0]);
                  reulermultinom(2,E,&rate[2],dt,&trans[2]);
                  reulermultinom(2,I,&rate[4],dt,&trans[4]);
                  // printf(\" t2 is: %f, t4 is: %f, t5 is: %f \", trans[2],trans[4],trans[5]);
                  
		  if (vac > S - trans[0] - trans[1]){
		  	vac = S - trans[0] - trans[1];
		  }

                  S += births   - trans[0] - trans[1] - vac;
                  // printf(\" births are: %f,  S is: %f \",births, S);
                  E += trans[0] - trans[2] - trans[3];
                  I += trans[2] - trans[4] - trans[5];
                  R = pop - S - E - I + vac;
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[4];           // true incidence
                  ")

initlz <- Csnippet("
  double m = pop/(S_0+E_0+I_0+R_0);
                   S = nearbyint(m*S_0);
                   E = nearbyint(m*E_0);
                   I = nearbyint(m*I_0);
                   R = nearbyint(m*R_0);
                   W = 0;
                   C = 0;
                   ")

##  Measurement Model  ##
dmeas <- Csnippet("
  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  if (cases > 0.0) {
                  // printf(\"%f \",C);
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  // printf(\" %f \",C);
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
                  ")

rmeas <- Csnippet("
  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  cases = rnorm(m,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")

##  Parameter Transformations  ##

## ----transforms----------------------------------------------------------
toEst <- Csnippet("
                  Tmu = log(mu);
                  Tsigma = log(sigma);
                  Tgamma = log(gamma);
                  Talpha = log(alpha);
                  Tiota = log(iota);
                  Trho = logit(rho);
                  // Tcohort = logit(cohort);
                  Tamplitude = logit(amplitude);
                  TsigmaSE = log(sigmaSE);
                  Tpsi = log(psi);
                  TR0 = log(R0);
                  Tvr = logit(vr);
                  to_log_barycentric (&TS_0, &S_0, 4);
                  ")

fromEst <- Csnippet("
                    Tmu = exp(mu);
                    Tsigma = exp(sigma);
                    Tgamma = exp(gamma);
                    Talpha = exp(alpha);
                    Tiota = exp(iota);
                    Trho = expit(rho);
                    // Tcohort = expit(cohort);
                    Tamplitude = expit(amplitude);
                    TsigmaSE = exp(sigmaSE);
                    Tpsi = exp(psi);
                    TR0 = exp(R0);
                    Tvr = expit(vr);
                    from_log_barycentric (&TS_0, &S_0, 4);
                    ")

mca %>% 
  pomp(t0=with(mca, time[1]),
       time="time",
       # params=theta,
       rprocess=euler.sim(rproc,delta.t=1/365.25),
       initializer=initlz,
       dmeasure=dmeas,
       rmeasure=rmeas,
       toEstimationScale=toEst,
       fromEstimationScale=fromEst,
       covar=covar,
       tcovar="time",
       zeronames=c("C","W"),
       statenames=c("S","E","I","R","C","W"),
       paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                    "rho","sigmaSE","psi", "amplitude",# "cohort", 
                    "S_0","E_0","I_0","R_0","vr")
  ) -> m1

m1 %>% as.data.frame() %>% 
  melt(id="time") %>%
  ggplot(aes(x=time,y=value))+
  geom_line()+
  facet_grid(variable~.,scales="free_y")

##--  Note we need to estimate parameters --##

run_level=3
measles_Np <-          c(1000,5e3,1e4)
measles_Nmif <-        c(10, 200,400)
measles_Neval <-  c(2,  10,  20)
measles_Nlocal <- c(10, 20, 40)
measles_Nglobal <-c(10, 20, 100)
measles_Nsim <-        c(50,100, 500)

##-- Test particle filter on measles paramters computed in He. et al  --##

## ----mles,include=FALSE--------------------------------------------------
read.csv(text=
         "town,loglik,loglik.sd,mu,delay,sigma,gamma,rho,R0,amplitude,alpha,iota,cohort,psi,S_0,E_0,I_0,R_0,sigmaSE
         Bedwellty,-1125.1,0.14,0.02,4,57.9,146,0.311,24.7,0.16,0.937,0.0396,0.351,0.951,0.0396,2.64e-05,2.45e-05,0.96,0.0611
         Birmingham,-3239.3,1.55,0.02,4,45.6,32.9,0.544,43.4,0.428,1.01,0.343,0.331,0.178,0.0264,8.96e-05,0.000335,0.973,0.0611
         Bradford,-2586.6,0.68,0.02,4,45.6,129,0.599,32.1,0.236,0.991,0.244,0.297,0.19,0.0365,7.41e-06,4.59e-06,0.964,0.0451
         Bristol,-2681.6,0.5,0.02,4,64.3,82.6,0.626,26.8,0.203,1.01,0.441,0.344,0.201,0.0358,9.62e-06,5.37e-06,0.964,0.0392
         Cardiff,-2364.9,0.73,0.02,4,39,143,0.602,34.4,0.223,0.996,0.141,0.267,0.27,0.0317,1.01e-05,9.21e-06,0.968,0.0539
         Consett,-1362.9,0.73,0.02,4,42.6,172,0.65,35.9,0.2,1.01,0.0731,0.31,0.406,0.0322,1.83e-05,1.97e-05,0.968,0.0712
         Dalton.in.Furness,-726.1,0.3,0.02,4,73.6,257,0.455,28.3,0.203,0.989,0.0386,0.421,0.818,0.0387,2.23e-05,2.36e-05,0.961,0.0779
         Halesworth,-318.6,0.51,0.02,4,49.6,210,0.754,33.1,0.381,0.948,0.00912,0.547,0.641,0.0526,1.99e-05,2.82e-05,0.947,0.0748
         Hastings,-1583.7,0.21,0.02,4,56.3,74.1,0.695,34.2,0.299,1,0.186,0.329,0.396,0.0233,5.61e-06,3.4e-06,0.977,0.0955
         Hull,-2729.4,0.39,0.02,4,42.1,73.9,0.582,38.9,0.221,0.968,0.142,0.275,0.256,0.0371,1.2e-05,1.13e-05,0.963,0.0636
         Leeds,-2918.6,0.23,0.02,4,40.7,35.1,0.666,47.8,0.267,1,1.25,0.592,0.167,0.0262,6.04e-05,3e-05,0.974,0.0778
         Lees,-548.1,1.1,0.02,4,45.6,244,0.612,29.7,0.153,0.968,0.0311,0.648,0.681,0.0477,2.66e-05,2.08e-05,0.952,0.0802
         Liverpool,-3403.1,0.34,0.02,4,49.4,39.3,0.494,48.1,0.305,0.978,0.263,0.191,0.136,0.0286,0.000184,0.00124,0.97,0.0533
         London,-3804.9,0.16,0.02,4,28.9,30.4,0.488,56.8,0.554,0.976,2.9,0.557,0.116,0.0297,5.17e-05,5.14e-05,0.97,0.0878
         Manchester,-3250.9,0.66,0.02,4,34.4,56.8,0.55,32.9,0.29,0.965,0.59,0.362,0.161,0.0489,2.41e-05,3.38e-05,0.951,0.0551
         Mold,-296.5,0.25,0.02,4,67.4,301,0.131,21.4,0.271,1.04,0.0145,0.436,2.87,0.064,2.61e-05,2.27e-05,0.936,0.0544
         Northwich,-1195.1,2.25,0.02,4,45.6,147,0.795,30.1,0.423,0.948,0.0602,0.236,0.402,0.0213,1.32e-05,1.58e-05,0.979,0.0857
         Nottingham,-2703.5,0.53,0.02,4,70.2,115,0.609,22.6,0.157,0.982,0.17,0.34,0.258,0.05,1.36e-05,1.41e-05,0.95,0.038
         Oswestry,-696.1,0.49,0.02,4,37.3,168,0.631,52.9,0.339,1.04,0.0298,0.263,0.476,0.0218,1.56e-05,1.61e-05,0.978,0.0699
         Sheffield,-2810.7,0.21,0.02,4,54.3,62.2,0.649,33.1,0.313,1.02,0.853,0.225,0.175,0.0291,6.04e-05,8.86e-05,0.971,0.0428
         ",stringsAsFactors=FALSE,strip.white=T) -> mles

## ----mle-----------------------------------------------------------------
paramnames <- c("R0","mu","sigma","gamma","alpha","iota",
                "rho","sigmaSE","psi","amplitude",#"cohort",
                "S_0","E_0","I_0","R_0", "vr")

mles$S_0 <- 0.03
mles$vr <- 0.8

# Data from https://www.cdph.ca.gov/data/statistics/Documents/VSC-2005-0101.pdf
#   gives us about 7 deaths per 1 thousand population
mles$mu <- 0.007 # death rate
mles$delay <- 0.5 # delay from entry into succeptibles: 
# mles$sigma <- 26 # rate of latent stage
##--  Carry out Particle Filtering to test that it works  --##
##--  Also use this step to select our intial parameters  --##
##--  We choose the parameters with the highest likelihood
##--    from those given by He et all --##
require(doParallel)
require(foreach)

mcopts <- list(set.seed=TRUE)
set.seed(396658101,kind="L'Ecuyer")

cl <- makeCluster(12)
registerDoParallel(cl)
clusterExport(cl,c("m1", "measles_Np", "run_level"))

stew(file=sprintf("pf1-%d.rda",run_level),{
    param_mles <- as.data.frame(list(town=vector("character", nrow(mles)),
                                     ell=vector("numeric", nrow(mles)),
                                     se=vector("numeric", nrow(mles))), stringsAsFactors = F)
    for(i in 1:nrow(mles)){
      temp.params <- unlist(extract(mles[i,], paramnames))
      temp.pf <- foreach(j=1:20,.packages='pomp',#.export=c("m1", "theta", "measles_Np", "run_level"),
                         .options.multicore=mcopts) %dopar% try(
                           pfilter(m1,params=temp.params, Np=measles_Np[run_level], save.states=T)
                         )
      temp.L <- logmeanexp(sapply(temp.pf, logLik), se=T)
      param_mles[i,1] <- mles[i,1]
      param_mles[i,2] <- temp.L[1]
      param_mles[i,3] <- temp.L[2]
      }
},seed=493536993,kind="L'Ecuyer")
stopCluster(cl)

##-- find argmax of our param_mles and extract those parameters --##
require(ggplot2)
ggplot(param_mles, aes(town, ell, ymin=ell-1.96*se, ymax=ell+1.96*se, color=town)) + 
  geom_errorbar() +
  xlab("Town's") +
  ylab("Log Likelihood Interval") + 
  ggtitle("Initial Parameter Likelihood for Ca Measles") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.2)) + 
  theme(legend.position="none")

##--  We go with Bedwellty  --##
init_param <- unlist(extract(mles[which(param_mles$town =="Bedwellty"),], paramnames))

m1 <- pomp(m1, params=init_param)
##--  Continue with iterated filtering to get better parameters --##
##--  This is a local search  --##

measles_rw.sd <- 0.02
measles_rw.sd.R0 <- 0.1
measles_rw.sd.simga <- 0.1
measles_cooling.fraction.50 <- 0.5

cl <- makeCluster(12)
registerDoParallel(cl)
clusterExport(cl,c("m1", 
                   "measles_Np",
                   "init_param",
                   "measles_Nmif",
                   "measles_Nlocal",
                   "run_level",
                   "measles_cooling.fraction.50",
                   "measles_rw.sd",
                   "measles_rw.sd.R0",
                   "measles_rw.sd.simga"))

stew(file=sprintf("local_search-%d.rda",run_level),{
  
  t_local <- system.time({
    mifs_local <- foreach(i=1:measles_Nlocal[run_level],
                          .packages='pomp',
                          .combine=c,
                          .options.multicore=mcopts) %dopar%  {
      mif2(
        m1,
        start=init_param,
        Np=measles_Np[run_level],
        Nmif= measles_Nmif[run_level],
        cooling.type="geometric",
        cooling.fraction.50=measles_cooling.fraction.50,
        transform=TRUE,
        rw.sd=rw.sd(
          R0=measles_rw.sd.R0,
          gamma=measles_rw.sd,
          alpha=measles_rw.sd,
          iota=measles_rw.sd,
          rho=measles_rw.sd,
          psi=measles_rw.sd,
          sigma=measles_rw.sd.simga,
          sigmaSE=measles_rw.sd,
          vr=measles_rw.sd,
          amplitude=measles_rw.sd
          # S_0=measles_rw.sd,
          # E_0=measles_rw.sd
        )
      )
      
    }
  })
  
},seed=900242057,kind="L'Ecuyer")

stopCluster(cl)


##--  Evaluate this local likelihood  --##
cl <- makeCluster(12)
registerDoParallel(cl)
clusterExport(cl, c("measles_Neval", "run_level", "m1", "mifs_local","measles_Np"))

stew(file=sprintf("lik_local-%d.rda",run_level),{
  t_local_eval <- system.time({
    liks_local <- foreach(i=1:measles_Nlocal[run_level],.packages='pomp',.combine=rbind) %dopar% {
      evals <- replicate(measles_Neval[run_level],
                         logLik(pfilter(m1,params=coef(mifs_local[[i]]),Np=measles_Np[run_level])))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")

stopCluster(cl)

results_local <- data.frame(logLik=liks_local[,1],logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))
summary(results_local$logLik,digits=5)


##--  Global Search --##

measles_box <- rbind(
  R0=c(29,40),
  gamma=c(60,169),
  alpha=c(.7,1),
  iota=c(0,.4),
  rho=c(.15,.65),
  psi=c(.15,.5),
  sigma=c(41,56),
  sigmaSE=c(.03,.09),
  vr=c(0.85,.95),
  amplitude=c(0.2, 0.5)
)

measles_fixed_params = c(mu=.007, S_0=.00477, E_0=2.66e-05, I_0=2.081e-05, R_0=.9522)

cl <- makeCluster(12)
registerDoParallel(cl)
clusterExport(cl, c("measles_Neval", "run_level", "m1", "mifs_local","measles_Np", "measles_box", "measles_fixed_params"))

stew(file=sprintf("box_eval-%d.rda",run_level),{
  
  t_global <- system.time({
    mifs_global <- foreach(i=1:measles_Nglobal[run_level],.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  mif2(
      mifs_local[[1]],
      start=c(apply(measles_box,1,function(x)runif(1,x[1],x[2])),measles_fixed_params)
    )
  })
},seed=1270401374,kind="L'Ecuyer")

stopCluster(cl)


##--  Evaluate our global search  --##

cl <- makeCluster(12)
registerDoParallel(cl)
clusterExport(cl, c("run_level", "measles_Np", "measles_Neval", "m1", "mifs_global"))

stew(file=sprintf("lik_global_eval-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:measles_Nglobal[run_level],.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(measles_Neval[run_level], logLik(pfilter(m1,params=coef(mifs_global[[i]]),Np=measles_Np[run_level])))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

stopCluster(cl)

results_global <- data.frame(logLik=liks_global[,1],logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))
summary(results_global$logLik,digits=5)
