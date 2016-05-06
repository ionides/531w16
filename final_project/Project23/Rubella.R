Rubella_data = read.csv("Rubella_MI.csv", header=T)
Rubella_data$Time = Rubella_data$Year + (Rubella_data$Month-1)/12
colnames(Rubella_data) = c("year","month","cases","births","pop","time")

library(pomp)
packageVersion("pomp")


Rubella_statenames <- c("SB1","SB2","SB3","SB4","SB5","SB6","IB","SO","V","IO")
Rubella_obsnames <- "cases"
Rubella_t0 <- 1969

Rubella_K <- 7
Rubella_tcovar <- Rubella_data$time
Rubella_bspline_basis <- periodic.bspline.basis(Rubella_tcovar,nbasis=Rubella_K,degree=3,period=1)
colnames(Rubella_bspline_basis)<- paste("xi",1:Rubella_K,sep="")
covartable <- data.frame(
  time=Rubella_tcovar,
  Rubella_bspline_basis,
  B=Rubella_data$births,
  P=predict(smooth.spline(x=1968:1980,y=Rubella_data$pop[12*(1:13)]),
            x=Rubella_tcovar)$y
)


Rubella_rp_names <- c("b1","b2","b3","b4","b5","b6","b7","psi","rho","tau","sigma_dem","sigma_env")

Rubella_ivp_names <- c("SO_0","IO_0")

Rubella_fp_names <- c("delta","K","SB1_0","SB2_0","SB3_0","SB4_0","SB5_0","SB6_0")

Rubella_paramnames <- c(Rubella_rp_names,Rubella_ivp_names,Rubella_fp_names)

Rubella_fixed_params <- c(delta=1/60,K=Rubella_K,SB1_0=Rubella_data$births[12],
                          SB2_0=Rubella_data$births[11],SB3_0=Rubella_data$births[10],
                          SB4_0=Rubella_data$births[9],SB5_0=Rubella_data$births[8],
                          SB6_0=Rubella_data$births[7])


Rubella_rprocess <- Csnippet("
  double lambda, beta, var_epsilon, p, q;
                           
beta = exp(dot_product( (int) K, &xi1, &b1));
lambda = (beta * (IO+IB) / P + psi);
var_epsilon = pow(sigma_dem,2)/ lambda +  pow(sigma_env,2);
lambda *= (var_epsilon < 1.0e-6) ? 1 : rgamma(1/var_epsilon,var_epsilon);
p = exp(- (delta+lambda)/12);
q = (1-p)*lambda/(delta+lambda);
SB1 = B;
SB2= SB1*p;
SB3=SB2*p;
SB4=SB3*p;
SB5=SB4*p;
SB6=SB5*p;
SO= (SB6)*p;
V = (V+SO)*p;
IB=(SB1+SB2+SB3+SB4+SB5+SB6)*q;
IO=(SO+V)*q;
                           ")



Rubella_dmeasure <- Csnippet("
double tol = 1.0e-25;
double mean_cases = rho*IO;
double sd_cases = sqrt(pow(tau*IO,2) + mean_cases);
if(cases > 0.0){
lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) - pnorm(cases-0.5,mean_cases,sd_cases,1,0) + tol; 
} else{
lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) + tol;
}
if (give_log) lik = log(lik);
                           ")


Rubella_rmeasure <- Csnippet("
  cases = rnorm(rho*IO, sqrt( pow(tau*IO,2) + rho*IO ) );
if (cases > 0.0) {
cases = nearbyint(cases);
} else {
cases = 0.0;
}
                           ")


Rubella_initializer <- Csnippet("
  SB1 = SB1_0;
  SB2 = SB2_0;
  SB3 = SB3_0;
  SB4 = SB4_0;
  SB5 = SB5_0;
  SB6 = SB6_0;
  IB = 0;
  V = 0;
  IO = IO_0 * P;
  SO = SO_0 * P;
                              ")


Rubella_toEstimationScale <- Csnippet("
 Tpsi = log(psi);
 Trho = logit(rho);
 Ttau = log(tau);
 Tsigma_dem = log(sigma_dem);
 Tsigma_env = log(sigma_env);
 TSO_0 =  logit(SO_0);
 TIO_0 = logit(IO_0);
                                    ")

Rubella_fromEstimationScale <- Csnippet("
 Tpsi = exp(psi);
 Trho = expit(rho);
 Ttau = exp(tau);
 Tsigma_dem = exp(sigma_dem);
 Tsigma_env = exp(sigma_env);
 TSO_0 =  expit(SO_0);
 TIO_0 = expit(IO_0);
                                      ")



#########
Rubella_params <- data.matrix(read.csv("Rubella_params.csv",row.names=NULL,header=TRUE))
Rubella_mle <- c(Rubella_params[which.max(Rubella_params[,"logLik"]),][Rubella_paramnames])
#######


Rubella <- pomp(
  data=subset(Rubella_data, 
              (time > Rubella_t0 + 0.01) & (time < 1980+11/12+0.01),	
              select=c("cases","time")),
  times="time",
  t0=Rubella_t0,
  params=Rubella_mle,
  rprocess = euler.sim(step.fun = Rubella_rprocess, delta.t=1/12),
  rmeasure= Rubella_rmeasure,
  dmeasure = Rubella_dmeasure,
  covar=covartable,
  tcovar="time",
  obsnames = Rubella_obsnames,
  statenames = Rubella_statenames,
  paramnames = Rubella_paramnames,
  covarnames = c("xi1","B","P"),
  initializer=Rubella_initializer,
  toEstimationScale=Rubella_toEstimationScale, 
  fromEstimationScale=Rubella_fromEstimationScale
)
plot(Rubella)


run_level = 2
Rubella_Np <-          c(100,5e3,1e4)
Rubella_Nmif <-        c(10, 200,400)
Rubella_Nreps_eval <-  c(2,  10,  20)
Rubella_Nreps_local <- c(10, 20, 40)
Rubella_Nreps_global <-c(10, 20, 100)
Rubella_Nsim <-        c(50,100, 500) 


require(doParallel)
registerDoParallel()

stew(file=sprintf("pf1-%d.rda",run_level),{
  t1 <- system.time(
    pf1 <- foreach(i=1:20,.packages='pomp',
                   .options.multicore=list(set.seed=TRUE)) %dopar% try(
                     pfilter(Rubella,Np=Rubella_Np[run_level])
                   )
  )
},seed=493536993,kind="L'Ecuyer")
(L1 <- logmeanexp(sapply(pf1,logLik),se=TRUE))



stew(sprintf("persistence-%d.rda",run_level),{
  t_sim <- system.time(
    sim <- foreach(i=1:Rubella_Nsim[run_level],.packages='pomp',
                   .options.multicore=list(set.seed=TRUE)) %dopar% 
      simulate(Rubella)
  )
},seed=493536993,kind="L'Ecuyer")

no_cases_data <- sum(obs(Rubella)==0)
no_cases_sim <- sum(sapply(sim,obs)==0)/length(sim)
fadeout1_sim <- sum(sapply(sim,function(Ru)states(Ru)["IB",]+states(Ru)["IO",]<1))/length(sim)
fadeout100_sim <- sum(sapply(sim,function(Ru)states(Ru)["IB",]+states(Ru)["IO",]<100))/length(sim)
imports_sim <- coef(Rubella)["psi"]*mean(sapply(sim,function(Ru) mean(states(Ru)["V",]+states(Ru)["SO",]+states(Ru)["SB1",]+states(Ru)["SB2",]+states(Ru)["SB3",]+states(Ru)["SB4",]+states(Ru)["SB5",]+states(Ru)["SB6",])))/12

mle_simulation <- simulate(Rubella,seed=127)
plot(mle_simulation)



Rubella_rw.sd_rp <- 0.02
Rubella_rw.sd_ivp <- 0.2
Rubella_cooling.fraction.50 <- 0.5

stew(sprintf("mif-%d.rda",run_level),{
  t2 <- system.time({
    m2 <- foreach(i=1:Rubella_Nreps_local[run_level],
                  .packages='pomp', .combine=c,
                  .options.multicore=list(set.seed=TRUE)) %dopar% try(
                    mif2(Rubella,
                         Np=Rubella_Np[run_level],
                         Nmif=Rubella_Nmif[run_level],
                         cooling.type="geometric",
                         cooling.fraction.50=Rubella_cooling.fraction.50,
                         transform=TRUE,
                         rw.sd=rw.sd(
                           b1=Rubella_rw.sd_rp,
                           b2=Rubella_rw.sd_rp,
                           b3=Rubella_rw.sd_rp,
                           b4=Rubella_rw.sd_rp,
                           b5=Rubella_rw.sd_rp,
                           b6=Rubella_rw.sd_rp,
                           b7=Rubella_rw.sd_rp,
                           psi=Rubella_rw.sd_rp,
                           rho=Rubella_rw.sd_rp,
                           tau=Rubella_rw.sd_rp,
                           sigma_dem=Rubella_rw.sd_rp,
                           sigma_env=Rubella_rw.sd_rp,
                           IO_0=ivp(Rubella_rw.sd_ivp),
                           SO_0=ivp(Rubella_rw.sd_ivp)
                         )
                    )
                  )
    
    lik_m2 <- foreach(i=1:Rubella_Nreps_local[run_level],.packages='pomp',
                      .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% 
                      {
                        logmeanexp(
                          replicate(Rubella_Nreps_eval[run_level],
                                    logLik(pfilter(Rubella,params=coef(m2[[i]]),Np=Rubella_Np[run_level]))
                          ),
                          se=TRUE)
                      }
  })
},seed=318817883,kind="L'Ecuyer")

r2 <- data.frame(logLik=lik_m2[,1],logLik_se=lik_m2[,2],t(sapply(m2,coef)))
summary(r2$logLik, digits=5)


pairs(~logLik+psi+rho+tau+sigma_dem+sigma_env,data=subset(r2,logLik>max(logLik)-2000))



######

Rubella_box <- rbind(
  b1=c(-2,8),
  b2=c(-2,8),
  b3=c(-2,8),
  b4=c(-2,8),
  b5=c(-2,8),
  b6=c(-2,8),
  b7=c(-2,8),
  psi=c(0,0.1),
  rho=c(0,0.1),
  tau=c(0,0.1),
  sigma_dem=c(0,0.5),
  sigma_env=c(0,1),
  SO_0=c(0,1),
  IO_0=c(0,0.01)
)


stew(file=sprintf("box_eval-%d.rda",run_level),{
  t3 <- system.time({
    m3 <- foreach(i=1:Rubella_Nreps_global[run_level],.packages='pomp',.combine=c,
                  .options.multicore=list(set.seed=TRUE)) %dopar%  
      mif2(
        m2[[1]],
        start=c(apply(Rubella_box,1,function(x)runif(1,x[1],x[2])),Rubella_fixed_params)
      )
    
    lik_m3 <- foreach(i=1:Rubella_Nreps_global[run_level],.packages='pomp',.combine=rbind,
                      .options.multicore=list(set.seed=TRUE)) %dopar% {
                        set.seed(87932+i)
                        logmeanexp(
                          replicate(Rubella_Nreps_eval[run_level],
                                    logLik(pfilter(Rubella,params=coef(m3[[i]]),Np=Rubella_Np[run_level]))
                          ), 
                          se=TRUE)
                      }
  })
},seed=290860873,kind="L'Ecuyer")


r3 <- data.frame(logLik=lik_m3[,1],logLik_se=lik_m3[,2],t(sapply(m3,coef)))
if(run_level>1) write.csv(rbind(r2,r3),file="Rubella_params.csv",row.names=FALSE)
summary(r3$logLik,digits=5)

pairs(~logLik+psi+rho+tau+sigma_dem+sigma_env,data=subset(r3,logLik>max(logLik)-200))

class(m3)
class(m3[[1]])
plot(m3[r3$logLik>max(r3$logLik)-1000])

loglik_convergence <- do.call(cbind,conv.rec(m3[r3$logLik>max(r3$logLik)-10],"loglik"))
matplot(loglik_convergence,type="l",lty=1,ylim=max(loglik_convergence,na.rm=T)+c(-10,0))
