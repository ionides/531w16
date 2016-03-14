# rm(list=ls())
require(pomp)
# setwd('/Users/Micaela/Documents/Biological_Rhythms_manuscript/Stevenson_2015/chickenpox_model')

#######################################################################################################
# Load Covars & Data
#######################################################################################################
# load the chickenpox and shingles data for Australia
dat<- read.csv('CPox_Shing_AUS06_15.csv',header=TRUE)
dat$time<- dat$YEAR + dat$MONTH/12
names(dat)<- c('year','month','chickenpox','shingles','time')
dat<- subset(dat,chickenpox>=0)

# load the google chickenpox data
goog<- read.csv('google_trends_data.csv',header=TRUE)
goog$mo<- NA
for(i in 1:dim(goog)[1]){
  month<- goog$Month[i] 
  if(month=='Jan'){goog$mo[i]<- 1}
  if(month=='Feb'){goog$mo[i]<- 2}
  if(month=='Mar'){goog$mo[i]<- 3}
  if(month=='Apr'){goog$mo[i]<- 4}
  if(month=='May'){goog$mo[i]<- 5}
  if(month=='Jun'){goog$mo[i]<- 6}
  if(month=='Jul'){goog$mo[i]<- 7}
  if(month=='Aug'){goog$mo[i]<- 8}
  if(month=='Sep'){goog$mo[i]<- 9}
  if(month=='Oct'){goog$mo[i]<- 10}
  if(month=='Nov'){goog$mo[i]<- 11}
  if(month=='Dec'){goog$mo[i]<- 12}
}
goog$time<- goog$Year + goog$mo/12
goog<- goog[order(goog$time,decreasing=FALSE),]
goog$time<- round(goog$time,digits=3)

# load the population data
pop<- read.csv('Australia_pop.csv',header=TRUE)
pop2010<- sum(subset(pop,year==2010)$pop) # data from Australian National University (http://adsri.anu.edu.au/demo-stats/aust)

# vaccine
vacc<- data.frame(year=2006:2010)
vacc$vaccine[which(vacc$year==2006)]<- 0.3  # account for vaccination from Heywood et al 2014 (http://www.scielosp.org/pdf/bwho/v92n8/0042-9686-bwho-92-08-593.pdf)
vacc$vaccine[which(vacc$year==2007)]<- 0.82  # account for vaccination
vacc$vaccine[which(vacc$year==2008)]<- 0.88  # account for vaccination
vacc$vaccine[which(vacc$year==2009)]<- 0.90  # account for vaccination
vacc$vaccine[which(vacc$year==2010)]<- 0.91
#######################################################################################################
# Functions we need to transform params
#######################################################################################################
expit=function(x){1.0/(1.0+exp(-x))}
logit= function(p) {log(p/(1-p))};

#######################################################################################################
# Set up a data frame containing the covariates (don't include the data)
#######################################################################################################
covar.table<- data.frame(timestep=seq.int(from=0,by=1,length=dim(dat)[1]+1))
covar.table$time<- c(2006,round(dat$time,digits=3))
covar.table$children <- pop2010 
# account for vaccination, still need data form 2010--2015
covar.table$vaccine<- predict(smooth.spline(vacc$year,vacc$vaccine),covar.table$time)$y

# plot(covar.table$time,covar.table$vaccine,type='l')
# points(vacc$year,vacc$vaccine)

for(i in 1:dim(covar.table)[1]){
  TIME<- covar.table$time[i]
  if(TIME<=2015){
    covar.table$googleLAG0[i]<- goog$Chicken.Pox.Australia[which(goog$time==TIME)]
    covar.table$googleLAG1[i]<- goog$Chicken.Pox.Australia[which(goog$time==TIME)-1] 
    covar.table$googleLAG2[i]<- goog$Chicken.Pox.Australia[which(goog$time==TIME)-2]
  }
  else{covar.table$googleLAG0[i]<- NA
    covar.table$googleLAG1[i]<- NA
    covar.table$googleLAG2[i]<- NA
  } 
}
covar.table[,c('googleLAG0','googleLAG1','googleLAG2')]
covar.table<- subset(covar.table,time<= 2015)

#######################################################################################################
# Set up data frame with the data to be fit, this dataframe should start at timestep==1
#######################################################################################################
monthly.cases<-data.frame(cases=subset(dat,time<=2015)$chickenpox,time=round(subset(dat,time<=2015)$time,digits=3))
monthly.cases$timestep<- 1:108

# # the google trends and the data align well
# par(plt=c(0.1,0.9,0.6,0.9))
# plot(monthly.cases$time,monthly.cases$cases,type='l')
# par(plt=c(0.1,0.9,0.1,0.5),new=TRUE)
# plot(covar.table$time,covar.table$googleLAG0,type='l',col='red')

#######################################################################################################
# Compile/load model from C-code
#######################################################################################################
if (is.loaded("chickenpox_meas_dens")) dyn.unload("chickenpox_seas_E.so")
system("R CMD SHLIB chickenpox_seas_E.c")
dyn.load("chickenpox_seas_E.so")

#######################################################################################################
# POMP Object
#######################################################################################################
my.model<- pomp(
              data=monthly.cases[c("timestep","cases")], # data are cases and associated timestep, starting w/ timestep==1
              times="timestep", # name of time variable in data dataframe
              rprocess=euler.sim("chickenpox_proc_sim",delta.t=1,PACKAGE="chickenpox_seas_E"), # chickenpox_proc_sim is the name of our process model in C-code; chickenpox.c is the name of our C file; see ?plugins for integration options
              rmeasure="chickenpox_meas_sim", # name of measurement model in C-code
              dmeasure="chickenpox_meas_dens", # name of measurement density model in C-code
              PACKAGE="chickenpox_seas_E", # name of C-file
              t0=0, # this is the time we initialize the model
              obsnames="cases", # name of observations in data dataframe
              statenames=c("FOI","I"), # process model state variable names, this vector must match the order of the states in (x[stateindex[i]]) in C-code
              covarnames=c("children","vaccine","googleLAG0","googleLAG1","googleLAG2"), # covariate names, must match order of (covar[covindex[i]]) in C-code
              paramnames=c('alpha','tau','beta1','beta2','beta3','beta.sd','omega','rho'), # param names, must match order of (p[parindex[0]]) in C-code, but notice the actual name doesn't have to match for instance here we use "alpha" in C-code we call it LOGALPHA, only the indexing matters not the exact name
              tcovar="timestep", # name of the time variable in the covariate data frame
              covar=covar.table, # name of the covariate dataframe
	     	      initializer=function(params,t0, ...){ # the initializer allows you to specify how to set up initial conditions using your parameter set "params"
    	      	  p <- expit(params)
  	      	    covars<- as.data.frame(covar.table)
  			        with(as.list(p),{
             	    	x0=c(FOI.0,I.0)
  					        x0[1]= 0;  
                    x0[2]= 0; 
                 		names(x0)=c("FOI","I");
                         x0
                    })
                }
            )
              
#######################################################################################################
#END standard code
#######################################################################################################
#######################################################################################################

#######################################################################################################
# test the model with a random parameter set
#######################################################################################################

# param.set<- read.csv('MIF_results_phase6_run1.csv',header=TRUE)
# param.set$omega<- log(2)
# param.set$beta3<- log(0.000001)
# param.set<- subset(param.set,LogLik==max(param.set$LogLik))
# param.set<- unlist(param.set)
# for(i in 1:3){
#   sim<- simulate(my.model,param=param.set)
#   dev.new()
#   plot(sim)
# }

# plot(as.vector(my.model@data),type='l',ylim=c(0,500))
# lines(as.vector(sim@data),col='red')

