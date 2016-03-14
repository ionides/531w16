// dear emacs, please treat this as -*- C++ -*-
// This code was last edited by Micaela Martinez-Bakker on June 19, 2015
#include <R.h>
#include <Rmath.h>
#include <math.h>

///////////////////////////////////////////////////////////////////
// SET UP SLOTS TO HOLD PARAMETERS, STATE VARIABLES, AND COVARIATES
///////////////////////////////////////////////////////////////////
// define parameter sets (if param is bounded between 0-1 use LOGIT otherwise LOG for params bounded >=0)
#define LOGTAU         (p[parindex[0]]) // observation dispersion parameter
#define LOGBETA1       (p[parindex[1]]) // scaler for cosine
#define LOGBETA3       (p[parindex[2]]) // baseline force of infection
#define LOGBETA_SD     (p[parindex[3]]) // standard deviation of gamma distributed process noise
#define LOGITRHO       (p[parindex[4]]) // mean report rate 
#define LOGOMEGA       (p[parindex[5]]) // the phase of the cosine function for seasonal forcing

// define process model state variables
#define FOI         (x[stateindex[0]]) // monthly force of infection, monthly per capita rate of infection for the 0-14 year old age group
#define I           (x[stateindex[1]]) // infected individuals

// define covariates
#define MONTH_YR    (covar[covindex[0]]) // time value in format "year + month/12", where Jan=1 and Dec=12
#define CHILDREN    (covar[covindex[1]]) // population size 0-14 yrs of age

// observation model state variables
#define CHICKENPOX       (y[obsindex[0]]) // observed CHICKENPOX cases

///////////////////////////////////////////////////////////////////
// DEFINE FUNCTIONS NEEDED FOR TRANSFORMING PARAMETERS
///////////////////////////////////////////////////////////////////
// define LOGIT and EXPIT functions

// static double logit (double p) {
//   return log(p/(1-p));
// }

static double expit (double x) {
  return 1.0/(1.0+exp(-x));
}

///////////////////////////////////////////////////////////////////
// DEFINE THE MEASUREMENT DENSITY FOR CALCULATING THE LIKELIHOOD
///////////////////////////////////////////////////////////////////
void null_chickenpox_meas_dens (double *lik, double *y, double *x, double *p, int give_log,
			int *obsindex, int *stateindex, int *parindex, int *covindex,
			int ncovar, double *covar, double t) {
  
  double report_rate, tau;
  double tol = 1.0e-18;
  // PUT PARS ON NATURAL SCALE
  tau = exp(LOGTAU);
  report_rate = expit(LOGITRHO);

  if(CHICKENPOX > 0.0){
   *lik = pnorm(CHICKENPOX+0.5,report_rate*I,tau*I,1,0) - pnorm(CHICKENPOX-0.5,report_rate*I,tau*I,1,0)+ tol;  
  } else{
   *lik = pnorm(CHICKENPOX+0.5,report_rate*I,tau*I,1,0)+ tol;  
  }
  if (give_log) *lik = log(*lik);
  if (!isfinite(*lik)) Rprintf("chickenpox_meas_dens %lg %lg %lg %lg %lg\n",CHICKENPOX,report_rate,tau,I,*lik);
}

///////////////////////////////////////////////////////////////////
// DEFINE THE MEASUREMENT MODEL
///////////////////////////////////////////////////////////////////
void null_chickenpox_meas_sim (double *y, double *x, double *p, 
		       int *obsindex, int *stateindex, int *parindex, int *covindex,
		       int ncovar, double *covar, double t) {
  
  double report_rate, tau;
  // PUT PARS ON NATURAL SCALE
  tau = exp(LOGTAU);
  report_rate = expit(LOGITRHO);

  
  CHICKENPOX = rnorm(report_rate*I,tau*I);
  if (CHICKENPOX >= 0) {
    CHICKENPOX = nearbyint(CHICKENPOX);
  	} else {CHICKENPOX = 0;}

}

///////////////////////////////////////////////////////////////////
// DEFINE THE PROCESS MODEL, THIS IS A STATISTICAL NON-MECHANISTIC MODEL
///////////////////////////////////////////////////////////////////
void null_chickenpox_proc_sim (double *x, const double *p, 
		     const int *stateindex, const int *parindex, const int *covindex,
		     int covdim, const double *covar, 
		     double t, double dt)
{
 double beta_sd;
 double epsilon;
 double beta1;
 double beta3;
 double scale;
 double omega;

 // PUT PARS ON NATURAL SCALE
 beta_sd = exp(LOGBETA_SD);
 beta1 = exp(LOGBETA1);
 beta3 = exp(LOGBETA3);
 omega = exp(LOGOMEGA);


  if (beta_sd > 0) {
    scale = pow(beta_sd,2);
    epsilon = rgamma(1/scale,scale); 
  	} else {epsilon = 1;}

  if(
     isfinite(epsilon)== FALSE ||
     isfinite(FOI)==FALSE ||
     isfinite(I)==FALSE ||
     isfinite(beta1)==FALSE
     )
    {
      Rprintf("non finite value in chickenpox_proc_sim\n");
      return;
    }

  FOI = (beta1*cos((2*M_PI/12)*(t+omega)) + beta3)*epsilon;
  if(FOI < 0){FOI = 0;}
  I = CHILDREN*FOI;

}