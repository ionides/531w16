sir_step <- Csnippet("
                     double dN_SL = rbinom(S,1-exp(-(1-k)*L*dt));
                     double dN_SA = rbinom(S,1-exp(-(k)*A*dt));
                     double dN_LA = rbinom(L,1-exp(-(x2)*A*dt));
                     double dN_RA = rbinom(R,1-exp(-(x7+c*lambda)*A*dt));
                     double dN_AO = rbinom(A,1-exp(-(x3-p*x3)*O*dt));
                     double dN_AB = rbinom(A,1-exp(-(p*x3)*B*dt));
                     double dN_EB = rbinom(E,1-exp(-(x5)*B*dt));
                     double dN_BE = rbinom(B,1-exp(-(1-phi)*E*dt));
                     double dN_BC = rbinom(B,1-exp(-(phi)*C*dt));
                     double dN_CR = rbinom(C,1-exp(-(alpha)*R*dt));
                     double dN_CE = rbinom(C,1-exp(-(1-alpha)*E*dt));
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

double dN_SL = rbinom(S,1-exp(-(1-k)*L*dt));
double dN_SA = 
  double dN_LA = 
  double dN_RA = 
  double dN_AO = 
  double dN_AB = 
  double dN_EB = 
  double dN_BE = 
  double dN_BC = 
  double dN_CR = 
  double dN_CE = 