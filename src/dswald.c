// Shifted Wald (Inverse Gaussian) Distribution

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double d_swald(double t, double alpha, double nu, double theta, int give_log)
{
    double d;

    if(give_log)
      d = log(alpha) +  (-.5) * (log(2) + log(M_PI) + 3*log(t-theta)) 
          + ( -pow((alpha-nu*(t-theta)),2) / (2*(t-theta)) );
    else
      d = alpha * pow(2*M_PI*pow((t-theta),3), -.5) 
          * exp( -pow(alpha-nu*(t-theta),2) / (2*(t-theta)));
  
    return(d);
}

SEXP dswald(SEXP t, SEXP alpha, SEXP nu, SEXP theta, SEXP give_log) {
  double d;
  SEXP value;
    
  d = d_swald(REAL(t)[0], REAL(alpha)[0], REAL(nu)[0], REAL(theta)[0], LOGICAL(give_log)[0]);
    
  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = d;
  UNPROTECT(1);
  return value;
}     
