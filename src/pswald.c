#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double p_swald(double t, double alpha, double nu, double theta, int lower_tail, int log_p)
{
    double p;
    double x;

    if(log_p)
      x = exp(t);
    else
      x = t;

    p = pnorm((nu*(x-theta)-alpha) / sqrt((x-theta)), 0,1,1,0) + 
        exp(2*alpha*nu) * pnorm(-(nu*(x-theta)+alpha) / sqrt((x-theta)), 0,1,1,0);
    
    return (lower_tail ? p : 1-p);
}

SEXP pswald(SEXP t, SEXP alpha, SEXP nu, SEXP theta, SEXP lower_tail, SEXP log_p) {
  double p;
  SEXP value;
    
  p = p_swald(REAL(t)[0], REAL(alpha)[0], REAL(nu)[0], REAL(theta)[0], LOGICAL(lower_tail)[0], LOGICAL(log_p)[0]);
    
  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = p;
  UNPROTECT(1);
  return value;
}     
