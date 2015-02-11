#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double dIG_d(double t, double alpha, double lambda, double nu)
{
    double d;

    d = alpha * pow(lambda / (2*M_PI*pow(t,3)),.5) * exp( -lambda * pow(nu*t-alpha,2) / (2*t));
  
    return(d);
}

SEXP  dIG(SEXP t, SEXP alpha, SEXP lambda, SEXP nu) {
  double d;
  SEXP value;
    
  d = dIG_d(REAL(t)[0], REAL(alpha)[0], REAL(lambda)[0], REAL(nu)[0]);
    
  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = d;
  UNPROTECT(1);
  return value;
}     