#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double dShiftedWald_d(double t, double alpha, double gamma, double theta)
{
    double d;

    d = alpha * pow(2*M_PI*((t-theta)^3), -.5) * exp( -pow(alpha-gamma*(t-theta),2) / (2*(t-theta)));
  
    return(d);
}

SEXP  dShiftedWald(SEXP t, SEXP alpha, SEXP gamma, SEXP theta) {
  double d;
  SEXP value;
    
  d = dShiftedWald_d(REAL(t)[0], REAL(alpha)[0], REAL(gamma)[0], REAL(theta)[0]);
    
  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = d;
  UNPROTECT(1);
  return value;
}     