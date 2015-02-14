// Driftrate as truncated Normal - Wald Distribution

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double dwald_trunc_d(double t, double lambda, double alpha, double v, double d, int give_log)
{
  double d;

  if(give_log)
  {
    d = log(alpha) + 
        0.5 * ( log(lambda) - log(2) - log(M_PI) - 3 * log(t) - log(lambda * t * v + 1)) -
        log(pnorm(d / sqrt(v), 0, 1, 1, 0) ) -
        (lambda * pow(d * t - alpha,2)) / (2 * t * (lambda * t * v + 1)) +
        log(pnorm((lambda * alpha * v + d) / (sqrt(lambda * t * pow(v,2) + v)), 0, 1, 1, 0 ));
  }
  else
  {
    d = alpha * sqrt( lambda / (2 * M_PI * pow(t, 3) * (lambda * t * v + 1)) ) *
        1 / pnorm(d / sqrt(v), 0, 1, 1, 0) *
        exp( - (lambda * pow(d * t - alpha, 2)) / (2 * t * (lambda * t * v + 1)) ) *
        pnorm( (lambda * alpha * v + d) / (sqrt(lambda * t * pow(v, 2) + v) ), 0, 1, 1, 0);
  }

  return d;
}

SEXP  dwald_trunc(SEXP t, SEXP lambda, SEXP alpha, SEXP v, SEXP delta SEXP give_log) { 
  double d;
  SEXP value;
    
  d = dwald_trunc_d(REAL(t)[0], REAL(lambda)[0], REAL(alpha)[0], REAL(v)[0], REAL(delta)[0], LOGICAL(give_log)[0]);
    
  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = d;
  UNPROTECT(1);
  return value;
}    
