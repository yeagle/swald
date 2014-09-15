#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double dwald_d(double q, double a, double v, double s)
{
  double d;

  return d;
}

SEXP dwald(SEXP q, SEXP a, SEXP v, SEXP s) {
  double d;
  SEXP value;

  d = dwald_d(REAL(q)[0], REAL(a)[0], REAL(v)[0], REAL(s)[0]);

  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = d;
  UNPROTECT(1);
  return value;
}
