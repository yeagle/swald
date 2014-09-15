#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double LaguerreL(int n, double a, double x) 
{
  double Lf;

  if (n == 0) Lf = 1;
  else if (n == 1) Lf = 1+a-x;
  else
  {
    Lf = (2*n-1+a-x)/n * LaguerreL(n-1,a,x) - LaguerreL(n-2,a,x);
  }

  Lf = 1;
  return Lf;
}

double dwald_d(double t, double lambda, double alpha, double tau, double kappa)
{
  double d;

  if (kappa == 1)
  {
    d = exp(- 1 / (2*t)) / (2 * pow(t,2));
  }
  else 
  {
    d = -sqrt(2) * pow(2, tau / 2 + 1 / 2) * alpha *
        exp(-1 / t * alpha * alpha * lambda / 2) * pow(kappa, -tau - 3) * pow(lambda, -tau / 2 - 1) * pow(t, -tau / 2 -
        3) * M_PI * (-sqrt(2) * LaguerreL(-tau / 2, 1 / 2, 1 / lambda * pow(alpha * kappa * lambda -
        1, 2) * pow(kappa, -2) / t / 2) * pow(lambda, 5 / 2) * sqrt(t) * cos(M_PI * tau / 2)
        * gamma(-tau / 2 + 2) * alpha * alpha * pow(kappa, 3) + sqrt(2) * LaguerreL(-tau / 2, 3 / 2, 1 /
        lambda * pow(alpha * kappa * lambda - 1, 2) * pow(kappa, -2) / t / 2) * pow(lambda, 5 / 2) *
        sqrt(t) * cos(M_PI * tau / 2) * gamma(-tau / 2 + 2) * alpha * alpha * pow(kappa, 3) - sqrt(2) *
        LaguerreL(-tau / 2, 1 / 2, 1 / lambda * pow(alpha * kappa * lambda - 1, 2) * pow(kappa, -2) / t /
        2) * pow(lambda, 3 / 2) * pow(t, 3 / 2) * cos(M_PI * tau / 2) * gamma(-tau / 2 + 2) *
        pow(kappa, 3) + 2 * sqrt(2) * LaguerreL(-tau / 2, 1 / 2, 1 / lambda * pow(alpha * kappa * lambda -
        1, 2) * pow(kappa, -2) / t / 2) * pow(lambda, 3 / 2) * sqrt(t) * cos(M_PI * tau / 2)
        * gamma(-tau / 2 + 2) * alpha * kappa * kappa - 2 * sqrt(2) * LaguerreL(-tau / 2, 3 / 2, 1 /
        lambda * pow(alpha * kappa * lambda - 1, 2) * pow(kappa, -2) / t / 2) * pow(lambda, 3 / 2) *
        sqrt(t) * cos(M_PI * tau / 2) * gamma(-tau / 2 + 2) * alpha * kappa * kappa - sqrt(2) *
        LaguerreL(-tau / 2, 1 / 2, 1 / lambda * pow(alpha * kappa * lambda - 1, 2) * pow(kappa, -2) / t /
        2) * sqrt(lambda) * sqrt(t) * cos(M_PI * tau / 2) * gamma(-tau / 2 + 2) * kappa + sqrt(2)
        * LaguerreL(-tau / 2, 3 / 2, 1 / lambda * pow(alpha * kappa * lambda - 1, 2) * pow(kappa,
        -2) / t / 2) * sqrt(lambda) * sqrt(t) * cos(M_PI * tau / 2) * gamma(-tau / 2 +
        2) * kappa + sin(M_PI * tau / 2) * gamma(-tau / 2 + 3 / 2) * LaguerreL(-tau / 2
        + 1 / 2, 1 / 2, 1 / lambda * pow(alpha * kappa * lambda - 1, 2) * pow(kappa,
        -2) / t / 2) * pow(alpha, 3) * pow(kappa, 3) * pow(lambda, 3) - sin(M_PI * tau
        / 2) * gamma(-tau / 2 + 3 / 2) * LaguerreL(-tau / 2 + 1 / 2, 3 /
        2, 1 / lambda * pow(alpha * kappa * lambda - 1, 2) * pow(kappa, -2) / t
        / 2) * pow(alpha, 3) * pow(kappa, 3) * pow(lambda, 3) + sin(M_PI * tau /
        2) * gamma(-tau / 2 + 3 / 2) * LaguerreL(-tau / 2 + 1 / 2, 1 /
        2, 1 / lambda * pow(alpha * kappa * lambda - 1, 2) * pow(kappa, -2) / t
        / 2) * alpha * pow(kappa, 3) * lambda * lambda * tau * t - sin(M_PI * tau / 2)
        * gamma(-tau / 2 + 3 / 2) * LaguerreL(-tau / 2 + 1 / 2, 1 /
        2, 1 / lambda * pow(alpha * kappa * lambda - 1, 2) * pow(kappa, -2) / t
        / 2) * alpha * pow(kappa, 3) * lambda * lambda * t - 3 * sin(M_PI * tau /
        2) * gamma(-tau / 2 + 3 / 2) * LaguerreL(-tau / 2 + 1 / 2, 1 / 2, 1 / lambda * pow(alpha * kappa *
          lambda - 1, 2) * pow(kappa, -2) / t / 2) * alpha * alpha
        * kappa * kappa * lambda * lambda + 3 * sin(M_PI * tau / 2)
        * gamma(-tau / 2 + 3 / 2) * LaguerreL(-tau / 2 + 1 / 2, 3 / 2, 1 / lambda * pow(alpha * kappa *
          lambda - 1, 2) * pow(kappa, -2) / t / 2) * alpha * alpha * kappa * kappa * lambda * lambda - sin(M_PI * tau / 2) *
        gamma(-tau / 2 + 3 / 2) * LaguerreL(-tau / 2 + 1 / 2, 1 / 2, 1 / lambda * pow(alpha * kappa * lambda - 1,
          2) * pow(kappa, -2) / t / 2) * kappa * kappa * lambda * tau * t + sin(M_PI * tau / 2) * gamma(-tau / 2 + 3 /
        2) * LaguerreL(-tau / 2 + 1 / 2, 1 / 2, 1 / lambda * pow(alpha * kappa * lambda - 1, 2) *
        pow(kappa, -2) / t / 2) * kappa * kappa * lambda * t + 3 * sin(M_PI * tau / 2) * gamma(-tau / 2 + 3
        / 2) * LaguerreL(-tau / 2 + 1 / 2, 1 / 2, 1 / lambda * pow(alpha * kappa * lambda - 1,
        2) * pow(kappa, -2) / t / 2) * alpha * kappa * lambda - 3 * sin(M_PI * tau / 2) *
        gamma(-tau / 2 + 3 / 2) * LaguerreL(-tau / 2 + 1 / 2, 3 / 2, 1 / lambda *
        pow(alpha * kappa * lambda - 1, 2) * pow(kappa, -2) / t / 2) * alpha * kappa * lambda -
        sin(M_PI * tau / 2) * gamma(-tau / 2 + 3 / 2) * LaguerreL(-tau / 2 + 1 / 2,
        1 / 2, 1 / lambda * pow(alpha * kappa * lambda - 1, 2) * pow(kappa, -2) / t /
        2) + sin(M_PI * tau / 2) * gamma(-tau / 2 + 3 / 2) * LaguerreL(-tau / 2
        + 1 / 2, 3 / 2, 1 / lambda * pow(alpha * kappa * lambda - 1,
        2) * pow(kappa, -2) / t / 2)) / gamma(tau) / sin(M_PI * tau /
        2) / gamma(-tau / 2 + 3 / 2) /
        cos(M_PI * tau / 2) / gamma(-tau / 2 + 2) / 16;
  }

  return d;
}

SEXP dwald(SEXP t, SEXP lambda, SEXP alpha, SEXP tau, SEXP kappa) {
  double d;
  SEXP value;

  d = dwald_d(REAL(t)[0], REAL(lambda)[0], REAL(alpha)[0], REAL(tau)[0], REAL(kappa)[0]);

  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = d;
  UNPROTECT(1);
  return value;
}
