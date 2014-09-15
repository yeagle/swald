#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double LaguerreL(double n, double a, double x) 
{
  double L;
  L = ((n+a+1)*beta(1+a,n+1));
  return L;
}

double dwald_d(double t, double lambda, double alpha, double tau, double kappa)
{
  double d;

  d = -sqrt(0.2e1) * pow(0.2e1, tau / 0.2e1 + 0.1e1 / 0.2e1) * alpha * exp(-0.1e1 / t * alpha * alpha * lambda / 0.2e1) * pow(kappa, -tau - 0.3e1) * pow(lambda, -tau / 0.2e1 - 0.1e1) * pow(t, -tau / 0.2e1 - 0.3e1) * 0.3141592654e1 * (-sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.5e1 / 0.2e1) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.2e1) * alpha * alpha * pow(kappa, 0.3e1) + sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.5e1 / 0.2e1) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.2e1) * alpha * alpha * pow(kappa, 0.3e1) - sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.3e1 / 0.2e1) * pow(t, 0.3e1 / 0.2e1) * cos(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.2e1) * pow(kappa, 0.3e1) + 0.2e1 * sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.3e1 / 0.2e1) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.2e1) * alpha * kappa * kappa - 0.2e1 * sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.3e1 / 0.2e1) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.2e1) * alpha * kappa * kappa - sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * sqrt(lambda) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.2e1) * kappa + sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * sqrt(lambda) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.2e1) * kappa + sin(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(alpha, 0.3e1) * pow(kappa, 0.3e1) * pow(lambda, 0.3e1) - sin(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(alpha, 0.3e1) * pow(kappa, 0.3e1) * pow(lambda, 0.3e1) + sin(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * pow(kappa, 0.3e1) * lambda * lambda * tau * t - sin(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * pow(kappa, 0.3e1) * lambda * lambda * t - 0.3e1 * sin(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * alpha * kappa * kappa * lambda * lambda + 0.3e1 * sin(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * alpha * kappa * kappa * lambda * lambda - sin(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * kappa * kappa * lambda * tau * t + sin(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * kappa * kappa * lambda * t + 0.3e1 * sin(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * kappa * lambda - 0.3e1 * sin(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * kappa * lambda - sin(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) + sin(0.3141592654e1 * tau / 0.2e1) * gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1)) / gamma(tau) / sin(0.3141592654e1 * tau / 0.2e1) / gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) / cos(0.3141592654e1 * tau / 0.2e1) / gamma(-tau / 0.2e1 + 0.2e1) / 0.16e2;

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
