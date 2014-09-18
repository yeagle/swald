#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <gsl/gsl_specfunc.h>

double r_gamma(double x)
{
  return exp(gamma(x));
}

double LaguerreL(double n, double a, double x) 
{
  double Lf;

  Lf = 1/((n+a+1)*beta(1+a,n+1)) * gsl_sf_hyperg_1F1(-n,a+1,x); 

  return Lf;
}

double erf(double x)
{
  double e;
  
  e = 2 * pnorm(x * sqrt(2), 0,1, 1,0) - 1;
  
  return e;
}

double dwald_gamma_d(double t, double lambda, double alpha, double tau, double kappa)
{
  double d;

  if (lambda == 1 && alpha == 1 && tau == 1 && kappa == 1 )
  {
    d = exp(- 1 / (2*t)) / (2 * pow(t,2));
  }
  else if (lambda == 1 && tau == 1) 
  {
    d = alpha*exp(-(2*alpha*kappa-1)/(2*t*pow(kappa,2)))* (1 + erf( (alpha*kappa-1) / (kappa*sqrt(2*t)) )) / (2*pow(t,2)*kappa);
  }
  else 
  {  
    d = -sqrt(0.2e1) * pow(0.2e1, tau / 0.2e1 + 0.1e1 / 0.2e1) * alpha * exp(-0.1e1 / t * alpha * alpha * lambda / 0.2e1) * pow(kappa, -tau - 0.3e1) * pow(lambda, -tau / 0.2e1 - 0.1e1) * pow(t, -tau / 0.2e1 - 0.3e1) * 0.3141592654e1 * (-sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.5e1 / 0.2e1) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * alpha * alpha * pow(kappa, 0.3e1) + sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.5e1 / 0.2e1) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * alpha * alpha * pow(kappa, 0.3e1) - sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.3e1 / 0.2e1) * pow(t, 0.3e1 / 0.2e1) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * pow(kappa, 0.3e1) + 0.2e1 * sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.3e1 / 0.2e1) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * alpha * kappa * kappa - 0.2e1 * sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.3e1 / 0.2e1) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * alpha * kappa * kappa - sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * sqrt(lambda) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * kappa + sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * sqrt(lambda) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * kappa + sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(alpha, 0.3e1) * pow(kappa, 0.3e1) * pow(lambda, 0.3e1) - sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(alpha, 0.3e1) * pow(kappa, 0.3e1) * pow(lambda, 0.3e1) + sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * pow(kappa, 0.3e1) * lambda * lambda * tau * t - sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * pow(kappa, 0.3e1) * lambda * lambda * t - 0.3e1 * sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * alpha * kappa * kappa * lambda * lambda + 0.3e1 * sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * alpha * kappa * kappa * lambda * lambda - sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * kappa * kappa * lambda * tau * t + sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * kappa * kappa * lambda * t + 0.3e1 * sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * kappa * lambda - 0.3e1 * sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * kappa * lambda - sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) + sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1)) / r_gamma(tau) / sin(0.3141592654e1 * tau / 0.2e1) / r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) / cos(0.3141592654e1 * tau / 0.2e1) / r_gamma(-tau / 0.2e1 + 0.2e1) / 0.16e2;
  }   
  return d;
}    
     
SEXP  dwald_gamma(SEXP t, SEXP lambda, SEXP alpha, SEXP tau, SEXP kappa) {
  double d;
  SEXP value;
    
  d = dwald_gamma_d(REAL(t)[0], REAL(lambda)[0], REAL(alpha)[0], REAL(tau)[0], REAL(kappa)[0]);
    
  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = d;
  UNPROTECT(1);
  return value;
}    
