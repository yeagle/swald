#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double LaguerreL(int n, double a, double x) 
{
  double Lf;

  if (n <= 0) Lf = 1;
  else if (n <= 1) Lf = 1+a-x;
  else
  {
    Lf = (2*n-1+a-x)/n * LaguerreL(n-1,a,x) - LaguerreL(n-2,a,x);
  }

  return Lf;
}


double erf(double x)
{
  double e;
  
  e = 2 * pnorm(x * sqrt(2)) - 1;
  
  return e;
}


double dwald_d(double t, double lambda, double alpha, double tau, double kappa)
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
    d = 0; 
    /* 
      -(1/16)*sqrt(2)*pow(2,((1/2)*tau+1/2))*alpha*exp(-(1/2)*pow(alpha,2)*lambda/t)*
    pow(kappa,(-tau-3))*pow(lambda,(-(1/2)*tau-1))*pow(t,(-(1/2)*tau-3))*M_PI*
    (
     -sqrt(2)*
     LaguerreL(-(1/2)*tau, 1/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     pow(lambda,(5/2))*sqrt(t)*cos((1/2)*M_PI*tau)*gamma(-(1/2)*tau+2)*
     pow(alpha,2)*pow(kappa,3)+sqrt(2)*
     LaguerreL(-(1/2)*tau, 3/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))* pow(lambda,5/2)*sqrt(t)*cos((1/2)*M_PI*tau)*gamma(-(1/2)*tau+2)*pow(alpha,2)*
     
     pow(kappa,3)-sqrt(2)*
     LaguerreL(-(1/2)*tau, 1/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     pow(lambda,3/2)*pow(t,3/2)*cos((1/2)*M_PI*tau)*gamma(-(1/2)*tau+2)*pow(kappa,3)+
     2*sqrt(2)*
     LaguerreL(-(1/2)*tau, 1/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     pow(lambda,3/2)*sqrt(t)*cos((1/2)*M_PI*tau)*gamma(-(1/2)*tau+2)*alpha*
     pow(kappa,2)-2*sqrt(2)*
     LaguerreL(-(1/2)*tau, 3/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     pow(lambda,3/2)*sqrt(t)*cos((1/2)*M_PI*tau)*gamma(-(1/2)*tau+2)*alpha*
     pow(kappa,2)-sqrt(2)*
     LaguerreL(-(1/2)*tau, 1/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     sqrt(lambda)*sqrt(t)*cos((1/2)*M_PI*tau)*gamma(-(1/2)*tau+2)*kappa+sqrt(2)*LaguerreL(-(1/2)*tau, 3/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     sqrt(lambda)*sqrt(t)*cos((1/2)*M_PI*tau)*gamma(-(1/2)*tau+2)*kappa+
     sin((1/2)*M_PI*tau)*gamma(-(1/2)*tau+3/2)*
     LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     pow(alpha,3)*pow(kappa,3)*pow(lambda,3)-sin((1/2)*M_PI*tau)*gamma(-(1/2)*tau+3/2)*
     LaguerreL(-(1/2)*tau+1/2, 3/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     pow(alpha,3)*pow(kappa,3)*pow(lambda,3)+sin((1/2)*M_PI*tau)*gamma(-(1/2)*tau+3/2)*
     LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))* alpha*pow(kappa,3)*pow(lambda,2)*tau*t-sin((1/2)*M_PI*tau)*gamma(-(1/2)*tau+3/2)*
     LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     alpha*pow(kappa,3)*pow(lambda,2)*t-3*sin((1/2)*M_PI*tau)*gamma(-(1/2)*tau+3/2)*
     
     LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     pow(alpha,2)*pow(kappa,2)*pow(lambda,2)+3*sin((1/2)*M_PI*tau)*gamma(-(1/2)*tau+
                                                                         3/2)*
     
     LaguerreL(-(1/2)*tau+1/2, 3/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     pow(alpha,2)*pow(kappa,2)*pow(lambda,2)-sin((1/2)*M_PI*tau)*gamma(-(1/2)*tau+3/2)*
     LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     pow(kappa,2)*lambda*tau*t+sin((1/2)*M_PI*tau)*gamma(-(1/2)*tau+3/2)*
     LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     pow(kappa,2)*lambda*t+3*sin((1/2)*M_PI*tau)*gamma(-(1/2)*tau+3/2)*
     LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     alpha*kappa*lambda-3*sin((1/2)*M_PI*tau)*gamma(-(1/2)*tau+3/2)*
     LaguerreL(-(1/2)*tau+1/2, 3/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     alpha*kappa*lambda-
     LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     gamma(-(1/2)*tau+3/2)*sin((1/2)*M_PI*tau)+
     LaguerreL(-(1/2)*tau+1/2, 3/2, (1/2)*pow(alpha*kappa*lambda-1,2)/(lambda*pow(kappa,2)*t))*
     gamma(-(1/2)*tau+3/2)*sin((1/2)*M_PI*tau)
     )/
    (gamma(tau)*sin((1/2)*M_PI*tau)*gamma(-(1/2)*tau+3/2)*cos((1/2)*M_PI*tau)*gamma(-(1/2)*tau+2));
    */
  }   
  return d;
}    
     
SEXP  dwald(SEXP t, SEXP lambda, SEXP alpha, SEXP tau, SEXP kappa) {
  double d;
  SEXP value;
    
  d = dwald_d(REAL(t)[0], REAL(lambda)[0], REAL(alpha)[0], REAL(tau)[0], REAL(kappa)[0]);
    
  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = d;
  UNPROTECT(1);
  return value;
}    
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
