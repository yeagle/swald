library(hypergeo)

LaguerreL <- function(n, a, x) {
  ((n+a+1)*beta(1+a,n+1))^-1 * genhypergeo(U=-n, L=a+1, z=x)
}

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

Maple = function(lambda, alpha, tau, t, kappa){
  # In
  #   t: RT [ms] >0	
  #   alpha: boundary separation >0  
  #   lambda: diffusion coefficient >0
  #   tau: scale parameter of the gamma distribution >0
  #   kappa: shape parameter of the gamma distribution >0	
  if (tau==1) {
  	alpha*exp(-(2*alpha*kappa-1)/(2*t*kappa^2))* (1 + erf( (alpha*kappa-1) / (kappa*sqrt(2*t)) )) / (2*t^2*kappa)
  }else
  {
-(1/16)*sqrt(2)*2^((1/2)*tau+1/2)*alpha*exp(-(1/2)*alpha^2*lambda/t)*kappa^(-tau-3)*lambda^(-(1/2)*tau-1)*t^(-(1/2)*tau-3)*pi*(-sqrt(2)*LaguerreL(-(1/2)*tau, 1/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*lambda^(5/2)*sqrt(t)*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*alpha^2*kappa^3+sqrt(2)*LaguerreL(-(1/2)*tau, 3/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*lambda^(5/2)*sqrt(t)*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*alpha^2*kappa^3-sqrt(2)*LaguerreL(-(1/2)*tau, 1/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*lambda^(3/2)*t^(3/2)*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*kappa^3+2*sqrt(2)*LaguerreL(-(1/2)*tau, 1/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*lambda^(3/2)*sqrt(t)*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*alpha*kappa^2-2*sqrt(2)*LaguerreL(-(1/2)*tau, 3/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*lambda^(3/2)*sqrt(t)*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*alpha*kappa^2-sqrt(2)*LaguerreL(-(1/2)*tau, 1/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*sqrt(lambda)*sqrt(t)*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*kappa+sqrt(2)*LaguerreL(-(1/2)*tau, 3/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*sqrt(lambda)*sqrt(t)*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*kappa+sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)*LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*alpha^3*kappa^3*lambda^3-sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)*LaguerreL(-(1/2)*tau+1/2, 3/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*alpha^3*kappa^3*lambda^3+sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)*LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*alpha*kappa^3*lambda^2*tau*t-sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)*LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*alpha*kappa^3*lambda^2*t-3*sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)*LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*alpha^2*kappa^2*lambda^2+3*sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)*LaguerreL(-(1/2)*tau+1/2, 3/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*alpha^2*kappa^2*lambda^2-sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)*LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*kappa^2*lambda*tau*t+sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)*LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*kappa^2*lambda*t+3*sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)*LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*alpha*kappa*lambda-3*sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)*LaguerreL(-(1/2)*tau+1/2, 3/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*alpha*kappa*lambda-LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*gamma(-(1/2)*tau+3/2)*sin((1/2)*pi*tau)+LaguerreL(-(1/2)*tau+1/2, 3/2, (1/2)*(alpha*kappa*lambda-1)^2/(lambda*kappa^2*t))*gamma(-(1/2)*tau+3/2)*sin((1/2)*pi*tau))/(gamma(tau)*sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2))
 }
}
