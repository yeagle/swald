phi <- function(x) {
  c <- 0
  for (i in 1:x) {
    c <- c + 1/i
  }
  return(c)
}

neg_gamma <- function(x) -pi / (-x*gamma(-x)*sin(pi*(-x)))

neg_int_gamma <- function(x) ((-1) ^ (-x) / factorial(-x)) * (phi(-x) + digamma(1))

LaguerreL <- function(n, a, x) {
  if (n+a+1 <= 0) {
    if (n+a+1 == 0) {
      c1 <- digamma(1)  # - Euler's constant
    } else {
      if ((n+a+1)%%1==0) {
        c1 <- neg_int_gamma(n+a+1)	
      } else {	
        c1 <- neg_gamma(n+a+1)	    	
      }
    }
  }else {
    c1 <- gamma(n+a+1)	    	
  }
    
  if (n+1 <= 0) {
    if (n+1 == 0) {
      c2 <- digamma(1)  # - Euler's constant
    } else {
      if ((n+1)%%1==0) {
        c2 <- neg_int_gamma(n+1)	
      } else {	
        c2 <- neg_gamma(n+1)	    	
      }
    }
  }else {
    c2 <- gamma(n+1)	    	
  }
  
  if (a+1 <= 0) {
    if (a+1 == 0) {
      c3 <- digamma(1)  # - Euler's constant
    } else {
      if ((a+1)%%1==0) {
        c3 <- neg_int_gamma(a+1)	
      } else {	
        c3 <- neg_gamma(a+1)	    	
      }
    }
  }else {
    c3 <- gamma(a+1)	    	
  }

  L <- c1/(c2*c3) * genhypergeo(U=-n, L=a+1, z=x)
  return(L)
}

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

dwald_r <- function(t, alpha, tau, kappa){
  dwald_gamma_r(t, alpha, tau, kappa)
}

dwald_gamma_r <- function(t, alpha, tau, kappa, give_log=FALSE){
  # In
  #   t: RT [ms] >0	
  #   alpha: boundary separation >0  
  #   tau: scale parameter of the gamma distribution >0
  #   kappa: shape parameter of the gamma distribution >0	
  if (give_log) {
    d <- dwald_gamma_r_log(t, alpha, tau, kappa)
  } else {  
    if (alpha == 1 && tau == 1 && kappa == 1 )
    {
      d = exp(- 1 / (2*t)) / (2 * t^2)
    }
    else if (tau == 1) 
    {
      d = alpha*exp(-(2*alpha*kappa-1)/(2*t*kappa^2))* 
          (1 + erf( (alpha*kappa-1) / (kappa*sqrt(2*t)) )) / (2*(t^2)*kappa)
    }
    else 
    {
      L1 <- LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*(alpha*kappa-1)^2/(kappa^2*t))
    
      L2 <- LaguerreL(-(1/2)*tau+1/2, 3/2, (1/2)*(alpha*kappa-1)^2/(kappa^2*t))
    
      L3 <- LaguerreL(-(1/2)*tau, 1/2, (1/2)*(alpha*kappa-1)^2/(kappa^2*t))
    
      L4 <- LaguerreL(-(1/2)*tau, 3/2, (1/2)*(alpha*kappa-1)^2/(kappa^2*t))
    
      C1 <- sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)
    
      C2 <- sqrt(2)*alpha^3*kappa^3*sqrt(t)
    
      d = -(1/16)*2^((1/2)*tau+1/2)*alpha*exp(-(1/2)*alpha^2/t)*kappa^(-tau-3)*
          t^(-(1/2)*tau-7/2)*pi*
    
      (
    
      C1*
      L1*
      C2 + 
    
      C1*
      L1*
      sqrt(2)*alpha*kappa^3*tau*t^(3/2) -
    
      C1*
      L2*
      C2 - 
    
      C1*
      L1*
      sqrt(2)*alpha*kappa^3*t^(3/2) - 
    
      3*C1*
      L1*
      sqrt(2)*alpha^2*kappa^2*sqrt(t) - 
    
      C1*
      L1*
      sqrt(2)*kappa^2*tau*t^(3/2) +
    
      3*C1*
      L2*
      sqrt(2)*alpha^2*kappa^2*sqrt(t)-
    
      2*
      L3*
      cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*alpha^2*kappa^3*t +
    
      2*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*
      L4*
      alpha^2*kappa^3*t +
    
      C1*
      L1*
      sqrt(2)*kappa^2*t^(3/2) -
    
      2*
      L3*
      cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*kappa^3*t^2 +
    
      3*C1*
      L1*
      sqrt(2)*alpha*kappa*sqrt(t) - 
    
      3*C1*
      L2*
      sqrt(2)*alpha*kappa*sqrt(t) + 
    
      4* 			
      L3*
      cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*alpha*kappa^2*t -
    
      4*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*
      L4*
      alpha*kappa^2*t -
    
      C1*
      L1*
      sqrt(2)*sqrt(t) +
    
      C1*
      L2*
      sqrt(2)*sqrt(t) -
    
      2*
      L3*
      cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*kappa*t +
    
      2*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*
      L4*
      kappa*t
    
      )/
    
      (gamma(tau)*C1*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2))
   }
 }
  return(d)
}

dwald_gamma_r_log <- function(t, alpha, tau, kappa){
  # In
  #   t: RT [ms] >0	
  #   alpha: boundary separation >0  
  #   tau: scale parameter of the gamma distribution >0
  #   kappa: shape parameter of the gamma distribution >0	
  if (alpha == 1 && tau == 1 && kappa == 1 )
  {
    d = (- 1 / (2*t)) - log(2) - 2*log(t)
  }
  else if (tau == 1) 
  {
    d = log(alpha) - ((2*alpha*kappa-1)/(2*t*kappa^2)) +  
        log(1 + erf( (alpha*kappa-1) / (kappa*sqrt(2*t)) )) - log(2) - 2*log(t) - log(kappa)
  }
  else 
  {
    L1 <- LaguerreL(-(1/2)*tau+1/2, 1/2, (1/2)*(alpha*kappa-1)^2/(kappa^2*t))
    
    L2 <- LaguerreL(-(1/2)*tau+1/2, 3/2, (1/2)*(alpha*kappa-1)^2/(kappa^2*t))
    
    L3 <- LaguerreL(-(1/2)*tau, 1/2, (1/2)*(alpha*kappa-1)^2/(kappa^2*t))
    
    L4 <- LaguerreL(-(1/2)*tau, 3/2, (1/2)*(alpha*kappa-1)^2/(kappa^2*t))
    
    C1 <- sin((1/2)*pi*tau)*gamma(-(1/2)*tau+3/2)
    
    C2 <- sqrt(2)*alpha^3*kappa^3*sqrt(t)
    
    d = log(-(1/16)*2^((1/2)*tau+1/2)*alpha*exp(-(1/2)*alpha^2/t)*kappa^(-tau-3)*
        t^(-(1/2)*tau-7/2)*pi *
    
    (
    
    C1*
    L1*
    C2 
    
    + 
    
    C1*
    L1*
    sqrt(2)*alpha*kappa^3*tau*t^(3/2) -
    
    C1*
    L2*
    C2 - 
    
    C1*
    L1*
    sqrt(2)*alpha*kappa^3*t^(3/2) - 
    
    3*C1*
    L1*
    sqrt(2)*alpha^2*kappa^2*sqrt(t) - 
    
    C1*
    L1*
    sqrt(2)*kappa^2*tau*t^(3/2) +
    
    3*C1*
    L2*
    sqrt(2)*alpha^2*kappa^2*sqrt(t)-
    
    2*
    L3*
    cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*alpha^2*kappa^3*t +
    
    2*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*
    L4*
    alpha^2*kappa^3*t +
    
    C1*
    L1*
    sqrt(2)*kappa^2*t^(3/2) -
    
    2*
    L3*
    cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*kappa^3*t^2 +
    
    3*C1*
    L1*
    sqrt(2)*alpha*kappa*sqrt(t) - 
    
    3*C1*
    L2*
    sqrt(2)*alpha*kappa*sqrt(t) + 
    
    4* 			
    L3*
    cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*alpha*kappa^2*t -
    
    4*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*
    L4*
    alpha*kappa^2*t -
    
    C1*
    L1*
    sqrt(2)*sqrt(t) +
    
    C1*
    L2*
    sqrt(2)*sqrt(t) -
    
    2*
    L3*
    cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*kappa*t +
    
    2*cos((1/2)*pi*tau)*gamma(-(1/2)*tau+2)*
    L4*
    kappa*t

       
    ) ) -
    
    log(gamma(tau)) - log(C1) - log(cos((1/2)*pi*tau)) - log(gamma(-(1/2)*tau+2))
 }
  return(d)
}

dIG_r <- function(t, alpha, lambda, nu) {
  alpha * (lambda / (2*pi*t^3)) ^ .5 * exp( -lambda * (nu*t-alpha)^2 / (2*t))
}

dswald_r <- function(t, alpha, gamma, theta) {
  alpha * (2*pi*((t-theta)^3)) ^ (-.5) * exp( -(alpha-gamma*(t-theta))^2 / (2*(t-theta)))
}
