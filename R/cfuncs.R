dwald <- function(t, alpha, tau, kappa){
  dwald_gamma(t, alpha, tau, kappa)
}

dIG <- function(t, alpha, lambda, nu){
  d <- vector("double", length=length(t))
  for (i in 1:length(t)) {
      d[i] <- .Call(dIG_c, t[i], alpha, lambda, nu)
  }
  return(d)
}

dswald <- function(t, alpha, gamma, theta){
  d <- vector("double", length=length(t))
  for (i in 1:length(t)) {
    d[i] <- .Call(dswald_c, t[i], alpha, gamma, theta)
  }
  return(d)
}

dwald_gamma <- function(t, alpha, tau, kappa){
  d <- vector("double", length=length(t))
  for (i in 1:length(t)) {
    d[i] <- .Call(dwald_gamma_c, t[i], alpha, tau, kappa)
  }
  return(d)
}


dwald_trunc <- function(t, lambda, alpha, v, d){
  w <- .Call(dwald_trunc_c, t, lambda, alpha, v, d)
  return(w)
} 

log_dwald_trunc <- function(t, lambda, alpha, v, d){
  w <- .Call(log_dwald_trunc_c, t, lambda, alpha, v, d)
  return(w)
} 
