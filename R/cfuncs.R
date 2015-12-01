dswald <- function(t, alpha, nu, theta, give_log=FALSE){
  d <- vector("double", length=length(t))
  for (i in 1:length(t)) {
    d[i] <- .Call(dswald_c, t[i], alpha, nu, theta, give_log)
  }
  return(d)
}

pswald <- function(t, alpha, nu, theta, lower.tail=TRUE, log.p=FALSE){
  d <- vector("double", length=length(t))
  for (i in 1:length(t)) {
    d[i] <- .Call(pswald_c, t[i], alpha, nu, theta, lower.tail, log.p)
  }
  return(d)
}
