rswald <- function(n, alpha, nu, theta, use.rwiener=TRUE) {
  beta.par <- 0.5
  if(use.rwiener)
    x <- rwiener(n, alpha/(1-beta.par), theta, beta.par, nu)
  else
    lambda <- 1
    x <- theta + rinvgauss(n, mean=alpha/nu, shape=lambda*alpha^2)
  return(x)
}

qswald <- function(p, alpha, nu, theta) {

  if (abs(p) > 1) return(NaN)

  pmin <- 0
  pmax <- 1
  pmid <- 0
  qmin <- 0
  qmax <- Inf
  q <- 1

  c <- 0
  repeat {
    c <- c+1
    if (p>=0) pmid = pswald.r(q, alpha, nu, theta);
    if (abs(p)<=pmid) { # near lower point
      pmax <- pmid
      qmax <- q
      q <- qmin + (qmax-qmin)/2
    }
    else { # near upper point
      pmin <- pmid
      qmin <- q
      if (is.finite(qmax)) 
        q <- qmin + (qmax-qmin)/2
      else
        q <- q*10
    }
    if(is.nan(pmid)) return(NaN)
    if(q>=1e+10) return(NaN)
    if (!(abs(abs(p)-pmid) > 1e-10 && c < 1000)) 
      break
  } 

  return(q)
}
