\name{swald}
\alias{dswald}
\alias{pswald}
\alias{qswald}
\alias{rswald}
\title{Shifted Wald distribution functions}
\description{
  Simple reaction time distribution functions,
  based on a shifted Wald (inverse Gaussian) distribution with density,
  CDF, quantile and random functions:
  \code{dswald} computes the shifted wald distribution density.
  \code{pswald} computes the CDF for the shifted wald distribution.
  \code{qswald} computes the quantile for a given CDF value.
  \code{rswald} generates random quantiles from a wiener process
  distribution (if use.rwiener is TRUE) or from a inverse gaussian
  distribution (if use.rwiener is FALSE).
}
\usage{
  dswald(t, alpha,nu,theta, give_log=FALSE)
  pswald(t, alpha,nu,theta, lower.tail=TRUE, log.p=FALSE)
  qswald(p, alpha,nu,theta)
  rswald(n, alpha,nu,theta, use.rwiener=TRUE)

}
\arguments{
  \item{t}{vector of reaction time quantiles.}
  \item{p}{vector of probabilities.}  
  \item{n}{number of observations.}  
  \item{alpha}{1st model parameter.}
  \item{nu}{2nd model parameter.}
  \item{theta}{non-decision time parameter.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].}
  \item{give_log}{logical; function returns log, if this argument is TRUE}
  \item{log.p}{logical; if TRUE probabilities are given as log(p)}
  \item{use.rwiener}{function returns log, if this argument is TRUE}
}

\examples{
dswald(.5,1,1,.1)
dswald.r(.5,1,1,.1)
}

\keyword{dwiener}
\keyword{pwiener}
\keyword{qwiener}
\keyword{rwiener}

\references{
R. Miller (todo)
}
