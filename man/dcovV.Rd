\name{distance covariance test}
\alias{dcovV.test}
\alias{dcovU.test}
\title{ Distance Covariance Test by Estimating Limit Distribution }
\description{
 Computes distance covariance and the p-value for a test of multivariate independence by estimating the eigenvalues and numerical integration of the limit distribution.
}
\usage{
dcovV.test(x, y, max.eigen = 250)
dcovU.test(x, y, max.eigen = 250)
}
\arguments{
  \item{x}{ data matrix of first sample}
  \item{y}{ data matrix of second sample}
  \item{max.eigen}{ integer, max number of eigenvalues}
}
\details{
This test requires the CompQuadForm package; if not present it will return NA for the p-value.

This large sample dcov test applies the limit distribution of of the V-statistic, which is a quadratic form
\deqn{\sum_{i=1}^\infty \lambda_i Z_i^2}
where \eqn{\lambda_i}{lambda} are positive constants that depend on the distributions and \eqn{Z_i}{Z_i} are iid standard normal random variables. Eigenvalues are estimated from the sample and numerical integration (Imhof's method) is used to compute the upper tail probability or p-value. For large samples, this should be faster than the permutation test. 

The sample sizes (number of rows) of the two samples must
agree, and samples must not contain missing values. 

See \code{\link{dcov}} for details.
See \code{\link{dcov.test}} for the permutation test 
based on the distance covariance statistic.
}
\value{
\code{dcovV.test} and \code{dcovU.test} return a list of class htest. The eigenvalues are returned
in the component \code{estimates}.
}
\note{
This is a large sample test. It is not recommended for sample size less than 100. For small samples \code{\link{dcov.test}} should be applied. If \eqn{d < 50}{d < 50} the call will be passed through to \code{\link{dcov.test}}. 
}
\seealso{
\code{\link{bcdcor}}  \code{\link{dcovU}}  \code{\link{pdcor}}
\code{\link{dcov.test}} \code{\link{dcor.ttest}}  \code{\link{pdcor.test}}
}
\references{
 Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007),
 Measuring and Testing Dependence by Correlation of Distances,
 \emph{Annals of Statistics}, Vol. 35 No. 6, pp. 2769-2794.
 \cr \url{http://dx.doi.org/10.1214/009053607000000505}

 Szekely, G.J. and Rizzo, M.L. (2009),
 Brownian Distance Covariance,
 \emph{Annals of Applied Statistics},
 Vol. 3, No. 4, 1236-1265.
 \cr \url{http://dx.doi.org/10.1214/09-AOAS312}
}
\author{ Maria L. Rizzo \email{mrizzo @ bgsu.edu} and
Gabor J. Szekely
}
\examples{
  dat <- MASS::Pima.tr
  x <- dat[["glu"]]
  y <- dat[["ped"]]
  dcov.test(x, y, R=499)
  dcovV.test(x, y)
  dcovU.test(x, y)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ distance covariance }
\concept{ energy statistics }

