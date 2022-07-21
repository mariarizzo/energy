\name{normal.test}
\alias{normal.test}
\alias{normal.e}
\title{Energy Test of Univariate Normality}
\description{
 Performs the energy test of univariate normality
 for the composite hypothesis Case 4, estimated parameters.
}
\usage{
normal.test(x, method=c("mc","limit"), R)
normal.e(x)
}
\arguments{
  \item{x}{ univariate data vector}
  \item{method}{ method for p-value}
  \item{R}{ number of replications if Monte Carlo method}
}
\details{
If \code{method="mc"} this test function applies the parametric
bootstrap method implemented in \code{\link{mvnorm.test}}. 

If \code{method="limit"}, the p-value of the test is computed from 
the asymptotic distribution of the test statistic under the null 
hypothesis. The asymptotic 
distribution is a quadratic form of centered Gaussian random variables, 
which has the form
\deqn{\sum_{k=1}^\infty \lambda_k Z_k^2,}
where \eqn{\lambda_k} are positive constants (eigenvalues) and
\eqn{Z_k} are iid standard normal variables. Eigenvalues are
pre-computed and stored internally. 
A p-value is computed using Imhof's method as implemented in the
\pkg{CompQuadForm} package. 

Note that the "limit" method is intended for moderately large
samples because it applies the asymptotic distribution.  

The energy test of normality was proposed
 and implemented by Szekely and Rizzo (2005). 
 See \code{\link{mvnorm.test}}
 for more details. 
}

\value{
\code{normal.e} returns the energy goodness-of-fit statistic for
a univariate sample. 

 \code{normal.test} returns a list with class \code{htest} containing
 \item{statistic}{observed value of the test statistic}
 \item{p.value}{p-value of the test}
 \item{estimate}{sample estimates: mean, sd}
 \item{data.name}{description of data}
}
\seealso{
  \code{\link{mvnorm.test}} and \code{\link{mvnorm.e}} for the
  energy test of multivariate normality and the test statistic
  for multivariate samples.
}
\references{
 Szekely, G. J. and Rizzo, M. L. (2005) A New Test for
 Multivariate Normality, \emph{Journal of Multivariate Analysis},
 93/1, 58-80,
 \doi{10.1016/j.jmva.2003.12.002}.

Mori, T. F., Szekely, G. J. and Rizzo, M. L. "On energy tests of normality." Journal of Statistical Planning and Inference 213 (2021): 1-15. 

 Rizzo, M. L. (2002). A New Rotation Invariant Goodness-of-Fit Test,
 Ph.D. dissertation, Bowling Green State University.

J. P. Imhof (1961). Computing the Distribution of Quadratic Forms in 
Normal Variables, \emph{Biometrika}, Volume 48, Issue 3/4,
419-426.
}
\author{ Maria L. Rizzo \email{mrizzo@bgsu.edu} and
Gabor J. Szekely
}
\examples{
  x <- iris[1:50, 1]
  normal.e(x)
  normal.test(x, R=199)
  normal.test(x, method="limit")
}
\keyword{ htest }
\concept{ goodness-of-fit}
\concept{ normal distribution}
\concept{ energy statistics }
