\name{dcorT}
\alias{dcorT.test}
\alias{dcorT}
\title{ Distance Correlation t-Test}
\description{
 Distance correlation t-test of multivariate independence for high dimension.}
\usage{
dcorT.test(x, y)
dcorT(x, y)
}
\arguments{
  \item{x}{ data or distances of first sample}
  \item{y}{ data or distances of second sample}
}
\details{
 \code{dcorT.test} performs a nonparametric t-test of
 multivariate independence in high dimension (dimension is close to
 or larger than sample size). As dimension goes to infinity, the
 asymptotic distribution of the test statistic is approximately Student t with \eqn{n(n-3)/2-1} degrees of freedom and for \eqn{n \geq 10} the statistic is approximately distributed as standard normal.

 The sample sizes (number of rows) of the two samples must
 agree, and samples must not contain missing values.

 The t statistic (dcorT) is a transformation of a bias corrected 
 version of distance correlation (see SR 2013 for details).

Large values (upper tail) of the dcorT statistic are significant.
}
\note{
\code{dcor.t} and \code{dcor.ttest} are deprecated.
}
\value{
\code{dcorT} returns the dcor t statistic, and
\code{dcorT.test} returns a list with class \code{htest} containing
   \item{     method}{ description of test}
   \item{  statistic}{ observed value of the test statistic}
   \item{  parameter}{ degrees of freedom}
   \item{   estimate}{ (bias corrected) squared dCor(x,y)}
   \item{    p.value}{ p-value of the t-test}
   \item{  data.name}{ description of data}
}
\seealso{
 \code{\link{bcdcor}} \code{\link{dcov.test}} \code{\link{dcor}} \code{\link{DCOR}}
}

\references{
 Szekely, G.J. and Rizzo, M.L. (2013). The distance correlation t-test of  independence in high dimension. \emph{Journal of Multivariate Analysis},  Volume 117, pp. 193-213. \cr
 \url{http://dx.doi.org/10.1016/j.jmva.2013.02.012}

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
\author{
Maria L. Rizzo \email{mrizzo @ bgsu.edu} and
Gabor J. Szekely
}
\examples{
 x <- matrix(rnorm(100), 10, 10)
 y <- matrix(runif(100), 10, 10)
 dcorT(x, y)
 dcorT.test(x, y)
}


\keyword{ htest }
\keyword{ multivariate }
\keyword{ nonparametric }
\concept{ independence }
\concept{ multivariate }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }

