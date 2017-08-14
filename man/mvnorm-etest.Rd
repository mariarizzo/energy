\name{mvnorm.etest}
\alias{mvnorm.etest}
\alias{mvnorm.e}
\alias{normal.e}
\title{E-statistic (Energy) Test of Multivariate Normality}
\description{
 Performs the E-statistic (energy) test of multivariate or univariate normality.
}
\usage{
 mvnorm.etest(x, R)
 mvnorm.e(x)
 normal.e(x)
}
\arguments{
  \item{x}{ data matrix of multivariate sample, or univariate data vector}
  \item{R}{ number of bootstrap replicates }
}
\details{
 If \code{x} is a matrix, each row is a multivariate observation. The
 data will be standardized to zero mean and identity covariance matrix
 using the sample mean vector and sample covariance matrix. If \code{x}
 is a vector, the univariate statistic \code{normal.e(x)} is returned.
 If the data contains missing values or the sample covariance matrix is
 singular, NA is returned.

 The \eqn{\mathcal{E}}{E}-test of multivariate normality was proposed
 and implemented by Szekely and Rizzo (2005). The test statistic for
 d-variate normality is given by
 \deqn{\mathcal{E} = n (\frac{2}{n} \sum_{i=1}^n E\|y_i-Z\| -
 E\|Z-Z'\| - \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n \|y_i-y_j\|),
 }{E = n((2/n) sum[1:n] E||y_i-Z|| - E||Z-Z'|| - (1/n^2) sum[1:n,1:n]
 ||y_i-y_j||),}
 where \eqn{y_1,\ldots,y_n} is the standardized sample,
 \eqn{Z, Z'} are iid standard d-variate normal, and
 \eqn{\| \cdot \|}{|| ||} denotes Euclidean norm.

   The \eqn{\mathcal{E}}{E}-test of multivariate (univariate) normality
  is implemented by parametric bootstrap with \code{R} replicates.
  
  If \code{R=0} the summary for the test gives the test statistic only (no p-value). 
}
\value{
 The value of the \eqn{\mathcal{E}}{E}-statistic for univariate
 normality is returned by \code{normal.e}.
 The value of the \eqn{\mathcal{E}}{E}-statistic for multivariate
 normality is returned by \code{mvnorm.e}.

 \code{mvnorm.etest} returns a list with class \code{htest} containing
 \item{method}{description of test}
 \item{statistic}{observed value of the test statistic}
 \item{p.value}{approximate p-value of the test}
 \item{data.name}{description of data}
}

\references{
 Szekely, G. J. and Rizzo, M. L. (2005) A New Test for
 Multivariate Normality, \emph{Journal of Multivariate Analysis},
 93/1, 58-80,
 \url{http://dx.doi.org/10.1016/j.jmva.2003.12.002}.

 Rizzo, M. L. (2002). A New Rotation Invariant Goodness-of-Fit Test,
 Ph.D. dissertation, Bowling Green State University.

 Szekely, G. J. (1989) Potential and Kinetic Energy in Statistics,
 Lecture Notes, Budapest Institute of Technology (Technical University).
}
\author{ Maria L. Rizzo \email{mrizzo @ bgsu.edu} and
Gabor J. Szekely
}
\examples{
 ## compute normality test statistics for iris Setosa data
 data(iris)
 mvnorm.e(iris[1:50, 1:4])
 normal.e(iris[1:50, 1])

 ## test if the iris Setosa data has multivariate normal distribution
 mvnorm.etest(iris[1:50,1:4], R = 199)

 ## test a univariate sample for normality
 x <- runif(50, 0, 10)
 mvnorm.etest(x, R = 199)
 }
\keyword{ multivariate }
\keyword{ htest }
\concept{ energy statistics }
