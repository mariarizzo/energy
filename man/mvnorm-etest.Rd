\name{mvnorm.test}
\alias{mvnorm.test}
\alias{mvnorm.etest}
\alias{mvnorm.e}
\title{E-statistic (Energy) Test of Multivariate Normality}
\description{
 Performs the E-statistic (energy) test of multivariate or univariate normality.
}
\usage{
mvnorm.test(x, R)
mvnorm.etest(x, R)
mvnorm.e(x)
}
\arguments{
  \item{x}{ data matrix of multivariate sample, or univariate data vector}
  \item{R}{ number of bootstrap replicates }
}
\details{
 If \code{x} is a matrix, each row is a multivariate observation. The
 data will be standardized to zero mean and identity covariance matrix
 using the sample mean vector and sample covariance matrix. If \code{x}
 is a vector, \code{mvnorm.e} returns the univariate statistic 
 \code{normal.e(x)}.
 If the data contains missing values or the sample covariance matrix is
 singular, \code{mvnorm.e} returns NA. 

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
}
\value{
 The value of the \eqn{\mathcal{E}}{E}-statistic for multivariate
 normality is returned by \code{mvnorm.e}.

 \code{mvnorm.test} returns a list with class \code{htest} containing
 \item{method}{description of test}
 \item{statistic}{observed value of the test statistic}
 \item{p.value}{approximate p-value of the test}
 \item{data.name}{description of data}
 
 \code{mvnorm.etest} is replaced by \code{mvnorm.test}. 
}
\seealso{
\code{\link{normal.test}} for the energy test of univariate
normality and \code{\link{normal.e}} for the statistic.
}
\note{If the data is univariate, the test statistic is formally
the same as the multivariate case, but a more efficient computational
formula is applied in \code{\link{normal.e}}.

\code{\link{normal.test}} also provides an optional method for the
test based on the asymptotic sampling distribution of the test
statistic.
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
 ## compute normality test statistic for iris Setosa data
 data(iris)
 mvnorm.e(iris[1:50, 1:4])

 ## test if the iris Setosa data has multivariate normal distribution
 mvnorm.etest(iris[1:50,1:4], R = 199)
}
\keyword{ multivariate }
\keyword{ htest }
\concept{ goodness-of-fit}
\concept{ normal distribution}
\concept{ energy statistics }
